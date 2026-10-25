// Block-domain KPM benchmark.
//
// Compares two ways of computing Chebyshev moments of a block-diagonal Hermitian
// operator (optionally with a sparse coupling between neighbouring blocks):
//   global: one CSR matrix of the whole system, full-length vectors, one sparse
//           product per Chebyshev step (the usual KPM layout);
//   local:  each block is a domain processed independently from start to end
//           (recursion, operator products and moment table in cache), then the
//           per-thread moment tables are summed. Exact only without coupling.
// Two observables: rank-1 moments mu_m = <r|T_m(H)|r> (density of states) and
// rank-2 moments mu_mn = <r|A T_m(H) B T_n(H)|r> (conductivity), where A and B
// are block-diagonal Hermitian operators.
//
// Usage: block_kpm <dos|cond> <global|local> <block> <dimension> <moments>
//                  <vectors> <coupling> <out_moments_file>
// Threads come from OMP_NUM_THREADS. Prints one JSON line.

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using cd = std::complex<double>;
using Clock = std::chrono::steady_clock;

namespace {

constexpr int kMaxBlock = 64;
constexpr int kMaxVectors = 8;

uint64_t Mix(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

double Uniform(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
  const uint64_t h = Mix(Mix(Mix(Mix(a) ^ b) ^ c) ^ d);
  return static_cast<double>(h >> 11) * 0x1.0p-53;
}

// Uninitialised array: pages are placed by the thread that first writes them (first touch).
template <class T>
class RawArray {
 public:
  RawArray() = default;
  explicit RawArray(size_t n) { resize(n); }
  RawArray(const RawArray&) = delete;
  RawArray& operator=(const RawArray&) = delete;
  RawArray(RawArray&& o) noexcept : p_(o.p_), n_(o.n_) { o.p_ = nullptr; o.n_ = 0; }
  ~RawArray() { std::free(p_); }
  void resize(size_t n) {
    std::free(p_);
    p_ = static_cast<T*>(std::malloc(n * sizeof(T)));
    n_ = n;
  }
  T* data() { return p_; }
  const T* data() const { return p_; }
  T& operator[](size_t i) { return p_[i]; }
  const T& operator[](size_t i) const { return p_[i]; }
  size_t size() const { return n_; }

 private:
  T* p_ = nullptr;
  size_t n_ = 0;
};

double Seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

struct Model {
  int64_t blocks = 0;
  int w = 0;
  int64_t n = 0;
  RawArray<cd> h, a, b;  // dense blocks, row-major, blocks * w * w
  // Global CSR of H (blocks plus coupling) and of A, B (blocks only).
  std::vector<int64_t> hPtr, aPtr;
  RawArray<int32_t> hCol, aCol;
  RawArray<cd> hVal, aVal, bVal;
  double scale = 1.0;
  int64_t couplingLinks = 0;
};

// Hermitian random block element (i, j) of operator `tag` in block `blk`.
cd BlockElement(int tag, int64_t blk, int i, int j, int w) {
  const int lo = std::min(i, j), hi = std::max(i, j);
  const double amp = 1.0 / std::sqrt(static_cast<double>(w));
  const double re = 2.0 * Uniform(tag, blk, lo, hi) - 1.0;
  if (i == j) return cd(re * amp, 0.0);
  const double im = 2.0 * Uniform(tag + 100, blk, lo, hi) - 1.0;
  const cd v(re * amp, im * amp);
  return i < j ? v : std::conj(v);
}

Model BuildModel(int w, int64_t dimension, bool coupled) {
  Model m;
  m.w = w;
  m.blocks = dimension / w;
  m.n = m.blocks * w;
  const int64_t ww = static_cast<int64_t>(w) * w;
  m.h.resize(m.blocks * ww);
  m.a.resize(m.blocks * ww);
  m.b.resize(m.blocks * ww);
#pragma omp parallel for schedule(static)
  for (int64_t blk = 0; blk < m.blocks; ++blk) {
    for (int i = 0; i < w; ++i) {
      for (int j = 0; j < w; ++j) {
        m.h[blk * ww + i * w + j] = BlockElement(1, blk, i, j, w);
        m.a[blk * ww + i * w + j] = BlockElement(2, blk, i, j, w);
        m.b[blk * ww + i * w + j] = BlockElement(3, blk, i, j, w);
      }
    }
  }

  // Coupling: every orbital of block k links to one random orbital of block k + 1.
  std::vector<int> extra(m.n, 0);
  std::vector<int32_t> partner(m.n, -1);
  std::vector<cd> link(m.n);
  if (coupled) {
    const double amp = 1.0 / std::sqrt(static_cast<double>(w));
    for (int64_t blk = 0; blk + 1 < m.blocks; ++blk) {
      for (int i = 0; i < w; ++i) {
        const int64_t row = blk * w + i;
        const int j = static_cast<int>(Uniform(7, blk, i, 0) * w) % w;
        const int64_t col = (blk + 1) * w + j;
        const double ph = 2.0 * M_PI * Uniform(8, blk, i, 0);
        partner[row] = static_cast<int32_t>(col);
        link[row] = amp * cd(std::cos(ph), std::sin(ph));
        ++extra[row];
        ++extra[col];
        ++m.couplingLinks;
      }
    }
  }

  m.hPtr.assign(m.n + 1, 0);
  m.aPtr.assign(m.n + 1, 0);
  for (int64_t r = 0; r < m.n; ++r) {
    m.hPtr[r + 1] = m.hPtr[r] + w + extra[r];
    m.aPtr[r + 1] = m.aPtr[r] + w;
  }
  m.hCol.resize(m.hPtr[m.n]);
  m.hVal.resize(m.hPtr[m.n]);
  m.aCol.resize(m.aPtr[m.n]);
  m.aVal.resize(m.aPtr[m.n]);
  m.bVal.resize(m.aPtr[m.n]);
  std::vector<int64_t> fill(m.n);
#pragma omp parallel for schedule(static)
  for (int64_t r = 0; r < m.n; ++r) {
    const int64_t blk = r / w;
    const int i = static_cast<int>(r % w);
    for (int j = 0; j < w; ++j) {
      const int64_t k = m.hPtr[r] + j;
      m.hCol[k] = static_cast<int32_t>(blk * w + j);
      m.hVal[k] = m.h[blk * ww + i * w + j];
      const int64_t ka = m.aPtr[r] + j;
      m.aCol[ka] = static_cast<int32_t>(blk * w + j);
      m.aVal[ka] = m.a[blk * ww + i * w + j];
      m.bVal[ka] = m.b[blk * ww + i * w + j];
    }
    fill[r] = m.hPtr[r] + w;
  }
  if (coupled) {
    for (int64_t r = 0; r < m.n; ++r) {
      if (partner[r] < 0) continue;
      const int64_t c = partner[r];
      m.hCol[fill[r]] = static_cast<int32_t>(c);
      m.hVal[fill[r]++] = link[r];
      m.hCol[fill[c]] = static_cast<int32_t>(r);
      m.hVal[fill[c]++] = std::conj(link[r]);
    }
  }

  // Gershgorin bound of H; the same scale is applied to the dense blocks.
  double bound = 0.0;
#pragma omp parallel for schedule(static) reduction(max : bound)
  for (int64_t r = 0; r < m.n; ++r) {
    double s = 0.0;
    for (int64_t k = m.hPtr[r]; k < m.hPtr[r + 1]; ++k) s += std::abs(m.hVal[k]);
    bound = std::max(bound, s);
  }
  m.scale = 1.0 / (1.01 * bound);
#pragma omp parallel for schedule(static)
  for (int64_t k = 0; k < static_cast<int64_t>(m.hVal.size()); ++k) m.hVal[k] *= m.scale;
#pragma omp parallel for schedule(static)
  for (int64_t k = 0; k < static_cast<int64_t>(m.h.size()); ++k) m.h[k] *= m.scale;
  return m;
}

cd RandomEntry(int vec, int64_t row) {
  const double ph = 2.0 * M_PI * Uniform(11, vec, row, 0);
  return cd(std::cos(ph), std::sin(ph));
}

// ---------------------------------------------------------------- global
// y[r*R+v] = alpha * (M x)[r*R+v] + beta * y[r*R+v]; returns sum_v <ref|y> when ref != nullptr.
template <int R>
cd CsrApply(const int64_t* ptr, const int32_t* col, const cd* val, int64_t n, const cd* x,
            cd* y, double alpha, double beta, const cd* ref) {
  double re = 0.0, im = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : re, im)
  for (int64_t r = 0; r < n; ++r) {
    cd acc[R];
    for (int v = 0; v < R; ++v) acc[v] = 0.0;
    for (int64_t k = ptr[r]; k < ptr[r + 1]; ++k) {
      const cd a = val[k];
      const cd* xc = x + static_cast<int64_t>(col[k]) * R;
      for (int v = 0; v < R; ++v) acc[v] += a * xc[v];
    }
    for (int v = 0; v < R; ++v) {
      const cd out = alpha * acc[v] + beta * y[r * R + v];
      y[r * R + v] = out;
      if (ref) {
        const cd d = std::conj(ref[r * R + v]) * out;
        re += d.real();
        im += d.imag();
      }
    }
  }
  return cd(re, im);
}

template <int R>
std::vector<cd> GlobalDos(const Model& m, int moments, double* tRec) {
  const int64_t len = m.n * R;
  RawArray<cd> r(len), x0(len), x1(len);
#pragma omp parallel for schedule(static)
  for (int64_t row = 0; row < m.n; ++row) {
    for (int v = 0; v < R; ++v) {
      r[row * R + v] = RandomEntry(v, row);
      x0[row * R + v] = r[row * R + v];
      x1[row * R + v] = 0.0;
    }
  }
  std::vector<cd> mu(moments);
  const auto t0 = Clock::now();
  double re = 0.0, im = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : re, im)
  for (int64_t k = 0; k < len; ++k) {
    const cd d = std::conj(r[k]) * x0[k];
    re += d.real();
    im += d.imag();
  }
  mu[0] = cd(re, im);
  if (moments > 1) {
    mu[1] = CsrApply<R>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), m.n, x0.data(), x1.data(),
                        1.0, 0.0, r.data());
  }
  cd* prev = x0.data();
  cd* cur = x1.data();
  for (int k = 2; k < moments; ++k) {
    mu[k] = CsrApply<R>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), m.n, cur, prev, 2.0, -1.0,
                        r.data());
    std::swap(prev, cur);
  }
  *tRec = Seconds(t0, Clock::now());
  return mu;
}

// Moment table mu[m][n] += sum_rows conj(L[m][row]) R[n][row]; vectors stored [moment][row].
void GlobalTable(const RawArray<cd>& left, const RawArray<cd>& right, int64_t n, int moments,
                 std::vector<cd>& mu) {
  const int threads = omp_get_max_threads();
  std::vector<std::vector<cd>> part(threads);
  constexpr int64_t chunk = 64;
  const int64_t chunks = (n + chunk - 1) / chunk;
#pragma omp parallel
  {
    std::vector<cd>& acc = part[omp_get_thread_num()];
    acc.assign(static_cast<size_t>(moments) * moments, 0.0);
    std::vector<cd> bl(chunk * moments), br(chunk * moments);
#pragma omp for schedule(static)
    for (int64_t c = 0; c < chunks; ++c) {
      const int64_t r0 = c * chunk;
      const int rows = static_cast<int>(std::min(chunk, n - r0));
      for (int mm = 0; mm < moments; ++mm) {
        const cd* lp = left.data() + static_cast<int64_t>(mm) * n + r0;
        const cd* rp = right.data() + static_cast<int64_t>(mm) * n + r0;
        for (int q = 0; q < rows; ++q) {
          bl[q * moments + mm] = std::conj(lp[q]);
          br[q * moments + mm] = rp[q];
        }
      }
      for (int q = 0; q < rows; ++q) {
        const cd* lq = bl.data() + q * moments;
        const cd* rq = br.data() + q * moments;
        for (int mm = 0; mm < moments; ++mm) {
          const cd l = lq[mm];
          cd* row = acc.data() + static_cast<int64_t>(mm) * moments;
          for (int nn = 0; nn < moments; ++nn) row[nn] += l * rq[nn];
        }
      }
    }
  }
  for (int t = 0; t < threads; ++t) {
    for (size_t k = 0; k < mu.size(); ++k) mu[k] += part[t][k];
  }
}

std::vector<cd> GlobalCond(const Model& m, int moments, int vectors, double* tRec, double* tTab) {
  const int64_t n = m.n;
  std::vector<cd> mu(static_cast<size_t>(moments) * moments, 0.0);
  RawArray<cd> left(static_cast<int64_t>(moments) * n), right(static_cast<int64_t>(moments) * n);
  RawArray<cd> r(n), tmp0(n), tmp1(n);
  *tRec = 0.0;
  *tTab = 0.0;
  for (int v = 0; v < vectors; ++v) {
#pragma omp parallel for schedule(static)
    for (int64_t row = 0; row < n; ++row) {
      r[row] = RandomEntry(v, row);
      for (int mm = 0; mm < moments; ++mm) {
        left[static_cast<int64_t>(mm) * n + row] = 0.0;
        right[static_cast<int64_t>(mm) * n + row] = 0.0;
      }
    }
    const auto t0 = Clock::now();
    // Left: L_m = T_m(H) A r.
    cd* l0 = left.data();
    CsrApply<1>(m.aPtr.data(), m.aCol.data(), m.aVal.data(), n, r.data(), l0, 1.0, 0.0, nullptr);
    if (moments > 1)
      CsrApply<1>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), n, l0, left.data() + n, 1.0, 0.0,
                  nullptr);
    for (int k = 2; k < moments; ++k) {
      cd* out = left.data() + static_cast<int64_t>(k) * n;
      std::memcpy(out, left.data() + static_cast<int64_t>(k - 2) * n, sizeof(cd) * n);
      CsrApply<1>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), n,
                  left.data() + static_cast<int64_t>(k - 1) * n, out, 2.0, -1.0, nullptr);
    }
    // Right: R_n = B T_n(H) r.
#pragma omp parallel for schedule(static)
    for (int64_t row = 0; row < n; ++row) {
      tmp0[row] = r[row];
      tmp1[row] = 0.0;
    }
    cd* prev = tmp0.data();
    cd* cur = tmp1.data();
    CsrApply<1>(m.aPtr.data(), m.aCol.data(), m.bVal.data(), n, prev, right.data(), 1.0, 0.0,
                nullptr);
    if (moments > 1) {
      CsrApply<1>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), n, prev, cur, 1.0, 0.0, nullptr);
      CsrApply<1>(m.aPtr.data(), m.aCol.data(), m.bVal.data(), n, cur, right.data() + n, 1.0, 0.0,
                  nullptr);
    }
    for (int k = 2; k < moments; ++k) {
      CsrApply<1>(m.hPtr.data(), m.hCol.data(), m.hVal.data(), n, cur, prev, 2.0, -1.0, nullptr);
      std::swap(prev, cur);
      CsrApply<1>(m.aPtr.data(), m.aCol.data(), m.bVal.data(), n, cur,
                  right.data() + static_cast<int64_t>(k) * n, 1.0, 0.0, nullptr);
    }
    const auto t1 = Clock::now();
    GlobalTable(left, right, n, moments, mu);
    const auto t2 = Clock::now();
    *tRec += Seconds(t0, t1);
    *tTab += Seconds(t1, t2);
  }
  return mu;
}

// ---------------------------------------------------------------- local
template <int R>
inline void BlockApply(const cd* blk, int w, const cd* x, cd* y) {
  for (int i = 0; i < w; ++i) {
    cd acc[R];
    for (int v = 0; v < R; ++v) acc[v] = 0.0;
    for (int j = 0; j < w; ++j) {
      const cd a = blk[i * w + j];
      for (int v = 0; v < R; ++v) acc[v] += a * x[j * R + v];
    }
    for (int v = 0; v < R; ++v) y[i * R + v] = acc[v];
  }
}

template <int R>
std::vector<cd> LocalDos(const Model& m, int moments, double* tRec) {
  constexpr int vectors = R;
  const int w = m.w;
  const int64_t ww = static_cast<int64_t>(w) * w;
  const int threads = omp_get_max_threads();
  std::vector<std::vector<cd>> part(threads);
  const auto t0 = Clock::now();
#pragma omp parallel
  {
    std::vector<cd>& acc = part[omp_get_thread_num()];
    acc.assign(moments, 0.0);
    cd r[kMaxBlock * kMaxVectors], x0[kMaxBlock * kMaxVectors], x1[kMaxBlock * kMaxVectors],
        y[kMaxBlock * kMaxVectors];
#pragma omp for schedule(static)
    for (int64_t blk = 0; blk < m.blocks; ++blk) {
      const cd* h = m.h.data() + blk * ww;
      const int len = w * vectors;
      for (int i = 0; i < w; ++i)
        for (int v = 0; v < vectors; ++v) r[i * vectors + v] = RandomEntry(v, blk * w + i);
      cd d = 0.0;
      for (int k = 0; k < len; ++k) {
        x0[k] = r[k];
        d += std::conj(r[k]) * x0[k];
      }
      acc[0] += d;
      if (moments > 1) {
        BlockApply<R>(h, w, x0, x1);
        d = 0.0;
        for (int k = 0; k < len; ++k) d += std::conj(r[k]) * x1[k];
        acc[1] += d;
      }
      cd* prev = x0;
      cd* cur = x1;
      for (int k = 2; k < moments; ++k) {
        BlockApply<R>(h, w, cur, y);
        d = 0.0;
        for (int q = 0; q < len; ++q) {
          prev[q] = 2.0 * y[q] - prev[q];
          d += std::conj(r[q]) * prev[q];
        }
        acc[k] += d;
        std::swap(prev, cur);
      }
    }
  }
  std::vector<cd> mu(moments, 0.0);
  for (int t = 0; t < threads; ++t)
    for (int k = 0; k < moments; ++k) mu[k] += part[t][k];
  *tRec = Seconds(t0, Clock::now());
  return mu;
}

std::vector<cd> LocalCond(const Model& m, int moments, int vectors, double* tAll) {
  const int w = m.w;
  const int64_t ww = static_cast<int64_t>(w) * w;
  const int threads = omp_get_max_threads();
  const size_t tableSize = static_cast<size_t>(moments) * moments;
  std::vector<std::vector<cd>> part(threads);
  const auto t0 = Clock::now();
#pragma omp parallel
  {
    std::vector<cd>& acc = part[omp_get_thread_num()];
    acc.assign(tableSize, 0.0);
    // Row-major per block: L[q][m] (conjugated) and R[q][n], q = orbital.
    std::vector<cd> lt(static_cast<size_t>(w) * moments), rt(static_cast<size_t>(w) * moments);
    cd r[kMaxBlock], p0[kMaxBlock], p1[kMaxBlock], y[kMaxBlock];
#pragma omp for schedule(static)
    for (int64_t blk = 0; blk < m.blocks; ++blk) {
      const cd* h = m.h.data() + blk * ww;
      const cd* a = m.a.data() + blk * ww;
      const cd* b = m.b.data() + blk * ww;
      for (int v = 0; v < vectors; ++v) {
        for (int i = 0; i < w; ++i) r[i] = RandomEntry(v, blk * w + i);
        // Left.
        BlockApply<1>(a, w, r, p0);
        for (int i = 0; i < w; ++i) lt[i * moments] = std::conj(p0[i]);
        if (moments > 1) {
          BlockApply<1>(h, w, p0, p1);
          for (int i = 0; i < w; ++i) lt[i * moments + 1] = std::conj(p1[i]);
        }
        cd* prev = p0;
        cd* cur = p1;
        for (int k = 2; k < moments; ++k) {
          BlockApply<1>(h, w, cur, y);
          for (int i = 0; i < w; ++i) {
            prev[i] = 2.0 * y[i] - prev[i];
            lt[i * moments + k] = std::conj(prev[i]);
          }
          std::swap(prev, cur);
        }
        // Right.
        for (int i = 0; i < w; ++i) p0[i] = r[i];
        BlockApply<1>(b, w, p0, y);
        for (int i = 0; i < w; ++i) rt[i * moments] = y[i];
        if (moments > 1) {
          BlockApply<1>(h, w, p0, p1);
          BlockApply<1>(b, w, p1, y);
          for (int i = 0; i < w; ++i) rt[i * moments + 1] = y[i];
        }
        prev = p0;
        cur = p1;
        for (int k = 2; k < moments; ++k) {
          BlockApply<1>(h, w, cur, y);
          for (int i = 0; i < w; ++i) prev[i] = 2.0 * y[i] - prev[i];
          std::swap(prev, cur);
          BlockApply<1>(b, w, cur, y);
          for (int i = 0; i < w; ++i) rt[i * moments + k] = y[i];
        }
        // Local table.
        for (int q = 0; q < w; ++q) {
          const cd* lq = lt.data() + q * moments;
          const cd* rq = rt.data() + q * moments;
          for (int mm = 0; mm < moments; ++mm) {
            const cd l = lq[mm];
            cd* row = acc.data() + static_cast<int64_t>(mm) * moments;
            for (int nn = 0; nn < moments; ++nn) row[nn] += l * rq[nn];
          }
        }
      }
    }
  }
  std::vector<cd> mu(tableSize, 0.0);
  for (int t = 0; t < threads; ++t)
    for (size_t k = 0; k < tableSize; ++k) mu[k] += part[t][k];
  *tAll = Seconds(t0, Clock::now());
  return mu;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 9) {
    std::fprintf(stderr,
                 "usage: %s <dos|cond> <global|local> <block> <dimension> <moments> <vectors> "
                 "<coupling 0|1> <out_moments_file>\n",
                 argv[0]);
    return 2;
  }
  const std::string mode = argv[1], method = argv[2];
  const int w = std::atoi(argv[3]);
  const int64_t dimension = std::atoll(argv[4]);
  const int moments = std::atoi(argv[5]);
  const int vectors = std::atoi(argv[6]);
  const bool coupled = std::atoi(argv[7]) != 0;
  if (w < 1 || w > kMaxBlock || vectors < 1 || vectors > kMaxVectors || moments < 2) {
    std::fprintf(stderr, "invalid block, vectors or moments\n");
    return 2;
  }
  const auto tb0 = Clock::now();
  const Model m = BuildModel(w, dimension, coupled);
  const double tBuild = Seconds(tb0, Clock::now());

  double tRec = 0.0, tTab = 0.0;
  std::vector<cd> mu;
  const auto t0 = Clock::now();
  if (mode == "dos" && method == "global") {
    switch (vectors) {
      case 1: mu = GlobalDos<1>(m, moments, &tRec); break;
      case 2: mu = GlobalDos<2>(m, moments, &tRec); break;
      case 4: mu = GlobalDos<4>(m, moments, &tRec); break;
      case 8: mu = GlobalDos<8>(m, moments, &tRec); break;
      default: std::fprintf(stderr, "global dos needs vectors in {1,2,4,8}\n"); return 2;
    }
  } else if (mode == "dos" && method == "local") {
    switch (vectors) {
      case 1: mu = LocalDos<1>(m, moments, &tRec); break;
      case 2: mu = LocalDos<2>(m, moments, &tRec); break;
      case 4: mu = LocalDos<4>(m, moments, &tRec); break;
      case 8: mu = LocalDos<8>(m, moments, &tRec); break;
      default: std::fprintf(stderr, "local dos needs vectors in {1,2,4,8}\n"); return 2;
    }
  } else if (mode == "cond" && method == "global") {
    mu = GlobalCond(m, moments, vectors, &tRec, &tTab);
  } else if (mode == "cond" && method == "local") {
    mu = LocalCond(m, moments, vectors, &tRec);
  } else {
    std::fprintf(stderr, "unknown mode or method\n");
    return 2;
  }
  const double tTotal = Seconds(t0, Clock::now());

  const double norm = 1.0 / (static_cast<double>(m.n) * vectors);
  for (auto& v : mu) v *= norm;
  if (FILE* f = std::fopen(argv[8], "wb")) {
    std::fwrite(mu.data(), sizeof(cd), mu.size(), f);
    std::fclose(f);
  }
  std::printf(
      "{\"mode\":\"%s\",\"method\":\"%s\",\"block\":%d,\"dimension\":%lld,\"moments\":%d,"
      "\"vectors\":%d,\"coupled\":%s,\"coupling_links\":%lld,\"h_nonzeros\":%lld,"
      "\"threads\":%d,\"build_s\":%.4f,\"total_s\":%.4f,\"recursion_s\":%.4f,\"table_s\":%.4f,"
      "\"mu0_re\":%.15e,\"mu1_re\":%.15e}\n",
      mode.c_str(), method.c_str(), w, static_cast<long long>(m.n), moments, vectors,
      coupled ? "true" : "false", static_cast<long long>(m.couplingLinks),
      static_cast<long long>(m.hPtr[m.n]), omp_get_max_threads(), tBuild, tTotal, tRec, tTab,
      mu[0].real(), mu[1].real());
  return 0;
}
