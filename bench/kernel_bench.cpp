// Item 5 measurement driver, not a gate. Tester records wall seconds, modeled
// FLOP/byte rates and output hash at each thread count. Setup and one warm-up
// are excluded. Rates use the documented algorithmic stream model, not DRAM PMCs.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <omp.h>
#include "core/algebra/block_diagonal_operator.hpp"
#include "core/algebra/build_profile.hpp"
#include "core/algebra/linalg.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/work_partition.hpp"
#include "core/chebyshev/contraction.hpp"

using namespace lsquant;

static std::uint64_t Hash(const void* data, std::size_t bytes) {
  const auto* p = static_cast<const unsigned char*>(data);
  std::uint64_t hash = 14695981039346656037ull;
  for (std::size_t i = 0; i < bytes; ++i) hash = (hash ^ p[i]) * 1099511628211ull;
  return hash;
}

template <typename Function>
static double Time(int repeats, Function run) {
  run();
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < repeats; ++i) run();
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

int main(int argc, char** argv) {
  try {
    if (argc != 5) throw std::invalid_argument("usage: stage8_kernels spmv|spmm|block|dot|contract SIZE BATCH REPEATS (block W=72)");
    const std::string mode = argv[1];
    const Index size = std::stoll(argv[2]);
    const int batch = std::stoi(argv[3]), repeats = std::stoi(argv[4]);
    if (size <= 0 || batch <= 0 || repeats <= 0) throw std::invalid_argument("positive sizes required");
    if (const char* value = std::getenv("LSQUANT_SPMV_MIN_NNZ")) SpmvThresholdSetting() = std::stoll(value);
    if (const char* value = std::getenv("LSQUANT_KERNEL_BATCH")) KernelBatchSetting() = std::stoi(value);
    if (SpmvThresholdSetting() < 0 || KernelBatchSetting() < 0) throw std::invalid_argument("negative tuning setting");
    double seconds = 0.0, flops = 0.0, bytes = 0.0;
    std::uint64_t hash = 0;
    if (mode == "contract") {
      const Index points = batch, m = size;
      std::vector<Scalar> mu(m * m), b(m * points);
      std::vector<double> a(m * points), output(points);
      for (Index i = 0; i < m * m; ++i) mu[i] = Scalar(std::sin(i * 0.1), std::cos(i * 0.2)) / double(m);
      for (Index i = 0; i < m * points; ++i) { a[i] = std::cos(i * 0.3); b[i] = Scalar(a[i], std::sin(i * 0.3)); }
      seconds = Time(repeats, [&] { ContractRank2Imag(mu.data(), m, m, a.data(), b.data(), points, output.data()); });
      flops = 8.0 * m * m * points;
      bytes = 16.0 * m * m + 24.0 * m * points + 8.0 * points;
      hash = Hash(output.data(), output.size() * sizeof(double));
    } else if (mode == "dot") {
      std::vector<Vector> a(batch, Vector(size)), b(batch, Vector(size));
      std::vector<const Vector*> left, right;
      for (int j = 0; j < batch; ++j) {
        for (Index i = 0; i < size; ++i) { a[j][i] = Scalar(std::sin(i * 0.1 + j), 0.2); b[j][i] = Scalar(0.3, std::cos(i * 0.2 + j)); }
        left.push_back(&a[j]); right.push_back(&b[j]);
      }
      std::vector<Scalar> output(batch * batch);
      seconds = Time(repeats, [&] {
#ifdef LSQUANT_NATIVE_BUILD
        output = linalg::MomentTable(left, right);
#else
        for (int j = 0; j < batch; ++j) for (int k = 0; k < batch; ++k) output[j * batch + k] = linalg::dot(a[j], b[k]);
#endif
      });
      flops = 8.0 * size * batch * batch;
      bytes = 32.0 * size * batch + 16.0 * batch * batch;
      hash = Hash(output.data(), output.size() * sizeof(Scalar));
    } else {
      SparseMatrix sparse;
      std::unique_ptr<LinearOperator> op;
      Index dimension = size, nnz = 0;
      if (mode == "block") {
        constexpr Index w = 72;
        dimension = size * w; nnz = size * w * w;
        std::vector<Scalar> blocks(nnz);
        for (Index k = 0; k < size; ++k) for (Index i = 0; i < w; ++i) for (Index j = 0; j < w; ++j) {
          blocks[(k * w + i) * w + j] = Scalar(i == j ? 2.0 : std::cos(double(i + j + k)) / w, 0.0);
        }
        op = std::make_unique<BlockDiagonalOperator>(size, w, std::move(blocks));
      } else if (mode == "spmv" || mode == "spmm") {
        std::vector<Index> rp(size + 1), cols;
        std::vector<Scalar> values;
        for (Index row = 0; row < size; ++row) {
          rp[row] = static_cast<Index>(cols.size());
          for (Index col = std::max<Index>(0, row - 4); col < std::min(size, row + 5); ++col) {
            cols.push_back(col); values.emplace_back(0.1 * (1 + row - col), 0.03 * (row - col));
          }
        }
        rp[size] = static_cast<Index>(cols.size()); nnz = rp[size];
        sparse = SparseMatrix::FromCsr(size, std::move(rp), std::move(cols), std::move(values));
        op = std::make_unique<SparseMatrixOperator>(sparse);
      } else throw std::invalid_argument("unknown mode");
      const int count = mode == "spmv" ? 1 : batch;
      std::vector<Vector> input(count, Vector(dimension)), output(count, Vector(dimension));
      for (int j = 0; j < count; ++j) for (Index i = 0; i < dimension; ++i) input[j][i] = Scalar(std::sin(i * 0.1 + j), std::cos(i * 0.2 + j));
      seconds = Time(repeats, [&] {
        if (mode == "spmv") op->Multiply(Scalar(1, 0), input[0], Scalar(0, 0), output[0]);
        else op->MultiplyBatch(Scalar(1, 0), input, Scalar(0, 0), output);
      });
      int sweeps = count;
#ifdef LSQUANT_NATIVE_BUILD
      if (mode != "block") sweeps = (count + CacheBatchSize(dimension, count) - 1) / CacheBatchSize(dimension, count);
      else sweeps = 1; // Lower-bound block stream; microbatch rereads may hit cache.
#endif
      flops = 8.0 * nnz * count;
      bytes = (mode == "block" ? 16.0 * nnz : 24.0 * nnz + 8.0 * (dimension + 1)) * sweeps + 48.0 * dimension * count;
      for (const auto& v : output) hash ^= Hash(v.data(), v.size() * sizeof(Scalar));
    }
    std::cout << std::setprecision(17) << "{\"mode\":\"" << mode << "\",\"profile\":\"" << BuildProfileName()
              << "\",\"threads\":" << omp_get_max_threads() << ",\"size\":" << size << ",\"batch\":" << batch
              << ",\"repeats\":" << repeats << ",\"wall_seconds\":" << seconds
              << ",\"modeled_flops_per_repeat\":" << flops << ",\"modeled_bytes_per_repeat\":" << bytes
              << ",\"modeled_gflops\":" << flops * repeats / seconds / 1e9
              << ",\"modeled_gbs\":" << bytes * repeats / seconds / 1e9 << ",\"hash\":" << hash << "}\n";
    return 0;
  } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
