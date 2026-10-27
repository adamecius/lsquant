# Handoff 0001 (Agi) — stage 1: the engine, rank 1 and rank 2 moments in real space

Roles: Agi implements, Codex (on TCN, `codex.md`) builds the harness, runs the gate
byte for byte and diagnoses, Fable checks and returns findings. Read `AGENTS.md`
first; sections 2, 3, 8, 9, 10 and 11 apply to every line below. Nothing here is
optional. A deviation goes in `report.md`, never silently in code.

Stage 0 (goldens from the legacy build, committed under `goldens/stage0/`, see
`handoff/0000-goldens/codex.md`) is closed before this stage starts.

## 0. Scope

The C++ engine of `AGENTS.md` section 3 (`Plan` of five instructions, basis-free),
with the real-space CSR operator, the trace-vector generators, the Chebyshev
recursion and batching, the rank-generic moment tensor, the CSR reader and the
`.chebmom` writer/reader. Deliverable: the two stage-1 plans (`dos`, `kubo`) on the
chain1d goldens G1 and G2 are **byte-identical** to the legacy files.

Explicitly out of scope, do not start, do not leave hooks or TODOs for:
`BlochTransform`, `BlockDiagonalMatrix`, `DiagonalMatrix`, any `MetricPolicy` beyond
nothing at all (there is no metric in this stage), `CoefficientTable`, pybind11,
parsers, CLI, FFT, any physics name in C++ (`VX`, `delta`, `Kubo`, `DOS` appear only
in `apps/stage1_moments.cpp` as command-line tokens and in file names).

## 1. Where things go

```
lsquant/                              new tree, C++17, namespace lsquant
  CMakeLists.txt                      STATIC library lsquant_core; executables stage1_moments, stage_1_check
  core/algebra/
    types.hpp                         Index, Scalar, Vector
    linear_operator.hpp               LinearOperator
    sparse_matrix.hpp .cpp            SparseMatrix, SparseMatrixOperator
    composite_operators.hpp .cpp      IdentityOperator, ScaledOperator, OperatorSum, OperatorProduct, AdjointOperator
    linalg.hpp .cpp                   axpy, scal, copy, dot, norm2
    state_generator.hpp .cpp          StateGenerator, RandomPhaseStates, LocalStates, UserStates
    spectral_bounds.hpp .cpp          SpectralBounds, LanczosBounds, GershgorinBounds
  core/chebyshev/
    spectral_window.hpp               SpectralWindow
    rescaled_operator.hpp .cpp        RescaledOperator
    chebyshev_recursion.hpp .cpp      ChebyshevRecursion
    chebyshev_vectors.hpp .cpp        ChebyshevVectors
    chebyshev_moments.hpp .cpp        ChebyshevMoments, ApplyKpmWeights
    plan.hpp                          ApplyOperator, Expand, ApplySeries, MirrorToBra, Trace, Instruction, Plan
    moment_engine.hpp .cpp            MomentEngine
  io/
    csr_text.hpp .cpp                 ReadCsrText, WriteCsrText
    chebmom_file.hpp .cpp             WriteChebmom, ReadChebmom, WriteLegacyChebmom1D, WriteLegacyChebmom2D
  apps/stage1_moments.cpp             the executable the gate runs (section 6)
tests/dev/stage_1_check.cpp           the single self-check (section 5), not a gate
```

The legacy tree (`src/`, `include/`, `test/`, `examples/`) is not modified. The root
`CMakeLists.txt` gains exactly two lines: `option(LSQUANT_BUILD_NEW_CORE "" ON)` and
the guarded `add_subdirectory(lsquant)`. Build: `cmake -S . -B build && cmake --build
build -j` from the repository root, on TCN.

## 2. Restrictions

1. C++17. Dependencies: Eigen3 and OpenMP only. No FFTW, GSL, Boost, MKL, nothing
   vendored, no `FetchContent`.
2. Under `core/` and `io/` there is no `std::cout`, `std::cerr`, `printf`, `getenv`,
   `std::ifstream`/`ofstream` (except inside `io/`), no `assert` that is compiled out
   in Release: errors throw `std::invalid_argument` or `std::runtime_error` with a
   message that names the function.
3. No global mutable state, no singletons, no static caches, no environment
   variables (`BATCH_SIZE`, `KPM_SEED`, `LSQUANT_*` are read only by
   `apps/stage1_moments.cpp` as command-line options).
4. `using Index = std::int64_t; using Scalar = std::complex<double>; using Vector =
   std::vector<Scalar>;` in `types.hpp`; nothing else is used for these.
5. Every public class, method and free function carries the comment of
   `agent_template/documentation_incode.md`: problem solved, approach, assumptions,
   units, complexity, pitfalls. A missing comment is a defect.
6. Names exactly as in this brief. PascalCase classes and methods, snake_case files
   and free functions, lower-case namespaces. No abbreviations.
7. Determinism: the same inputs and seed give bit-identical moments at any OpenMP
   thread count. The sum over trace vectors is sequential and in generator order;
   parallelism only inside the matvec (row parallel, each row summed in column order
   by one thread).
8. Byte identity with the legacy goldens is the gate, so the legacy accumulation
   order is reproduced exactly (section 4, "conventions"). Where the legacy does
   something numerically odd, the odd thing is reproduced and noted in `report.md`;
   it is not improved here.
9. Agi does not edit `goldens/`, `scripts/check_goldens.sh`, `figures/`, and does not
   run `check_goldens.sh`. Agi does not write "done".

## 3. Interfaces

```cpp
namespace lsquant {
using Index  = std::int64_t;  using Scalar = std::complex<double>;  using Vector = std::vector<Scalar>;

// ---- core/algebra ----------------------------------------------------------------------
class LinearOperator {
public:
  virtual ~LinearOperator() = default;
  virtual Index Dimension() const = 0;
  virtual void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const = 0;  // y = a*O*x + b*y
  void Multiply(const Vector& x, Vector& y) const;                                  // a = 1, b = 0
  virtual bool IsIdentity() const { return false; }
  virtual bool IsHermitian() const = 0;
};

class SparseMatrix {                                   // CSR, Eigen::SparseMatrix<Scalar, RowMajor, Index>
public:
  static SparseMatrix FromCsr(Index dim, std::vector<Index> rowPointers,
                              std::vector<Index> columns, std::vector<Scalar> values);
  static SparseMatrix FromTriplets(Index dim, const std::vector<Index>& rows,
                              const std::vector<Index>& cols, const std::vector<Scalar>& vals);
  Index Dimension() const;  Index NonZeros() const;
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const;
  SparseMatrix Rescaled(Scalar a, Scalar b) const;     // returns a*M + b*I as a NEW matrix, same Eigen expression as legacy Rescale
  SparseMatrix Adjoint() const;
  bool IsHermitian(double tol = 1e-12) const;
  const Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index>& Eigen() const;
};
class SparseMatrixOperator : public LinearOperator { /* non-owning view of a SparseMatrix */ };
class IdentityOperator : public LinearOperator { /* IsIdentity() true; Multiply does y = a*x + b*y */ };
class ScaledOperator   : public LinearOperator { /* alpha * O */ };
class OperatorSum      : public LinearOperator { /* O1 + O2 */ };
class OperatorProduct  : public LinearOperator { /* O1 * O2, O2 applied first */ };
class AdjointOperator  : public LinearOperator { /* O^+; for SparseMatrixOperator holds the adjoint CSR */ };

namespace linalg {                                      // thin wrappers over Eigen::Map, same calls as legacy linalg
  void axpy(Scalar a, const Vector& x, Vector& y);  void scal(Scalar a, Vector& x);  void copy(const Vector& x, Vector& y);
  Scalar dot(const Vector& x, const Vector& y);         // Eigen Map(x).dot(Map(y)): conj(x) . y
  double norm2(const Vector& x);
}

class StateGenerator {
public:
  virtual ~StateGenerator() = default;
  virtual Index Dimension() const = 0;  virtual Index NumberOfStates() const = 0;
  virtual bool Next(Vector& out) = 0;   virtual void Reset() = 0;
  virtual std::string Label() const = 0;               // "default" | "exact" | "local" | "user"
};
class RandomPhaseStates : public StateGenerator {      // ctor (dimension, numberOfStates, seed = 12345)
  // Exact legacy convention (include/quantum_states.hpp, case RANDOM_STATE):
  //   norm = sqrt(N);  for i in 0..N-1: phi = 2*pi*uni(rng);  x_i = Scalar(cos(phi)/norm, sin(phi)/norm)
  // with std::mt19937 rng(seed) and std::uniform_real_distribution<double> uni(0.0, 1.0),
  // the distribution object created once per generator. Label "default".
};
class LocalStates : public StateGenerator {            // ctor (dimension, std::vector<Index> sites); Label "exact" iff sites == 0..N-1, else "local"
};
class UserStates  : public StateGenerator {            // ctor (std::vector<Vector>); Label "user"
};

struct SpectralBounds { double lower, upper; };
SpectralBounds LanczosBounds(const LinearOperator& H, int iterations = 40, unsigned seed = 12345);
SpectralBounds GershgorinBounds(const SparseMatrix& H, double pad = 0.10);

// ---- core/chebyshev --------------------------------------------------------------------
struct SpectralWindow {
  double lower, upper;
  double HalfWidth() const { return 0.5 * (upper - lower); }
  double Center() const    { return 0.5 * (upper + lower); }
  double Scale() const     { return 1.0 / HalfWidth(); }            // legacy ScaleFactor, CUTOFF = 1
  double Shift() const     { return -Center() / HalfWidth(); }      // legacy ShiftFactor
};

class RescaledOperator : public LinearOperator {
  // H~ = Scale*H + Shift*I. Holds its OWN rescaled CSR copy built with SparseMatrix::Rescaled
  // (the legacy rescales in place and then multiplies; a matrix-free a*(Hx) + b*x rounds
  // differently and breaks byte identity). The caller's H is never modified.
public:
  RescaledOperator(const SparseMatrix& H, const SpectralWindow& window);
};

class ChebyshevRecursion {                             // T_{m+1} = 2 H~ T_m - T_{m-1}
public:
  explicit ChebyshevRecursion(const LinearOperator& rescaledH);
  void Start(const Vector& phi);                       // T_0 = phi; T_1 = H~ phi  (legacy SetInitVectors)
  void Step();                                         // legacy Iterate: multiply(2, T_m, -1, T_{m-1}); swap
  const Vector& Current() const;  const Vector& Previous() const;  int Order() const;
};

class ChebyshevVectors {                               // batch of T_m phi for m in [first, first + Size)
public:
  ChebyshevVectors(Index dimension, int batchSize);
  void Fill(ChebyshevRecursion& recursion);            // legacy IterateAll: copy Current, Step, ...
  void ApplyToAll(const LinearOperator& A);            // legacy Vectors::Multiply(OP): copy then OP.Multiply
  Vector& At(int m);  const Vector& At(int m) const;  int Size() const;  int FirstOrder() const;
};

class ChebyshevMoments {
public:
  explicit ChebyshevMoments(std::vector<int> dims);    // rank = dims.size(); 0, 1 or 2 in this stage
  int Rank() const;  const std::vector<int>& Dims() const;
  Scalar& At(const std::vector<int>& index);  Scalar& At(int m);  Scalar& At(int m, int n);   // + const
  std::vector<Scalar>& Data();                         // row-major
  SpectralWindow window;  Index systemSize = 0;  int numberOfStates = 0;
  std::string stateLabel, word;  bool kpmWeighted = false;  bool symmetrised = false;
};
// Legacy post-processing of CorrelationExpansionMoments, verbatim (see section 4): rank 2 only.
void ApplyKpmWeights(ChebyshevMoments& mu, int numberOfStates, bool symmetrise);

// ---- plan.hpp: the five instructions (AGENTS.md section 3) ----------------------------
struct ApplyOperator { int operatorIndex; };
struct Expand        { int numMoments; };
struct ApplySeries   { std::vector<Scalar> coeffs; };  // present in the type; NOT implemented in this stage: Run throws std::runtime_error("ApplySeries: not implemented in stage 1")
struct MirrorToBra   {};
struct Trace         { std::vector<Scalar> weights; }; // size 1 (scalar weight) or the tensor size
using Instruction = std::variant<ApplyOperator, Expand, ApplySeries, MirrorToBra, Trace>;
struct Plan { std::vector<Instruction> code; std::string word; };

class MomentEngine {
public:
  MomentEngine(const LinearOperator& rescaledH, std::vector<const LinearOperator*> operators,
               StateGenerator& states, int batchSize);
  ChebyshevMoments Run(const Plan& plan);
};
}
```

## 4. Semantics of `MomentEngine::Run` and the legacy conventions it must reproduce

Read `src/chebyshev_solver.cpp` (`SpectralMoments`, `CorrelationExpansionMoments`,
`ComputeMomTable`), `src/chebyshev_moments.cpp` (`SetInitVectors`, `Iterate`),
`src/chebyshev_vectors.cpp` (`IterateAll`, `Multiply`), `src/eigen_sparse_matrix.cpp`
(`Multiply`, `Rescale`) and `src/eigen_linear_algebra.cpp` (`vdot`) before coding.
The gate is byte identity; these are the facts it depends on.

**State machine.** For each trace vector `r` (in generator order): the ket side
starts as the single live vector `r`; instructions act on the ket side until
`MirrorToBra`, after which they act on the bra side, which starts as a fresh copy
of `r`. Each side holds a list of live vectors and a list of opened indices.

- `ApplyOperator(k)`: every live vector `v <- O_k v` (identity: no operation). On
  a batch, this is `ChebyshevVectors::ApplyToAll` (legacy: copy into a temporary,
  then `OP.Multiply(temp, v)`).
- `Expand(M)`: for every live vector `v`, run `ChebyshevRecursion` from `v` and
  replace `v` by the `M` vectors `T_0 v .. T_{M-1} v`, opening one index of length
  `M` on that side. In this stage a side may hold at most one open index (a second
  `Expand` on the same side throws `std::runtime_error("Expand: rank > 2 not
  implemented in stage 1")`). Batching: `M` is processed in blocks of `batchSize`
  vectors exactly as `CorrelationExpansionMoments` does (`mL`, `mR` loops), so that
  the dot products are accumulated in the legacy block order.
- `MirrorToBra`: switch sides. At most once per plan.
- `Trace(w)`: for the opened indices, `mu[i_bra..., i_ket...] += w[...] *
  dot(bra_{i_bra}, ket_{i_ket})`, with `dot = conj(bra) . ket` (`linalg::dot`,
  i.e. the legacy `vdot(left, right)`). `w` of size 1 is a scalar weight; otherwise
  `w` has the tensor shape and is applied per element. The accumulation into `mu`
  is `mu += w * dot` per trace vector (not `w * sum`), because that is what the
  legacy does for rank 1. Rank 0 (no index opened): the tensor has one element.

**The two stage-1 plans** (built by `apps/stage1_moments.cpp`, no Python yet) and the
legacy code each reproduces:

```
dos  (legacy SpectralMoments, OP = identity, Moments1D):
     ket: Expand(M); Trace(w) with w[m] = (2/N_states) * (0.5 if m == 0 else 1)
     => mu[m] += w[m] * dot(r, T_m r)      per state, m = 0..M-1 in order     [legacy: chebMoms(m) += scal*vdot(Phi, Chebyshev0)]
     ChebyshevMoments: rank 1, kpmWeighted = true, symmetrised = false, word "1 delta"

kubo (legacy CorrelationExpansionMoments(OPL = VX, OPR = VX), Moments2D):
     ket: Expand(M), ApplyOperator(VX)            -> right_n = VX T_n r         [legacy chevVecR: SetInitVectors(PhiR); IterateAll; Multiply(OPR)]
     MirrorToBra
     bra: ApplyOperator(VX), Expand(M)            -> left_m  = T_m VX r         [legacy chevVecL: SetInitVectors(OPL, PhiL): Chebyshev0 = OPL*T0; IterateAll]
     Trace(1)                                     -> mu[m, n] += dot(left_m, right_n)   [legacy ComputeMomTable + AddSubMatrix, raw sums over states]
     then ApplyKpmWeights(mu, N_states, symmetrise = true), verbatim legacy:
        for mL in 0..M-1: for mR in mL..M-1:
           scal = 4.0/N_states; if mL == 0: scal *= 0.5; if mR == 0: scal *= 0.5
           tmp = scal * (mu(mL,mR) + conj(mu(mR,mL))) / 2.0
           mu(mL,mR) = tmp; mu(mR,mL) = conj(tmp)
     ChebyshevMoments: rank 2, kpmWeighted = true, symmetrised = true, word "VX delta VX delta"
```

Details that decide byte identity:

- `Expand` for the bra in `kubo` starts the recursion from `VX r` (legacy
  `SetInitVectors(OP, T0)`: `Chebyshev0 = OP * T0; Chebyshev1 = H~ * Chebyshev0`),
  while the ket starts from `r` and applies `VX` to the finished batch. Both orders
  must be exactly these.
- The legacy rounds `M` up to a multiple of `batchSize` (default `batchSize = M`
  when `BATCH_SIZE` is unset) and truncates back; the gate runs with the default.
- Legacy `Rescale(a, b)`: `matrix = a*matrix + b*I` as an Eigen sparse expression;
  `SparseMatrix::Rescaled` uses the same expression on a copy.
- Legacy `Multiply(a, x, b, y)`: reproduce the exact Eigen expression and the
  per-row summation of `src/eigen_sparse_matrix.cpp` (it has a row-parallel OpenMP
  kernel that is bit-identical to the Eigen product at any thread count; the test
  `test/spmv_threaded_bitexact_test.cpp` documents it).
- Legacy `vdot(x, y)` is `Eigen::Map(x).dot(Eigen::Map(y))`; `linalg::dot` is the
  same call. `MomentTable` is not used in this stage: `Trace` does the per-pair
  `dot` in the legacy loop order (`m` outer, `n` inner). A GEMM table comes in
  stage 3 behind the same `Trace`.
- Legacy `Moments2D::saveIn` writes the header `size " " bandwidth "  " center`
  (two spaces before the center), the dims line `"64 64 "` (trailing space), and
  values with `precision(std::numeric_limits<double>::digits10)` as `re " " im`,
  so `-0` appears. `Moments1D::saveIn`: `size " " bandwidth " " center`, then the
  moment count, then the pairs. `WriteLegacyChebmom1D/2D` reproduce these byte for
  byte; `cmp` is the judge.
- `BandWidth` in the legacy header is `upper - lower` and `BandCenter` is
  `(upper + lower)/2`, written as computed by `(b[1]-b[0])*1.0` and
  `(b[1]+b[0])*0.5`.

## 5. Implementation steps and the one self-check

1. `types.hpp`, `SparseMatrix`, `SparseMatrixOperator`, `linalg`, `io/csr_text`.
2. Composite operators.
3. `StateGenerator` family.
4. `SpectralBounds`.
5. `SpectralWindow`, `RescaledOperator`.
6. `ChebyshevRecursion`, `ChebyshevVectors`.
7. `ChebyshevMoments`, `ApplyKpmWeights`, `io/chebmom_file` (new format of section
   7 plus the two legacy writers and a reader for all three layouts).
8. `plan.hpp`, `MomentEngine::Run` with the state machine of section 4.
9. `apps/stage1_moments.cpp` (section 6).
10. `tests/dev/stage_1_check.cpp`, one executable, no framework: a random complex
    Hermitian `64 x 64` matrix from `std::mt19937(7)` as CSR; dense reference with
    `Eigen::SelfAdjointEigenSolver`; prints the max abs error of the `dos` and
    `kubo` plans by exact trace against `sum_k T_m(eps~_k)` and
    `sum_k sum_l A_kl T_m(eps~_k) A_lk T_n(eps~_l)` with the section-4 weights, and
    whether the tensors are bit-identical at `OMP_NUM_THREADS=1` and `8` (it
    re-executes itself with the environment set). Expected: below 1e-11 and
    identical. This is Agi's tool while implementing; it is not a gate and is not
    extended.

## 6. The gate executable

```
stage1_moments --operators <dir> --label chain1d --window <lower> <upper>
               --trace exact|random --seed 12345 --states <N> --moments 64
               --word dos|kubo --batch <M> --out <dir>
```

Reads `<dir>/<label>.HAM.CSR` (and `.VX.CSR` for `kubo`), builds the plan of
section 4, runs it, and writes with the legacy layout and the legacy file names:

- `dos`:  `<out>/SpectralOp1<label>KPM_M<M>_state<label_of_states>.chebmom1D`
- `kubo`: `<out>/NonEqOpVX-VX<label>KPM_M<M>x<M>_state<label_of_states>.chebmom2D`

where the state label is `exact` for `--trace exact` and `default` for random
phase (legacy `gen.StateLabel()`). It also writes the same tensor in the new
`.chebmom` format next to it. Exit code non-zero on any exception, with the
message on `stderr`. This is the only program the gate runs; Codex's harness calls
it on `goldens/stage0/chain1d/operators` with `--window -2 2 --trace exact
--states 2048 --moments 64` and compares with `cmp`.

## 7. The `.chebmom` format (new, rank generic)

```
chebmom 1
rank 2
dims 64 64
window -2 2
system_size 2048
states exact 2048
word VX delta VX delta
kpm_weighted 1
symmetrised 1
data
re im            # product(dims) lines, row-major, precision 17
```

## 8. Report

`handoff/0001-core-algebra-and-chebyshev/report.md`, in this order: files added;
exact build commands and their output on TCN; the output of `stage_1_check` at 1
and 8 threads; the output of `stage1_moments` for both words on
`goldens/stage0/chain1d/operators`; every deviation from this brief with its
reason; open questions. No "done", no gate run, no edits outside the listed paths.
