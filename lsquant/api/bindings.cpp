/**
 * @file bindings.cpp
 * @brief pybind11 module `lsquant._core` exposing the C++ Chebyshev engine to Python.
 *
 * Binds `SparseMatrix`, the `LinearOperator` hierarchy, `StateGenerator` implementations,
 * `SpectralBounds` and `GershgorinBounds`, `SpectralWindow`, `RescaledOperator`, the five
 * bytecode instructions and `Plan`, `MomentEngine`, `ChebyshevMoments` (with a zero-copy
 * NumPy `complex128` view), file I/O routines, rank-1 and rank-2 series contractions
 * (`ContractRank1`, `ContractRank2Imag`, `ContractRank2Real`, `ContractRank2RealReal`,
 * `WriteCurveText`), and the gate-only `_core.compat` submodule.
 */
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "compat/legacy_kernels.hpp"
#include "core/algebra/bloch_transform.hpp"
#include "core/algebra/block_diagonal_operator.hpp"
#include "core/algebra/composite_operators.hpp"
#include "core/algebra/linear_operator.hpp"
#include "core/algebra/metric_policy.hpp"
#include "core/algebra/overlap_appliers.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/spectral_bounds.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/algebra/types.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/contraction.hpp"
#include "core/algebra/fft_moment_contraction.hpp"
#include "core/algebra/work_partition.hpp"
#include "core/chebyshev/moment_engine.hpp"
#include "core/chebyshev/plan.hpp"
#include "core/chebyshev/rescaled_operator.hpp"
#include "core/chebyshev/spectral_window.hpp"
#include "io/chebmom_file.hpp"
#include "io/csr_text.hpp"
#include "io/curve_text.hpp"

namespace py = pybind11;

namespace lsquant {
namespace {

/**
 * @brief Validate that a Python object is a 1D C-contiguous NumPy array of exact element type `T`.
 *
 * Rejects non-arrays, non-1D arrays, non-contiguous slices, or mismatched dtypes without
 * implicit conversion (`forcecast` off), raising a Python `TypeError` that names `argName`.
 *
 * @tparam T            Expected C++ element type (`Index` or `Scalar`).
 * @param obj           Input Python handle.
 * @param argName       Parameter name included in the `TypeError` message on failure.
 * @param expectedDtype Human-readable NumPy dtype name (`"int64"` or `"complex128"`).
 * @param expectedKind  NumPy dtype kind character (`'i'` for signed integer, `'c'` for complex).
 * @param outSize       Output element count.
 * @return              Const pointer to the first element of the validated buffer.
 */
template <typename T>
const T* ValidateCsrBuffer1D(py::handle obj, const char* argName, const char* expectedDtype,
                             char expectedKind, Index& outSize) {
  if (!py::isinstance<py::array>(obj)) {
    throw py::type_error(std::string("SparseMatrix.from_csr: argument '") + argName +
                         "' must be a 1D C-contiguous " + expectedDtype + " numpy.ndarray");
  }
  auto arr = py::reinterpret_borrow<py::array>(obj);
  const bool dtypeOk = (arr.dtype().kind() == expectedKind) &&
                       (arr.itemsize() == static_cast<py::ssize_t>(sizeof(T)));
  const bool shapeOk = (arr.ndim() == 1);
  const bool contigOk =
      ((arr.flags() & py::array::c_style) != 0) &&
      (arr.ndim() == 1 &&
       (arr.size() <= 1 || arr.strides(0) == static_cast<py::ssize_t>(sizeof(T))));
  if (!dtypeOk || !shapeOk || !contigOk) {
    throw py::type_error(std::string("SparseMatrix.from_csr: argument '") + argName +
                         "' must be a 1D C-contiguous " + expectedDtype + " numpy.ndarray");
  }
  outSize = static_cast<Index>(arr.size());
  return static_cast<const T*>(arr.data());
}

/**
 * @brief Validate a 2D `float64` NumPy array and contract it with a strided real moment vector.
 *
 * @param momentsPtr   Pointer to the first real moment element.
 * @param momentStride Element stride between consecutive moments.
 * @param numMoments   Number of moments `M`.
 * @param table        2D `float64` NumPy array of shape `(M, N)`.
 * @return             1D `float64` NumPy array of length `N`.
 */
py::array_t<double> ContractRank1Buffer(const double* momentsPtr, Index momentStride,
                                        Index numMoments, const py::array_t<double>& table) {
  if (table.ndim() != 2) {
    throw std::invalid_argument("ContractRank1: table must be a 2D float64 array of shape (M, N)");
  }
  if (static_cast<Index>(table.shape(0)) != numMoments) {
    throw std::invalid_argument(
        "ContractRank1: table.shape[0] must equal the number of moments M");
  }
  const Index numPoints = static_cast<Index>(table.shape(1));
  const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
  if (table.strides(0) % elemSize != 0 || table.strides(1) % elemSize != 0) {
    throw std::invalid_argument("ContractRank1: unaligned table strides");
  }
  const Index tableStrideM = static_cast<Index>(table.strides(0) / elemSize);
  const Index tableStrideN = static_cast<Index>(table.strides(1) / elemSize);

  py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
  ContractRank1(momentsPtr, momentStride, numMoments, table.data(), tableStrideM, tableStrideN,
                numPoints, result.mutable_data());
  return result;
}

template <typename TableAType, typename TableBType>
Index ValidateRank2Tables(const char* funcName, Index m0, Index m1, const TableAType& tableA,
                          const TableBType& tableB) {
  if (tableA.ndim() != 2 || tableB.ndim() != 2) {
    throw std::invalid_argument(std::string(funcName) + ": tables must be 2D arrays");
  }
  if (static_cast<Index>(tableA.shape(0)) != m0) {
    throw std::invalid_argument(std::string(funcName) +
                                ": table_a.shape[0] must match moments axis 0");
  }
  if (static_cast<Index>(tableB.shape(0)) != m1) {
    throw std::invalid_argument(std::string(funcName) +
                                ": table_b.shape[0] must match moments axis 1");
  }
  if (tableA.shape(1) != tableB.shape(1)) {
    throw std::invalid_argument(std::string(funcName) +
                                ": table_a.shape[1] and table_b.shape[1] must match");
  }
  return static_cast<Index>(tableA.shape(1));
}

}  // namespace
}  // namespace lsquant

PYBIND11_MODULE(_core, m) {
  using namespace lsquant;

  m.doc() = "C++ execution engine and zero-copy tensor bindings for LSQUANT.";
  m.attr("build_profile") = BuildProfileName();
  // Configure once before execution, never while an engine is running.
  m.def("ConfigurePerformance", [](Index threshold, int batch) {
    if (threshold < 0 || batch < 0) throw std::invalid_argument("ConfigurePerformance: negative setting");
    SpmvThresholdSetting() = threshold;
    KernelBatchSetting() = batch;
  }, py::arg("spmv_min_nnz") = 16384, py::arg("kernel_batch") = 0);

  m.def("FftMomentContraction", [](const py::array_t<Scalar, py::array::c_style>& moments,
                                    Index points) {
    if (moments.ndim() != 2) throw std::invalid_argument("FftMomentContraction: expected rank 2");
    const auto values = FftMomentContraction(moments.data(), moments.shape(0), moments.shape(1), points);
    py::array_t<Scalar> out({py::ssize_t(2), static_cast<py::ssize_t>(points)});
    std::copy(values.begin(), values.end(), out.mutable_data());
    return out;
  }, py::arg("moments"), py::arg("points"));

  py::class_<SparseMatrix>(m, "SparseMatrix")
      .def(py::init<>())
      .def_static(
          "from_csr",
          [](Index dim, py::object rowPointers, py::object columns, py::object values) {
            Index rpSize = 0;
            Index colSize = 0;
            Index valSize = 0;
            const Index* rpPtr =
                ValidateCsrBuffer1D<Index>(rowPointers, "row_pointers", "int64", 'i', rpSize);
            const Index* colPtr =
                ValidateCsrBuffer1D<Index>(columns, "columns", "int64", 'i', colSize);
            const Scalar* valPtr =
                ValidateCsrBuffer1D<Scalar>(values, "values", "complex128", 'c', valSize);
            return SparseMatrix::FromCsr(dim, rpSize, rpPtr, colSize, colPtr, valSize, valPtr);
          },
          py::arg("dim"), py::arg("row_pointers"), py::arg("columns"), py::arg("values"))
      .def_property_readonly("dimension", &SparseMatrix::Dimension)
      .def_property_readonly("non_zeros", &SparseMatrix::NonZeros)
      .def("hermitian_fold", &SparseMatrix::HermitianFold)
      .def("is_hermitian", &SparseMatrix::IsHermitian, py::arg("tol") = 1e-12);

  py::class_<LinearOperator>(m, "LinearOperator")
      .def_property_readonly("dimension", &LinearOperator::Dimension)
      .def_property_readonly("is_identity", &LinearOperator::IsIdentity)
      .def_property_readonly("is_hermitian", &LinearOperator::IsHermitian)
      .def(
          "multiply",
          [](const LinearOperator& op,
             const py::array_t<Scalar, py::array::c_style | py::array::forcecast>& x) {
            if (x.ndim() != 1) {
              throw std::invalid_argument("LinearOperator.multiply: x must be a 1D complex128 array");
            }
            const Index dim = op.Dimension();
            if (static_cast<Index>(x.size()) != dim) {
              throw std::invalid_argument("LinearOperator.multiply: x dimension mismatch");
            }
            Vector xVec(x.data(), x.data() + x.size());
            Vector yVec(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
            op.Multiply(Scalar(1.0, 0.0), xVec, Scalar(0.0, 0.0), yVec);
            py::array_t<Scalar> out(static_cast<py::ssize_t>(dim));
            if (dim > 0) {
              std::copy(yVec.begin(), yVec.end(), out.mutable_data());
            }
            return out;
          },
          py::arg("x"))
      .def(
          "multiply",
          [](const LinearOperator& op, Scalar a,
             const py::array_t<Scalar, py::array::c_style | py::array::forcecast>& x, Scalar b,
             py::array_t<Scalar, py::array::c_style>& y) {
            if (x.ndim() != 1 || y.ndim() != 1) {
              throw std::invalid_argument(
                  "LinearOperator.multiply: x and y must be 1D complex128 arrays");
            }
            const Index dim = op.Dimension();
            if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
              throw std::invalid_argument("LinearOperator.multiply: vector dimension mismatch");
            }
            Vector xVec(x.data(), x.data() + x.size());
            Vector yVec(y.data(), y.data() + y.size());
            op.Multiply(a, xVec, b, yVec);
            if (dim > 0) {
              std::copy(yVec.begin(), yVec.end(), y.mutable_data());
            }
          },
          py::arg("a"), py::arg("x"), py::arg("b"), py::arg("y"));

  py::class_<SparseMatrixOperator, LinearOperator>(m, "SparseMatrixOperator")
      .def(py::init<const SparseMatrix&>(), py::arg("matrix"), py::keep_alive<1, 2>());

  py::class_<IdentityOperator, LinearOperator>(m, "IdentityOperator")
      .def(py::init<Index>(), py::arg("dim"));

  py::class_<DiagonalOperator, LinearOperator>(m, "DiagonalOperator")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style | py::array::forcecast>& vals) {
             if (vals.ndim() != 1) {
               throw std::invalid_argument(
                   "DiagonalOperator: values must be a 1D complex128/float64 array");
             }
             return DiagonalOperator(std::vector<Scalar>(vals.data(), vals.data() + vals.size()));
           }),
           py::arg("values"))
      .def(py::init<std::vector<Scalar>>(), py::arg("values"))
      .def("adjoint", &DiagonalOperator::Adjoint);

  py::class_<ScaledOperator, LinearOperator>(m, "ScaledOperator")
      .def(py::init<Scalar, const LinearOperator&>(), py::arg("alpha"), py::arg("op"),
           py::keep_alive<1, 3>());

  py::class_<OperatorSum, LinearOperator>(m, "OperatorSum")
      .def(py::init<const LinearOperator&, const LinearOperator&>(), py::arg("first"),
           py::arg("second"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<OperatorProduct, LinearOperator>(m, "OperatorProduct")
      .def(py::init<const LinearOperator&, const LinearOperator&>(), py::arg("first"),
           py::arg("second"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<AdjointOperator, LinearOperator>(m, "AdjointOperator")
      .def(py::init<const LinearOperator&>(), py::arg("op"), py::keep_alive<1, 2>())
      .def(py::init<const SparseMatrix&>(), py::arg("matrix"));

  py::class_<SandwichOperator, LinearOperator>(m, "SandwichOperator")
      .def(py::init<const LinearOperator&, const LinearOperator&>(), py::arg("B"), py::arg("V"),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<LowdinOperator, LinearOperator>(m, "LowdinOperator")
      .def(py::init<const LinearOperator&, const LinearOperator&>(), py::arg("inverse_sqrt"),
           py::arg("operator"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<MetricInverseProduct, LinearOperator>(m, "MetricInverseProduct")
      .def(py::init<const LinearOperator&, const LinearOperator&>(), py::arg("inverse"),
           py::arg("operator"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<HermitianTripleProduct, LinearOperator>(m, "HermitianTripleProduct")
      .def(py::init<const LinearOperator&, const LinearOperator&, const LinearOperator&, Scalar>(),
           py::arg("first"), py::arg("middle"), py::arg("last"),
           py::arg("coefficient") = Scalar(0.5, 0.0),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>());

  py::class_<MetricPolicy>(m, "MetricPolicy")
      .def_property_readonly("is_orthogonal", &MetricPolicy::IsOrthogonal)
      .def_property_readonly("overlap", &MetricPolicy::Overlap,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("inverse_overlap", &MetricPolicy::InverseOverlap,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("inverse_sqrt_overlap", &MetricPolicy::InverseSqrtOverlap,
                             py::return_value_policy::reference_internal);

  py::class_<OrthogonalMetric, MetricPolicy>(m, "OrthogonalMetric")
      .def(py::init<Index>(), py::arg("dim") = 0);

  py::class_<NonOrthogonalMetric, MetricPolicy>(m, "NonOrthogonalMetric")
      .def(py::init<const SparseMatrix&, const LinearOperator&, const LinearOperator&>(),
           py::arg("s"), py::arg("inverse"), py::arg("inverse_sqrt"), py::keep_alive<1, 2>(),
           py::keep_alive<1, 3>(), py::keep_alive<1, 4>());

  py::class_<ExactInverseOverlap, LinearOperator>(m, "ExactInverseOverlap")
      .def(py::init<const SparseMatrix&, Index>(), py::arg("s"), py::arg("max_dimension") = 4096);

  py::class_<ExactInverseSqrtOverlap, LinearOperator>(m, "ExactInverseSqrtOverlap")
      .def(py::init<const SparseMatrix&, Index>(), py::arg("s"), py::arg("max_dimension") = 4096);

  py::class_<IschoInverseOverlap, LinearOperator>(m, "IschoInverseOverlap")
      .def(py::init<const SparseMatrix&, double, int>(), py::arg("s"),
           py::arg("tolerance") = 1e-13, py::arg("max_iterations") = 200, py::keep_alive<1, 2>())
      .def_property_readonly("last_iterations", &IschoInverseOverlap::LastIterations)
      .def_property_readonly("last_error", &IschoInverseOverlap::LastError);

  py::class_<ChebyshevInverseOverlap, LinearOperator>(m, "ChebyshevInverseOverlap")
      .def(py::init<const SparseMatrix&, SpectralBounds, int, double>(), py::arg("s"),
           py::arg("bounds_of_s"), py::arg("order"), py::arg("coeff_tol") = 1e-12,
           py::keep_alive<1, 2>())
      .def_property_readonly("effective_order", &ChebyshevInverseOverlap::EffectiveOrder);

  py::class_<ChebyshevInverseSqrtOverlap, LinearOperator>(m, "ChebyshevInverseSqrtOverlap")
      .def(py::init<const SparseMatrix&, SpectralBounds, int, double>(), py::arg("s"),
           py::arg("bounds_of_s"), py::arg("order"), py::arg("coeff_tol") = 1e-12,
           py::keep_alive<1, 2>())
      .def_property_readonly("effective_order", &ChebyshevInverseSqrtOverlap::EffectiveOrder);

  py::class_<BlockDiagonalOperator, LinearOperator>(m, "BlockDiagonalOperator")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style | py::array::forcecast>& arr) {
             if (arr.ndim() != 3 || arr.shape(1) != arr.shape(2)) {
               throw std::invalid_argument(
                   "BlockDiagonalOperator: blocks array must have shape (N_k, W, W)");
             }
             const Index nk = static_cast<Index>(arr.shape(0));
             const Index w = static_cast<Index>(arr.shape(1));
             std::vector<Scalar> vec(arr.data(), arr.data() + arr.size());
             return BlockDiagonalOperator(nk, w, std::move(vec));
           }),
           py::arg("blocks"))
      .def(py::init([](Index numBlocks, Index blockSize,
                       const py::array_t<Scalar, py::array::c_style | py::array::forcecast>& arr) {
             std::vector<Scalar> vec(arr.data(), arr.data() + arr.size());
             return BlockDiagonalOperator(numBlocks, blockSize, std::move(vec));
           }),
           py::arg("num_blocks"), py::arg("block_size"), py::arg("blocks"))
      .def(py::init<Index, Index, std::vector<Scalar>>(), py::arg("num_blocks"),
           py::arg("block_size"), py::arg("blocks"))
      .def_property_readonly("num_blocks", &BlockDiagonalOperator::NumBlocks)
      .def_property_readonly("block_size", &BlockDiagonalOperator::BlockSize)
      .def("adjoint", &BlockDiagonalOperator::Adjoint);

  py::enum_<BlochDirection>(m, "BlochDirection")
      .value("RealToBloch", BlochDirection::RealToBloch)
      .value("BlochToReal", BlochDirection::BlochToReal);

  py::class_<BlochTransform, LinearOperator>(m, "BlochTransform")
      .def(py::init([](Index n1, Index n2, Index n3, Index numOrbitals,
                       py::object atomPhasesObj, BlochDirection direction) {
             std::vector<Scalar> phases;
             if (!atomPhasesObj.is_none()) {
               auto arr = py::array_t<Scalar, py::array::c_style | py::array::forcecast>::ensure(
                   atomPhasesObj);
               if (!arr) {
                 throw std::invalid_argument(
                     "BlochTransform: atom_phases must be a complex128 array");
               }
               phases.assign(arr.data(), arr.data() + arr.size());
             }
             return BlochTransform(n1, n2, n3, numOrbitals, std::move(phases), direction);
           }),
           py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("num_orbitals"),
           py::arg("atom_phases") = py::none(),
           py::arg("direction") = BlochDirection::RealToBloch)
      .def_property_readonly("n1", &BlochTransform::N1)
      .def_property_readonly("n2", &BlochTransform::N2)
      .def_property_readonly("n3", &BlochTransform::N3)
      .def_property_readonly("num_cells", &BlochTransform::NumCells)
      .def_property_readonly("num_orbitals", &BlochTransform::NumOrbitals)
      .def_property_readonly("direction", &BlochTransform::GetDirection)
      .def("adjoint", &BlochTransform::Adjoint);

  py::class_<StateGenerator>(m, "StateGenerator")
      .def_property_readonly("dimension", &StateGenerator::Dimension)
      .def_property_readonly("number_of_states", &StateGenerator::NumberOfStates)
      .def_property_readonly("label", &StateGenerator::Label)
      .def("reset", &StateGenerator::Reset)
      .def("next", [](StateGenerator& states) -> py::object {
        Vector state;
        if (!states.Next(state)) {
          return py::none();
        }
        py::array_t<Scalar> result(state.size());
        std::copy(state.begin(), state.end(), result.mutable_data());
        return result;
      }, "Return a copy of the next state, or None; reset before replaying.");

  py::class_<RandomPhaseStates, StateGenerator>(m, "RandomPhaseStates")
      .def(py::init<Index, Index, unsigned>(), py::arg("dim"), py::arg("num_states"),
           py::arg("seed") = 12345u);

  py::class_<LocalStates, StateGenerator>(m, "LocalStates")
      .def(py::init([](Index dim, const py::array_t<Index, py::array::c_style>& sites) {
             if (sites.ndim() != 1) {
               throw std::invalid_argument("LocalStates: sites must be a 1D int64 array");
             }
             std::vector<Index> siteVec(sites.data(), sites.data() + sites.size());
             return LocalStates(dim, std::move(siteVec));
           }),
           py::arg("dim"), py::arg("sites"))
      .def(py::init<Index, std::vector<Index>>(), py::arg("dim"), py::arg("sites"));

  py::class_<BlockLocalStates, StateGenerator>(m, "BlockLocalStates")
      .def(py::init<Index, Index>(), py::arg("num_blocks"), py::arg("block_size"))
      .def(py::init<Index, Index, Index, Index>(), py::arg("n1"), py::arg("n2"), py::arg("n3"),
           py::arg("block_size"))
      .def_static("block_index", &BlockLocalStates::BlockIndex, py::arg("m1"), py::arg("m2"),
                  py::arg("m3"), py::arg("n1"), py::arg("n2"));

  py::class_<TransformedStates, StateGenerator>(m, "TransformedStates")
      .def(py::init<const LinearOperator&, StateGenerator&>(), py::arg("transform"),
           py::arg("inner"), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

  py::class_<UserStates, StateGenerator>(m, "UserStates")
      .def(py::init([](const std::vector<py::array_t<Scalar, py::array::c_style>>& arrays) {
             std::vector<Vector> states;
             states.reserve(arrays.size());
             for (const auto& arr : arrays) {
               if (arr.ndim() != 1) {
                 throw std::invalid_argument("UserStates: each state must be a 1D complex128 array");
               }
               states.emplace_back(arr.data(), arr.data() + arr.size());
             }
             return UserStates(std::move(states));
           }),
           py::arg("states"))
      .def(py::init<std::vector<Vector>>(), py::arg("states"));

  py::class_<SpectralBounds>(m, "SpectralBounds")
      .def(py::init<double, double>(), py::arg("lower") = 0.0, py::arg("upper") = 0.0)
      .def_readwrite("lower", &SpectralBounds::lower)
      .def_readwrite("upper", &SpectralBounds::upper);

  m.def("GershgorinBounds", &GershgorinBounds, py::arg("matrix"), py::arg("pad") = 0.10);
  m.def("LanczosBounds", &LanczosBounds, py::arg("operator"), py::arg("iterations") = 40,
        py::arg("seed") = 12345u);

  py::class_<SpectralWindow>(m, "SpectralWindow")
      .def(py::init<>())
      .def(py::init<double, double>(), py::arg("lower"), py::arg("upper"))
      .def_readwrite("lower", &SpectralWindow::lower)
      .def_readwrite("upper", &SpectralWindow::upper)
      .def_property_readonly("scale", &SpectralWindow::Scale)
      .def_property_readonly("shift", &SpectralWindow::Shift)
      .def_property_readonly("half_width", &SpectralWindow::HalfWidth)
      .def_property_readonly("center", &SpectralWindow::Center);

  py::class_<RescaledOperator, LinearOperator>(m, "RescaledOperator")
      .def(py::init<const SparseMatrix&, const SpectralWindow&>(), py::arg("matrix"),
           py::arg("window"), py::keep_alive<1, 2>())
      .def(py::init<const LinearOperator&, const SpectralWindow&>(), py::arg("operator"),
           py::arg("window"), py::keep_alive<1, 2>())
      .def_property_readonly("window", &RescaledOperator::Window);

  py::class_<GeneralisedRescaledOperator, LinearOperator>(m, "GeneralisedRescaledOperator")
      .def(py::init<const SparseMatrix&, const SparseMatrix&, const LinearOperator&,
                    const SpectralWindow&>(),
           py::arg("h"), py::arg("s"), py::arg("inverse"), py::arg("window"),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
      .def_property_readonly("window", &GeneralisedRescaledOperator::Window);

  py::enum_<Register>(m, "Register")
      .value("Ket", Register::Ket)
      .value("Bra", Register::Bra)
      .value("Aux", Register::Aux);

  py::class_<ApplyOperator>(m, "ApplyOperator")
      .def(py::init<int>(), py::arg("operator_index"))
      .def_readwrite("operator_index", &ApplyOperator::operatorIndex);

  py::class_<Expand>(m, "Expand")
      .def(py::init<int>(), py::arg("num_moments"))
      .def_readwrite("num_moments", &Expand::numMoments);

  py::class_<ApplySeries>(m, "ApplySeries")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style>& coeffs, Scalar scale) {
             if (coeffs.ndim() != 1) {
               throw std::invalid_argument("ApplySeries: coeffs must be a 1D complex128 array");
             }
             return ApplySeries{
                 std::vector<Scalar>(coeffs.data(), coeffs.data() + coeffs.size()), scale};
           }),
           py::arg("coeffs"), py::arg("scale") = Scalar(1.0, 0.0))
      .def(py::init<std::vector<Scalar>, Scalar>(), py::arg("coeffs"),
           py::arg("scale") = Scalar(1.0, 0.0))
      .def_readwrite("coeffs", &ApplySeries::coeffs)
      .def_readwrite("scale", &ApplySeries::scale);

  py::class_<ApplyDerivedSeries>(m, "ApplyDerivedSeries")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style>& coeffs, int operatorIndex,
                       Scalar scale) {
             if (coeffs.ndim() != 1) {
               throw std::invalid_argument(
                   "ApplyDerivedSeries: coeffs must be a 1D complex128 array");
             }
             return ApplyDerivedSeries{
                 std::vector<Scalar>(coeffs.data(), coeffs.data() + coeffs.size()), operatorIndex,
                 scale};
           }),
           py::arg("coeffs"), py::arg("operator_index"), py::arg("scale") = Scalar(1.0, 0.0))
      .def(py::init<std::vector<Scalar>, int, Scalar>(), py::arg("coeffs"),
           py::arg("operator_index"), py::arg("scale") = Scalar(1.0, 0.0))
      .def_readwrite("coeffs", &ApplyDerivedSeries::coeffs)
      .def_readwrite("operator_index", &ApplyDerivedSeries::operatorIndex)
      .def_readwrite("scale", &ApplyDerivedSeries::scale);

  py::class_<MirrorToBra>(m, "MirrorToBra").def(py::init<>());

  py::class_<EvolveRegister>(m, "EvolveRegister")
      .def(py::init([](Register target, const py::array_t<Scalar, py::array::c_style>& coeffs,
                       Scalar scale) {
             if (coeffs.ndim() != 1) {
               throw std::invalid_argument("EvolveRegister: coeffs must be a 1D complex128 array");
             }
             return EvolveRegister{
                 target, std::vector<Scalar>(coeffs.data(), coeffs.data() + coeffs.size()), scale};
           }),
           py::arg("target"), py::arg("coeffs"), py::arg("scale") = Scalar(1.0, 0.0))
      .def(py::init<Register, std::vector<Scalar>, Scalar>(), py::arg("target"), py::arg("coeffs"),
           py::arg("scale") = Scalar(1.0, 0.0))
      .def_readwrite("target", &EvolveRegister::target)
      .def_readwrite("coeffs", &EvolveRegister::coeffs)
      .def_readwrite("scale", &EvolveRegister::scale);

  py::class_<EvolveDerived>(m, "EvolveDerived")
      .def(py::init([](Register target, Register source,
                       const py::array_t<Scalar, py::array::c_style>& coeffs, int operatorIndex,
                       Scalar scale) {
             if (coeffs.ndim() != 1) {
               throw std::invalid_argument("EvolveDerived: coeffs must be a 1D complex128 array");
             }
             return EvolveDerived{
                 target, source,
                 std::vector<Scalar>(coeffs.data(), coeffs.data() + coeffs.size()), operatorIndex,
                 scale};
           }),
           py::arg("target"), py::arg("source"), py::arg("coeffs"), py::arg("operator_index"),
           py::arg("scale") = Scalar(1.0, 0.0))
      .def(py::init<Register, Register, std::vector<Scalar>, int, Scalar>(), py::arg("target"),
           py::arg("source"), py::arg("coeffs"), py::arg("operator_index"),
           py::arg("scale") = Scalar(1.0, 0.0))
      .def_readwrite("target", &EvolveDerived::target)
      .def_readwrite("source", &EvolveDerived::source)
      .def_readwrite("coeffs", &EvolveDerived::coeffs)
      .def_readwrite("operator_index", &EvolveDerived::operatorIndex)
      .def_readwrite("scale", &EvolveDerived::scale);

  py::class_<CopyRegister>(m, "CopyRegister")
      .def(py::init<Register, Register>(), py::arg("from_reg"), py::arg("to_reg"))
      .def_readwrite("from_reg", &CopyRegister::from)
      .def_readwrite("to_reg", &CopyRegister::to);

  py::class_<TimeAxis>(m, "TimeAxis")
      .def(py::init<int, std::vector<AxisStep>>(), py::arg("num_steps"),
           py::arg("evolve") = std::vector<AxisStep>{})
      .def_readwrite("num_steps", &TimeAxis::numSteps)
      .def_readwrite("evolve", &TimeAxis::evolve);

  py::class_<Trace>(m, "Trace")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style>& weights) {
             if (weights.ndim() != 1) {
               throw std::invalid_argument("Trace: weights must be a 1D complex128 array");
             }
             return Trace{std::vector<Scalar>(weights.data(), weights.data() + weights.size())};
           }),
           py::arg("weights"))
      .def(py::init<std::vector<Scalar>>(), py::arg("weights"))
      .def_readwrite("weights", &Trace::weights);

  py::class_<Plan>(m, "Plan")
      .def(py::init<>())
      .def(py::init<std::vector<Instruction>, std::string>(),
           py::arg("code") = std::vector<Instruction>{}, py::arg("word") = "")
      .def_readwrite("code", &Plan::code)
      .def_readwrite("word", &Plan::word);

  py::class_<MomentEngine>(m, "MomentEngine")
      .def(py::init<const LinearOperator&, std::vector<const LinearOperator*>, StateGenerator&,
                    int>(),
           py::arg("rescaled_h"), py::arg("operators"), py::arg("states"), py::arg("batch_size"),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
      // All registered operators and state generators are native C++.
      // The complete recursion/trace runs without a Python callback or GIL.
      .def("run", &MomentEngine::Run, py::arg("plan"),
           py::call_guard<py::gil_scoped_release>());

  py::class_<ChebyshevMoments>(m, "ChebyshevMoments")
      .def(py::init<std::vector<int>>(), py::arg("dims") = std::vector<int>{})
      .def_property_readonly("rank", &ChebyshevMoments::Rank)
      .def_property_readonly("dims", &ChebyshevMoments::Dims)
      .def_property_readonly(
          "data",
          [](py::object selfObj) {
            auto& self = selfObj.cast<ChebyshevMoments&>();
            std::vector<py::ssize_t> shape;
            shape.reserve(self.Dims().size());
            for (int d : self.Dims()) {
              shape.push_back(static_cast<py::ssize_t>(d));
            }
            std::vector<py::ssize_t> strides(shape.size());
            py::ssize_t stride = static_cast<py::ssize_t>(sizeof(Scalar));
            for (py::ssize_t i = static_cast<py::ssize_t>(shape.size()) - 1; i >= 0; --i) {
              strides[static_cast<std::size_t>(i)] = stride;
              stride *= shape[static_cast<std::size_t>(i)];
            }
            return py::array_t<Scalar>(shape, strides, self.Data().data(), selfObj);
          })
      .def_readwrite("window", &ChebyshevMoments::window)
      .def_readwrite("system_size", &ChebyshevMoments::systemSize)
      .def_readwrite("number_of_states", &ChebyshevMoments::numberOfStates)
      .def_readwrite("state_label", &ChebyshevMoments::stateLabel)
      .def_readwrite("word", &ChebyshevMoments::word)
      .def_readonly("build_profile", &ChebyshevMoments::buildProfile)
      .def_readwrite("kpm_weighted", &ChebyshevMoments::kpmWeighted)
      .def_readwrite("symmetrised", &ChebyshevMoments::symmetrised)
      .def_readwrite("time_step", &ChebyshevMoments::timeStep)
      .def_readwrite("number_of_time_steps", &ChebyshevMoments::numberOfTimeSteps);

  m.def("ApplyKpmWeights", &ApplyKpmWeights, py::arg("moments"), py::arg("num_states"),
        py::arg("symmetrise"));

  m.def("ReadCsrText", &ReadCsrText, py::arg("path"), py::arg("legacy_fold") = false);
  m.def("ReadChebmom", &ReadChebmom, py::arg("path"));
  m.def("WriteChebmom", &WriteChebmom, py::arg("path"), py::arg("moments"));
  m.def("WriteLegacyChebmom1D", &WriteLegacyChebmom1D, py::arg("path"), py::arg("moments"));
  m.def("WriteLegacyChebmom2D", &WriteLegacyChebmom2D, py::arg("path"), py::arg("moments"));
  m.def("WriteLegacyChebmomTD", &WriteLegacyChebmomTD, py::arg("path"), py::arg("moments"));

  m.def(
      "ContractRank1",
      [](const py::array_t<double>& moments, const py::array_t<double>& table) {
        if (moments.ndim() != 1) {
          throw std::invalid_argument("ContractRank1: moments must be a 1D float64 array");
        }
        const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
        if (moments.strides(0) % elemSize != 0) {
          throw std::invalid_argument("ContractRank1: unaligned moments stride");
        }
        const Index momentStride = static_cast<Index>(moments.strides(0) / elemSize);
        const Index numMoments = static_cast<Index>(moments.shape(0));
        return ContractRank1Buffer(moments.data(), momentStride, numMoments, table);
      },
      py::arg("moments"), py::arg("table"));

  m.def(
      "ContractRank1",
      [](const ChebyshevMoments& moments, const py::array_t<double>& table) {
        if (moments.Rank() != 1) {
          throw std::invalid_argument("ContractRank1: ChebyshevMoments must have rank 1");
        }
        const Index numMoments = static_cast<Index>(moments.Dims()[0]);
        const auto* realPtr = reinterpret_cast<const double*>(moments.Data().data());
        return ContractRank1Buffer(realPtr, 2, numMoments, table);
      },
      py::arg("moments"), py::arg("table"));

  m.def(
      "ContractRank2Imag",
      [](const py::array_t<Scalar, py::array::c_style>& moments,
         const py::array_t<double, py::array::c_style>& tableA,
         const py::array_t<Scalar, py::array::c_style>& tableB) {
        if (moments.ndim() != 2) {
          throw std::invalid_argument("ContractRank2Imag: moments must be a 2D complex128 array");
        }
        const Index m0 = static_cast<Index>(moments.shape(0));
        const Index m1 = static_cast<Index>(moments.shape(1));
        const Index numPoints =
            ValidateRank2Tables("ContractRank2Imag", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2Imag(moments.data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                          result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "ContractRank2Imag",
      [](const ChebyshevMoments& moments, const py::array_t<double, py::array::c_style>& tableA,
         const py::array_t<Scalar, py::array::c_style>& tableB) {
        if (moments.Rank() != 2) {
          throw std::invalid_argument("ContractRank2Imag: ChebyshevMoments must have rank 2");
        }
        const Index m0 = static_cast<Index>(moments.Dims()[0]);
        const Index m1 = static_cast<Index>(moments.Dims()[1]);
        const Index numPoints =
            ValidateRank2Tables("ContractRank2Imag", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2Imag(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                          result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "ContractRank2Real",
      [](const py::array_t<Scalar, py::array::c_style>& moments,
         const py::array_t<Scalar, py::array::c_style>& tableA,
         const py::array_t<Scalar, py::array::c_style>& tableB) {
        if (moments.ndim() != 2) {
          throw std::invalid_argument("ContractRank2Real: moments must be a 2D complex128 array");
        }
        const Index m0 = static_cast<Index>(moments.shape(0));
        const Index m1 = static_cast<Index>(moments.shape(1));
        const Index numPoints =
            ValidateRank2Tables("ContractRank2Real", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2Real(moments.data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                          result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "ContractRank2Real",
      [](const ChebyshevMoments& moments, const py::array_t<Scalar, py::array::c_style>& tableA,
         const py::array_t<Scalar, py::array::c_style>& tableB) {
        if (moments.Rank() != 2) {
          throw std::invalid_argument("ContractRank2Real: ChebyshevMoments must have rank 2");
        }
        const Index m0 = static_cast<Index>(moments.Dims()[0]);
        const Index m1 = static_cast<Index>(moments.Dims()[1]);
        const Index numPoints =
            ValidateRank2Tables("ContractRank2Real", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2Real(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                          result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "ContractRank2RealReal",
      [](const py::array_t<Scalar, py::array::c_style>& moments,
         const py::array_t<double, py::array::c_style>& tableA,
         const py::array_t<double, py::array::c_style>& tableB) {
        if (moments.ndim() != 2) {
          throw std::invalid_argument(
              "ContractRank2RealReal: moments must be a 2D complex128 array");
        }
        const Index m0 = static_cast<Index>(moments.shape(0));
        const Index m1 = static_cast<Index>(moments.shape(1));
        const Index numPoints =
            ValidateRank2Tables("ContractRank2RealReal", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2RealReal(moments.data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                              result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "ContractRank2RealReal",
      [](const ChebyshevMoments& moments, const py::array_t<double, py::array::c_style>& tableA,
         const py::array_t<double, py::array::c_style>& tableB) {
        if (moments.Rank() != 2) {
          throw std::invalid_argument("ContractRank2RealReal: ChebyshevMoments must have rank 2");
        }
        const Index m0 = static_cast<Index>(moments.Dims()[0]);
        const Index m1 = static_cast<Index>(moments.Dims()[1]);
        const Index numPoints =
            ValidateRank2Tables("ContractRank2RealReal", m0, m1, tableA, tableB);
        py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
        ContractRank2RealReal(moments.Data().data(), m0, m1, tableA.data(), tableB.data(), numPoints,
                              result.mutable_data());
        return result;
      },
      py::arg("moments"), py::arg("table_a"), py::arg("table_b"));

  m.def(
      "WriteCurveText",
      [](const std::string& path, const py::array_t<double>& x, const py::array_t<double>& y) {
        if (x.ndim() != 1 || y.ndim() != 1) {
          throw std::invalid_argument("WriteCurveText: x and y must be 1D float64 arrays");
        }
        if (x.shape(0) != y.shape(0)) {
          throw std::invalid_argument("WriteCurveText: x and y must have the same length");
        }
        const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
        if (x.strides(0) % elemSize != 0 || y.strides(0) % elemSize != 0) {
          throw std::invalid_argument("WriteCurveText: unaligned array strides");
        }
        const Index xStride = static_cast<Index>(x.strides(0) / elemSize);
        const Index yStride = static_cast<Index>(y.strides(0) / elemSize);
        WriteCurveText(path, x.data(), xStride, y.data(), yStride,
                       static_cast<Index>(x.shape(0)));
      },
      py::arg("path"), py::arg("x"), py::arg("y"));

  py::module_ compatMod = m.def_submodule(
      "compat", "Gate-only compatibility kernels reproducing legacy libstdc++ rounding.");
  compatMod.def(
      "legacy_green_tables",
      [](const py::array_t<double, py::array::c_style | py::array::forcecast>& xGrid,
         Index numMoments0, py::object numMoments1Obj) {
        if (xGrid.ndim() != 1) {
          throw std::invalid_argument("legacy_green_tables: x_grid must be a 1D float64 array");
        }
        const Index m1 = numMoments1Obj.is_none() ? numMoments0 : numMoments1Obj.cast<Index>();
        if (numMoments0 <= 0 || m1 <= 0) {
          throw std::invalid_argument("legacy_green_tables: moment counts must be positive");
        }
        const Index numPoints = static_cast<Index>(xGrid.shape(0));
        py::array_t<Scalar> greenR({static_cast<py::ssize_t>(numMoments0),
                                    static_cast<py::ssize_t>(numPoints)});
        py::array_t<Scalar> dgreenR({static_cast<py::ssize_t>(m1),
                                     static_cast<py::ssize_t>(numPoints)});
        compat::ComputeLegacyGreenTables(xGrid.data(), numPoints, numMoments0, m1,
                                         greenR.mutable_data(), dgreenR.mutable_data());
        return std::make_pair(std::move(greenR), std::move(dgreenR));
      },
      py::arg("x_grid"), py::arg("num_moments_0"), py::arg("num_moments_1") = py::none());
  compatMod.def(
      "bessel_series_coefficients",
      [](double x, double tol) {
        std::vector<Scalar> coeffs = compat::ComputeLegacyBesselCoefficients(x, tol);
        py::array_t<Scalar> out(static_cast<py::ssize_t>(coeffs.size()));
        if (!coeffs.empty()) {
          std::copy(coeffs.begin(), coeffs.end(), out.mutable_data());
        }
        return out;
      },
      py::arg("x"), py::arg("tol") = 1e-15);
}
