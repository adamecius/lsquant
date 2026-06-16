#ifndef LSQUANT_HAMILTONIAN_OP
#define LSQUANT_HAMILTONIAN_OP

#include <vector>
#include <complex>
#include <cstddef>
#include <stdexcept>

#include "sparse_matrix.hpp"   // SparseMatrixType (the live adapter wraps it)

// Phase 6: the KPM recursion's contact with H is this interface. The stored sparse matrix is one
// implementation (SparseHamiltonianOp); an on-the-fly ("matrix-free") analytic model would be a
// second, with no change to the recursion.
//
// DESIGN RULE: virtual dispatch is at "apply to a whole vector" granularity (multiply), NEVER
// per matrix element -- so the hot Chebyshev recursion pays one virtual call per matvec. The
// descriptor/sidecar (identity, units, bounds) never enters here; this is numeric buffers only.
//
// STATUS: SparseHamiltonianOp is the live implementation (Phase 6A). The matrix-free placeholder
// below throws until implemented (it must never be silently a no-op).

namespace lsquant
{
	class HamiltonianOp
	{
	public:
		typedef std::complex<double>   scalar;
		typedef std::vector<scalar>    vector_t;

		virtual ~HamiltonianOp() {}

		// Hilbert-space dimension (rank of H).
		virtual std::size_t dim() const = 0;

		// y = a * H * x + b * y  -- the recursion's ONLY contact with H (matches
		// SparseMatrixType::Multiply(a,x,b,y) so the sparse adapter is a trivial forward).
		virtual void multiply(scalar a, const vector_t& x, scalar b, vector_t& y) const = 0;

		// Convenience: y = H * x  (forwards to the 4-arg form -- one virtual dispatch, mirrors
		// SparseMatrixType's 2-arg Multiply so the 2-arg call sites convert with no arg changes).
		void multiply(const vector_t& x, vector_t& y) const
		{ this->multiply(scalar(1.0, 0.0), x, scalar(0.0, 0.0), y); }

		// True if H is already rescaled into the Chebyshev domain [-1,1] (H~). The owner applies
		// the (H - BandCenter)/HalfWidth rescale; bounds come from the operator descriptor.
		virtual bool is_rescaled() const { return false; }
	};

	// The live HamiltonianOp: a trivial forward to a stored SparseMatrixType (already rescaled in
	// place by the owning Moments). Non-owning -- the Moments keeps the matrix alive.
	class SparseHamiltonianOp : public HamiltonianOp
	{
		SparseMatrixType* m_;          // already rescaled in place by the owner
		bool rescaled_;
	public:
		explicit SparseHamiltonianOp(SparseMatrixType& m, bool rescaled)
			: m_(&m), rescaled_(rescaled) {}

		std::size_t dim() const override { return m_->rank(); }

		void multiply(scalar a, const vector_t& x, scalar b, vector_t& y) const override
		{ m_->Multiply(a, x, b, y); }                  // forwards to sparse_matrix.hpp vector_t& overload

		bool is_rescaled() const override { return rescaled_; }

		using HamiltonianOp::multiply;   // keep the inherited 2-arg convenience visible
	};

	// PLACEHOLDER (Phase 6): build H*x on the fly from an analytic tight-binding model
	// (the deferred _hr.dat consumer), with no stored matrix. Not implemented yet -- throws rather
	// than pretending to work, so an accidental use is loud, not silent.
	class MatrixFreeHamiltonianOp : public HamiltonianOp
	{
	public:
		std::size_t dim() const override
		{ throw std::logic_error("MatrixFreeHamiltonianOp: matrix-free apply not implemented yet (Phase 6)"); }

		void multiply(scalar, const vector_t&, scalar, vector_t&) const override
		{ throw std::logic_error("MatrixFreeHamiltonianOp: matrix-free apply not implemented yet (Phase 6)"); }
	};
}

#endif // LSQUANT_HAMILTONIAN_OP
