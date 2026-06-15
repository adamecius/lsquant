#ifndef LSQUANT_HAMILTONIAN_OP
#define LSQUANT_HAMILTONIAN_OP

#include <vector>
#include <complex>
#include <cstddef>
#include <stdexcept>

// Phase 6 SCAFFOLDING (interface only -- the core is not migrated onto it yet).
//
// Today the KPM core takes a concrete SparseMatrixType*. Phase 6 will have it consume this
// HamiltonianOp interface instead, so the stored-sparse matrix and the on-the-fly ("matrix-free")
// analytic models (Graphene_NNN / Device) become two implementations of one contract rather than
// forked solver/moments/vectors files.
//
// DESIGN RULE: virtual dispatch is at "apply to a whole vector" granularity (multiply), NEVER
// per matrix element -- so the hot Chebyshev recursion pays one virtual call per matvec. The
// descriptor/sidecar (identity, units, bounds) never enters here; this is numeric buffers only.
//
// STATUS: the sparse adapter and the matrix-free implementation are deferred to Phase 6. The
// matrix-free placeholder below exists so callers can compile against the final interface and so
// the boundary is explicit; it throws until implemented (it must never be silently a no-op).

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

		// True if H is already rescaled into the Chebyshev domain [-1,1] (H~). The owner applies
		// the (H - BandCenter)/HalfWidth rescale; bounds come from the operator descriptor.
		virtual bool is_rescaled() const { return false; }
	};

	// PLACEHOLDER (Phase 6): build H*x on the fly from an analytic tight-binding model
	// (Graphene_NNN / Device), with no stored matrix. Not implemented yet -- throws rather than
	// pretending to work, so an accidental use is loud, not silent.
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
