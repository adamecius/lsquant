#ifndef LSQUANT_METRIC_POLICY
#define LSQUANT_METRIC_POLICY

#include <vector>
#include <complex>
#include <stdexcept>

// Phase 6 SCAFFOLDING (interface only -- not wired into the core yet).
//
// Non-orthogonality becomes a POLICY rather than a forked code path (the current `_nonOrth`
// files). An orthogonal basis uses the identity metric (today's behaviour). A non-orthogonal
// basis carries an overlap matrix S, and the Chebyshev recursion must work against S^{-1}
// (solve S y = x each step; iterative CG when S is matrix-free).
//
// STATUS: OrthogonalMetric is the real, trivial default (identity). NonOrthogonalMetric is a
// PLACEHOLDER (throws) -- the S-solve lands in Phase 6, with the CG tolerance as that phase's
// only human checkpoint, guarded by kpm-nonOrth-precision-test.

namespace lsquant
{
	class MetricPolicy
	{
	public:
		typedef std::complex<double>   scalar;
		typedef std::vector<scalar>    vector_t;

		virtual ~MetricPolicy() {}

		// Apply S^{-1} in place (orthogonal basis: a no-op).
		virtual void apply_inverse_metric(vector_t& v) const = 0;

		virtual bool is_orthogonal() const = 0;
	};

	// Orthogonal basis: S = I, so S^{-1} is the identity. The real default behaviour.
	class OrthogonalMetric : public MetricPolicy
	{
	public:
		void apply_inverse_metric(vector_t& /*v*/) const override { /* identity: nothing to do */ }
		bool is_orthogonal() const override { return true; }
	};

	// PLACEHOLDER (Phase 6): non-orthogonal basis -- solve S y = v (CG when S is matrix-free).
	// Not implemented yet -- throws so an accidental use cannot silently behave as orthogonal.
	class NonOrthogonalMetric : public MetricPolicy
	{
	public:
		void apply_inverse_metric(vector_t& /*v*/) const override
		{ throw std::logic_error("NonOrthogonalMetric: S^{-1} solve not implemented yet (Phase 6)"); }
		bool is_orthogonal() const override { return false; }
	};
}

#endif // LSQUANT_METRIC_POLICY
