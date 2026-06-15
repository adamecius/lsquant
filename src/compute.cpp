#include "compute.hpp"

#include <iostream>
#include <string>
#include <array>

#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"

namespace lsquant
{
	int compute_noneq(const std::string& label, const std::string& op_right,
	                  const std::string& op_left, int num_moments, const std::string& state_file)
	{
		chebyshev::Moments2D chebMoms(num_moments, num_moments);

		SparseMatrixType OP[3];
		OP[0].SetID("HAM");
		OP[1].SetID(op_right);
		OP[2].SetID(op_left);

		SparseMatrixBuilder builder;
		std::array<double,2> spectral_bounds;
		for (int i = 0; i < 3; i++)
		{
			const std::string input = "operators/" + label + "." + OP[i].ID() + ".CSR";
			builder.setSparseMatrix(&OP[i]);
			builder.BuildOPFromCSRFile(input);
			if (i == 0) // Hamiltonian: bounds from the descriptor if present, else BOUNDS/Gershgorin
				spectral_bounds = chebyshev::utility::SpectralBounds(OP[0], "operators/" + label + ".HAM.desc");
		}

		chebMoms.SystemLabel(label);
		chebMoms.BandWidth ( (spectral_bounds[1] - spectral_bounds[0]) * 1.0);
		chebMoms.BandCenter( (spectral_bounds[1] + spectral_bounds[0]) * 0.5);
		chebMoms.SetAndRescaleHamiltonian(OP[0]);
		chebMoms.Print();

		qstates::generator gen;
		if (!state_file.empty())
			gen = qstates::LoadStateFile(state_file);

		chebyshev::CorrelationExpansionMoments(OP[1], OP[2], chebMoms, gen);

		const std::string outputfilename =
			"NonEqOp" + op_right + "-" + op_left + label + "KPM_M" +
			std::to_string(num_moments) + "x" + std::to_string(num_moments) +
			"_state" + gen.StateLabel() + ".chebmom2D";
		chebMoms.saveIn(outputfilename);
		std::cout << "Saved moments in " << outputfilename << std::endl;
		return 0;
	}
}
