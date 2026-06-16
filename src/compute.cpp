#include "compute.hpp"

#include <iostream>
#include <string>
#include <array>
#include <sstream>

#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"
#include "util/report.hpp"

namespace lsquant
{
	int compute_noneq(const std::string& label, const std::string& op_right,
	                  const std::string& op_left, int num_moments, const std::string& state_file,
	                  std::string* out_filename)
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

		{
			LSQ_SCOPED("compute_noneq/CorrelationExpansion");
			chebyshev::CorrelationExpansionMoments(OP[1], OP[2], chebMoms, gen);
		}

		const std::string outputfilename =
			"NonEqOp" + op_right + "-" + op_left + label + "KPM_M" +
			std::to_string(num_moments) + "x" + std::to_string(num_moments) +
			"_state" + gen.StateLabel() + ".chebmom2D";
		chebMoms.saveIn(outputfilename);
		if (out_filename) *out_filename = outputfilename;
		LSQ_LOG_INFO("saved moments in " << outputfilename);
		return 0;
	}

	int compute_spectral(const std::string& label, const std::string& op,
	                     int num_moments, const std::string& state_file,
	                     std::string* out_filename)
	{
		chebyshev::Moments1D chebMoms(num_moments);

		const int NOPS = 2;
		SparseMatrixType OP[NOPS];
		OP[0].SetID("HAM");
		OP[1].SetID(op);

		SparseMatrixBuilder builder;
		std::array<double,2> spectral_bounds;
		for (int i = 0; i < NOPS; i++)
			if (OP[i].isIdentity())
				std::cout << "The operator " << OP[i].ID() << " is treated as the identity" << std::endl;
			else
			{
				const std::string input = "operators/" + label + "." + OP[i].ID() + ".CSR";
				builder.setSparseMatrix(&OP[i]);
				builder.BuildOPFromCSRFile(input);
				if (i == 0) // Hamiltonian: obtain the energy bounds automatically
					spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
			}

		chebMoms.SystemLabel(label);
		chebMoms.BandWidth ( (spectral_bounds[1] - spectral_bounds[0]) * 1.0);
		chebMoms.BandCenter( (spectral_bounds[1] + spectral_bounds[0]) * 0.5);
		chebMoms.SetAndRescaleHamiltonian(OP[0]);
		chebMoms.Print();

		qstates::generator gen;
		if (!state_file.empty())
			gen = qstates::LoadStateFile(state_file);

		{
			LSQ_SCOPED("compute_spectral/SpectralMoments");
			chebyshev::SpectralMoments(OP[1], chebMoms, gen);
		}

		const std::string prefix = "SpectralOp" + OP[1].ID();
		const std::string outputfilename =
			prefix + label + "KPM_M" + std::to_string(num_moments) +
			"_state" + gen.StateLabel() + ".chebmom1D";
		LSQ_LOG_INFO("saving the moments in " << outputfilename);
		chebMoms.saveIn(outputfilename);
		if (out_filename) *out_filename = outputfilename;
		return 0;
	}

	int compute_msd(const std::string& label, const std::string& vop,
	                int num_moments, int num_times, double tmax, const std::string& state_file,
	                std::string* out_filename)
	{
		chebyshev::MomentsTD chebMoms(num_moments, num_times);

		SparseMatrixType OP[2];
		OP[0].SetID("HAM");
		OP[1].SetID(vop);

		SparseMatrixBuilder builder;
		std::array<double,2> spectral_bounds;
		for (int i = 0; i < 2; i++)
			if (OP[i].isIdentity())
				std::cout << "The operator " << OP[i].ID() << " is treated as the identity" << std::endl;
			else
			{
				const std::string input = "operators/" + label + "." + OP[i].ID() + ".CSR";
				builder.setSparseMatrix(&OP[i]);
				builder.BuildOPFromCSRFile(input);
				if (i == 0) // Hamiltonian: obtain the energy bounds automatically
					spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
			}

		chebMoms.SystemLabel(label);
		chebMoms.BandWidth ( (spectral_bounds[1] - spectral_bounds[0]) * 1.0);
		chebMoms.BandCenter( (spectral_bounds[1] + spectral_bounds[0]) * 0.5);
		chebMoms.TimeDiff( tmax / (num_times - 1) );
		chebMoms.SetAndRescaleHamiltonian(OP[0]);
		chebMoms.SetAndRescaleConmuter(OP[1]);
		chebMoms.Print();

		qstates::generator gen;
		if (!state_file.empty())
			gen = qstates::LoadStateFile(state_file);

		{
			LSQ_SCOPED("compute_msd/MeanSquareDisplacement");
			chebyshev::MeanSquareDisplacement(chebMoms, gen);
		}

		// Default-format tmax (20.0 -> "20", 0.5 -> "0.5") so the filename token matches the legacy
		// argv string and the .chebmomTD name is preserved byte-for-byte.
		std::ostringstream tmaxtok; tmaxtok << tmax;
		std::string prefix = "Correlation" + OP[1].ID();
		if (OP[1].isIdentity())
			prefix = "Correlations" + OP[1].ID();
		const std::string outputfilename =
			prefix + label + "KPM_M" + std::to_string(num_moments) +
			"_state" + tmaxtok.str() + gen.StateLabel() + ".chebmomTD";
		LSQ_LOG_INFO("saving convergence data in " << outputfilename);
		chebMoms.saveIn(outputfilename);
		if (out_filename) *out_filename = outputfilename;
		return 0;
	}
}
