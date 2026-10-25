/**
 * @file chebmom_file.cpp
 * @brief Implementation of `.chebmom`, `.chebmom1D`, and `.chebmom2D` file I/O.
 */
#include "io/chebmom_file.hpp"

#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace lsquant {

namespace {
// Item 3: tag native on the existing first header line; golden bytes stay unchanged.
// Tester: round-trip all four layouts across profiles and reject mixed-profile comparisons.
void WriteProfileTag(std::ostream& out, const ChebyshevMoments& moments) {
  if (moments.buildProfile != "golden") {
    out << " # build_profile=" << moments.buildProfile;
  }
}

std::string ReadProfileTag(const std::string& headerTail) {
  const std::string marker = "# build_profile=";
  const auto pos = headerTail.find(marker);
  if (pos == std::string::npos) return "golden";
  std::istringstream value(headerTail.substr(pos + marker.size()));
  std::string profile;
  value >> profile;
  if (profile != "native" && profile != "golden") {
    throw std::runtime_error("ReadChebmom: unknown build profile: " + profile);
  }
  return profile;
}
}  // namespace

void WriteChebmom(const std::string& filename, const ChebyshevMoments& moments) {
  std::ofstream out(filename);
  if (!out.is_open()) {
    throw std::runtime_error("WriteChebmom: cannot open file for writing: " + filename);
  }
  out.precision(17);
  out << "chebmom 1";
  WriteProfileTag(out, moments);
  out << "\n";
  out << "rank " << moments.Rank() << "\n";
  out << "dims";
  for (int d : moments.Dims()) {
    out << " " << d;
  }
  out << "\n";
  out << "window " << moments.window.lower << " " << moments.window.upper << "\n";
  out << "system_size " << moments.systemSize << "\n";
  const std::string label = moments.stateLabel.empty() ? "default" : moments.stateLabel;
  out << "states " << label << " " << moments.numberOfStates << "\n";
  out << "word " << moments.word << "\n";
  out << "kpm_weighted " << (moments.kpmWeighted ? 1 : 0) << "\n";
  out << "symmetrised " << (moments.symmetrised ? 1 : 0) << "\n";
  out << "data\n";
  for (const auto& val : moments.Data()) {
    out << val.real() << " " << val.imag() << "\n";
  }
  if (!out) {
    throw std::runtime_error("WriteChebmom: error writing file: " + filename);
  }
}

ChebyshevMoments ReadChebmom(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw std::runtime_error("ReadChebmom: cannot open file: " + filename);
  }

  std::string firstToken;
  if (!(in >> firstToken)) {
    throw std::runtime_error("ReadChebmom: empty file: " + filename);
  }

  if (firstToken == "chebmom") {
    int version = 0;
    if (!(in >> version) || version != 1) {
      throw std::runtime_error("ReadChebmom: unsupported chebmom version in: " + filename);
    }
    std::string profileTail;
    std::getline(in, profileTail);

    std::string key;
    int rank = 0;
    if (!(in >> key >> rank) || key != "rank" || rank < 0 || rank > 2) {
      throw std::runtime_error("ReadChebmom: invalid rank line in: " + filename);
    }

    if (!(in >> key) || key != "dims") {
      throw std::runtime_error("ReadChebmom: expected dims line in: " + filename);
    }
    std::vector<int> dims(static_cast<std::size_t>(rank));
    for (int i = 0; i < rank; ++i) {
      if (!(in >> dims[static_cast<std::size_t>(i)])) {
        throw std::runtime_error("ReadChebmom: failed reading dims in: " + filename);
      }
    }

    ChebyshevMoments moments(dims);
    moments.buildProfile = ReadProfileTag(profileTail);

    if (!(in >> key >> moments.window.lower >> moments.window.upper) || key != "window") {
      throw std::runtime_error("ReadChebmom: invalid window line in: " + filename);
    }
    if (!(in >> key >> moments.systemSize) || key != "system_size") {
      throw std::runtime_error("ReadChebmom: invalid system_size line in: " + filename);
    }
    if (!(in >> key >> moments.stateLabel >> moments.numberOfStates) || key != "states") {
      throw std::runtime_error("ReadChebmom: invalid states line in: " + filename);
    }

    std::string line;
    std::getline(in >> std::ws, line);
    if (line.rfind("word", 0) != 0) {
      throw std::runtime_error("ReadChebmom: invalid word line in: " + filename);
    }
    if (line.size() > 5) {
      moments.word = line.substr(5);
    } else {
      moments.word.clear();
    }

    int kpmWeightedFlag = 0;
    if (!(in >> key >> kpmWeightedFlag) || key != "kpm_weighted") {
      throw std::runtime_error("ReadChebmom: invalid kpm_weighted line in: " + filename);
    }
    moments.kpmWeighted = (kpmWeightedFlag != 0);

    int symmetrisedFlag = 0;
    if (!(in >> key >> symmetrisedFlag) || key != "symmetrised") {
      throw std::runtime_error("ReadChebmom: invalid symmetrised line in: " + filename);
    }
    moments.symmetrised = (symmetrisedFlag != 0);

    if (!(in >> key) || key != "data") {
      throw std::runtime_error("ReadChebmom: expected data marker in: " + filename);
    }

    for (auto& elem : moments.Data()) {
      double re = 0.0;
      double im = 0.0;
      if (!(in >> re >> im)) {
        throw std::runtime_error("ReadChebmom: truncated data section in: " + filename);
      }
      elem = Scalar(re, im);
    }
    return moments;
  }

  // Legacy .chebmom1D or .chebmom2D layout
  Index systemSize = 0;
  try {
    systemSize = std::stoll(firstToken);
  } catch (const std::exception&) {
    throw std::runtime_error("ReadChebmom: unrecognized header token in: " + filename);
  }

  double bandWidth = 0.0;
  double bandCenter = 0.0;
  if (!(in >> bandWidth >> bandCenter)) {
    throw std::runtime_error("ReadChebmom: invalid legacy header in: " + filename);
  }

  std::string restOfLine1;
  std::getline(in, restOfLine1);
  std::istringstream line1Tail(restOfLine1);
  int tdMaxTimeStep = 0;
  double tdTimeDiff = 0.0;
  const bool isLegacyTD = static_cast<bool>(line1Tail >> tdMaxTimeStep >> tdTimeDiff);

  std::string dimsLine;
  while (std::getline(in, dimsLine)) {
    if (!dimsLine.empty()) {
      break;
    }
  }
  std::istringstream dimsStream(dimsLine);
  std::vector<int> dims;
  int d = 0;
  while (dimsStream >> d) {
    dims.push_back(d);
  }
  if (dims.empty() || dims.size() > 2) {
    throw std::runtime_error("ReadChebmom: invalid legacy dimensions line in: " + filename);
  }

  ChebyshevMoments moments(dims);
  moments.buildProfile = ReadProfileTag(restOfLine1);
  moments.systemSize = systemSize;
  moments.window.lower = bandCenter - 0.5 * bandWidth;
  moments.window.upper = bandCenter + 0.5 * bandWidth;
  moments.kpmWeighted = true;
  if (isLegacyTD) {
    moments.symmetrised = false;
    moments.numberOfTimeSteps = tdMaxTimeStep;
    moments.timeStep = tdTimeDiff;
  } else {
    moments.symmetrised = (moments.Rank() == 2);
  }

  const std::size_t slashPos = filename.find_last_of("/\\");
  const std::string base =
      (slashPos == std::string::npos) ? filename : filename.substr(slashPos + 1);
  const std::size_t statePos = base.rfind("_state");
  if (statePos != std::string::npos) {
    const std::size_t start = statePos + 6;
    const std::size_t dotPos = base.find('.', start);
    moments.stateLabel = (dotPos == std::string::npos) ? base.substr(start)
                                                       : base.substr(start, dotPos - start);
  }

  for (auto& elem : moments.Data()) {
    double re = 0.0;
    double im = 0.0;
    if (!(in >> re >> im)) {
      throw std::runtime_error("ReadChebmom: truncated legacy moment data in: " + filename);
    }
    elem = Scalar(re, im);
  }
  return moments;
}

void WriteLegacyChebmom1D(const std::string& filename, const ChebyshevMoments& moments) {
  if (moments.Rank() != 1) {
    throw std::invalid_argument("WriteLegacyChebmom1D: requires a rank-1 tensor");
  }
  std::ofstream outputfile(filename);
  if (!outputfile.is_open()) {
    throw std::runtime_error("WriteLegacyChebmom1D: cannot open file: " + filename);
  }

  const double bandWidth = (moments.window.upper - moments.window.lower) * 1.0;
  const double bandCenter = (moments.window.upper + moments.window.lower) * 0.5;

  outputfile.precision(std::numeric_limits<double>::digits10);
  outputfile << moments.systemSize << " " << bandWidth << " " << bandCenter;
  WriteProfileTag(outputfile, moments);
  outputfile << std::endl;
  outputfile << moments.Dims()[0] << std::endl;
  for (const auto& mom : moments.Data()) {
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  }
  if (!outputfile) {
    throw std::runtime_error("WriteLegacyChebmom1D: error writing file: " + filename);
  }
}

void WriteLegacyChebmom2D(const std::string& filename, const ChebyshevMoments& moments) {
  if (moments.Rank() != 2) {
    throw std::invalid_argument("WriteLegacyChebmom2D: requires a rank-2 tensor");
  }
  std::ofstream outputfile(filename);
  if (!outputfile.is_open()) {
    throw std::runtime_error("WriteLegacyChebmom2D: cannot open file: " + filename);
  }

  const double bandWidth = (moments.window.upper - moments.window.lower) * 1.0;
  const double bandCenter = (moments.window.upper + moments.window.lower) * 0.5;

  outputfile.precision(std::numeric_limits<double>::digits10);
  outputfile << moments.systemSize << " " << bandWidth << "  " << bandCenter;
  WriteProfileTag(outputfile, moments);
  outputfile << std::endl;
  for (int x : moments.Dims()) {
    outputfile << x << " ";
  }
  outputfile << std::endl;
  for (const auto& mom : moments.Data()) {
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  }
  if (!outputfile) {
    throw std::runtime_error("WriteLegacyChebmom2D: error writing file: " + filename);
  }
}

void WriteLegacyChebmomTD(const std::string& filename, const ChebyshevMoments& moments) {
  if (moments.Rank() != 2) {
    throw std::invalid_argument("WriteLegacyChebmomTD: requires a rank-2 tensor");
  }
  std::ofstream outputfile(filename);
  if (!outputfile.is_open()) {
    throw std::runtime_error("WriteLegacyChebmomTD: cannot open file: " + filename);
  }

  const double bandWidth = (moments.window.upper - moments.window.lower) * 1.0;
  const double bandCenter = (moments.window.upper + moments.window.lower) * 0.5;
  const int numMoms = moments.Dims()[0];
  const int maxTimeStep =
      (moments.numberOfTimeSteps > 0) ? moments.numberOfTimeSteps : moments.Dims()[1];

  outputfile.precision(std::numeric_limits<double>::digits10);
  outputfile << moments.systemSize << " " << bandWidth << " " << bandCenter << " " << maxTimeStep
             << " " << moments.timeStep << " ";
  WriteProfileTag(outputfile, moments);
  outputfile << std::endl;
  outputfile << numMoms << " " << maxTimeStep << " " << std::endl;
  for (const auto& mom : moments.Data()) {
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  }
  if (!outputfile) {
    throw std::runtime_error("WriteLegacyChebmomTD: error writing file: " + filename);
  }
}

}  // namespace lsquant
