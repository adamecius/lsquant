// Standalone `lsquant_inspect` (kept for compatibility; the `lsquant inspect` subcommand shares
// the same lsquant::inspect_path). Phase 7.
#include <iostream>
#include <string>
#include "inspect.hpp"

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: " << argv[0] << " <run.json | operator.desc>\n"; return 2; }
    return lsquant::inspect_path(argv[1]);
}
