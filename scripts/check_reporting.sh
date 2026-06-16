#!/usr/bin/env bash
# Pre-PR reporting self-check (advisory) -- see CLAUDE.md "Reporting requirement".
#
# Flags, in C++ changed vs the merge base (or a given ref), patterns that mean a
# new implementation skipped the reporting layer:
#   * raw std::cout / std::cerr / printf  -> should be LSQ_LOG_* (util/log.hpp)
#   * bare `throw std::runtime_error`/exit -> should be LSQ_THROW(lsquant::*Error)
#
# Scope: added/modified lines only, under src/ and include/, EXCLUDING legacy/ and
# the reporting layer's own implementation. Advisory: exit 0 always; prints a
# report so the author (human or agent) can justify or fix before the PR.
#
# Usage: scripts/check_reporting.sh [base-ref]   (default: origin/main, else HEAD~1)
set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

base="${1:-}"
if [ -z "$base" ]; then
  if git rev-parse --verify -q origin/main >/dev/null; then base="origin/main";
  else base="HEAD~1"; fi
fi
git rev-parse --verify -q "$base" >/dev/null || { echo "check_reporting: base '$base' not found; skipping"; exit 0; }

# Added lines (leading '+', not the +++ header) in tracked C++ outside exclusions.
# Two-dot diff vs base captures both committed and working-tree changes since base.
diff="$(git diff "$base" -- 'src/*.cpp' 'src/*.hpp' 'include/*.hpp' \
        ':(exclude)legacy/**' ':(exclude)src/util/**' ':(exclude)include/util/**' 2>/dev/null)"

added="$(printf '%s\n' "$diff" | grep -E '^\+' | grep -vE '^\+\+\+')"

flag() { printf '%s\n' "$added" | grep -nE "$1" || true; }

io_hits="$(flag '(std::c(out|err)|[^a-zA-Z_]printf)\b')"
throw_hits="$(flag 'throw[[:space:]]+std::(runtime_error|logic_error|exception)|[^a-zA-Z_]exit[[:space:]]*\(')"

echo "== reporting self-check (base: $base) =="
status=0
if [ -n "$io_hits" ]; then
  echo "[advisory] raw I/O in new/changed code -- prefer LSQ_LOG_* (util/log.hpp):"
  printf '%s\n' "$io_hits" | sed 's/^/    /'
  status=1
fi
if [ -n "$throw_hits" ]; then
  echo "[advisory] bare throw/exit in new/changed code -- prefer LSQ_THROW(lsquant::*Error):"
  printf '%s\n' "$throw_hits" | sed 's/^/    /'
  status=1
fi
if [ "$status" -eq 0 ]; then
  echo "OK: no raw I/O or bare throws introduced outside the reporting layer."
else
  echo
  echo "Review the lines above against CLAUDE.md before opening the PR."
fi
exit 0
