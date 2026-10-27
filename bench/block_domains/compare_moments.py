"""Compare moment files written by block_kpm (complex128, raw)."""
import sys

import numpy as np

a = np.fromfile(sys.argv[1], dtype=np.complex128)
b = np.fromfile(sys.argv[2], dtype=np.complex128)
print(f"max_abs_diff_over_max_abs={np.max(np.abs(a - b)) / np.max(np.abs(a)):.3e}")
