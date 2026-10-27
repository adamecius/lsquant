# Spin Hall Conductivity from SIESTA

This tutorial shows how to compute the intrinsic spin Hall conductivity of a two-dimensional material starting from a SIESTA calculation.

## Prerequisites

- A converged SIESTA calculation with spin–orbit coupling (files `.fdf`, `.TSHS` or `.HSX`, and `.ion.xml`).
- LSQUANT built with the `sisl` Python package available.

## Example: Bi(111) bilayer

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run.json --out /tmp/bi_shc
python3 -m lsquant inspect /tmp/bi_shc
```

The calculation evaluates

```
Tr[JXSZ delta(E-H) VY dG-(E)/dE]
```

where `JXSZ` is the spin current operator. The result is the spin Hall conductivity $\sigma^z_{xy}(E)$.

## Key settings in the JSON file

- `"model": {"kind": "siesta", ...}` loads the SIESTA Hamiltonian and overlap.
- `"basis": {"name": "eigen"}` is used for the reference calculation on the unit cell.
- `"metric": {"policy": "nonorthogonal", "route": "lowdin"}` handles the non-orthogonal basis.
- Velocity operators can be chosen as `bare`, `berry_connection`, or `covariant`.

## Expected result

A spin Hall conductivity that is roughly constant (quantized or large) inside the gap or near the band edges, depending on the material.

## Next

- Change the energy window and number of moments.
- Try the PdSe₂ monolayer example in `examples/08_spin_hall_pdse2_monolayer/`.
