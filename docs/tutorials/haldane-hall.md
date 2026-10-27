# Haldane Model — Hall Conductivity

This tutorial computes the quantized Hall conductivity of the Haldane model.

The Haldane model is a tight-binding model on the honeycomb lattice that realizes a Chern insulator. In the topological phase the Hall conductivity is quantized:

$$
\sigma_{xy} = \frac{e^2}{h}
$$

## Run the calculation

```bash
python3 -m lsquant run examples/03_haldane_hall/run_haldane.json --out /tmp/haldane
python3 -m lsquant inspect /tmp/haldane
python3 examples/03_haldane_hall/plot.py /tmp/haldane
```

The output shows the Hall conductivity as a function of Fermi energy. You should observe a clear quantized plateau at $e^2/h$ inside the gap, together with the decomposition into Fermi-sea and Fermi-surface contributions.

## What the calculation does

The observable word used is

```
Tr[VX delta(E-H) VY dG-(E)/dE]
```

This is the Kubo–Bastin formula for the Hall conductivity. LSQUANT evaluates it with the Kernel Polynomial Method on the real-space representation of the Haldane Hamiltonian.

## Next

- Change the broadening or the number of moments in the JSON file to see the effect on the plateau.
- Compare with the longitudinal conductivity computed in the same run.
