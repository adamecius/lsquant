# How to Write a JSON Run File

Every calculation in LSQUANT is defined by a single JSON file. This page shows the minimal structure and the most common options.

## Minimal example (density of states)

```json
{
  "model": {
    "kind": "csr",
    "dir": "examples/01_dos_1d_chain/operators",
    "label": "chain",
    "operators": []
  },
  "basis": { "name": "real" },
  "word": "Tr[delta(E-H)]",
  "expansion": {
    "moments": 128,
    "window": [-2.5, 2.5],
    "trace": "exact"
  },
  "reconstruction": {
    "kernel": "jackson",
    "broadening_meV": 10.0,
    "formulas": ["dos"]
  },
  "output": { "dir": "." }
}
```

## Common fields

| Field | Meaning | Typical values |
|-------|---------|----------------|
| `model.kind` | How the Hamiltonian is loaded | `"csr"`, `"siesta"` |
| `basis.name` | Representation | `"real"`, `"bloch"`, `"eigen"` |
| `word` | Physical quantity to compute | See [Observable Catalogue](../reference/observable-words.md) |
| `expansion.moments` | Number of Chebyshev moments | 64–512 |
| `expansion.window` | Energy window [min, max] | Must cover the spectrum |
| `reconstruction.broadening_meV` | Broadening (meV) | 5–50 |
| `reconstruction.formulas` | What to reconstruct | `["dos"]`, `["bastin"]`, `["greenwood"]`, … |

## Run it

```bash
python3 -m lsquant run my_calculation.json --out /tmp/result
python3 -m lsquant inspect /tmp/result
```

## Tips

- Start with a small number of moments and a larger broadening, then refine.
- For SIESTA models set `"kind": "siesta"` and point to the `.fdf` / `.TSHS` files.
- Use `"basis": {"name": "eigen"}` on small systems when you want an exact reference calculation.
