# LSQUANT Documentation Structure (Diátaxis)

This document outlines the recommended structure for LSQUANT documentation.

## tutorials/
Learning-oriented content. Hands-on lessons that build skills.

- **getting-started.md**: Build LSQUANT and run your first DOS calculation on a 1D chain.
- **haldane-hall.md**: Compute the quantized Hall conductivity of the Haldane model.
- **spin-hall-siesta.md**: Load a SIESTA model and compute intrinsic spin Hall conductivity.

## how-to/
Goal-oriented guides for specific tasks.

- **run-json-file.md**: How to write and execute a JSON run file.
- **reconstruct-curve.md**: How to reconstruct physical quantities from moments.
- **use-nonorthogonal.md**: How to handle non-orthogonal bases (Löwdin / generalised).
- **performance-tuning.md**: How to choose golden vs native builds and NUMA settings.

## reference/
Information-oriented technical reference.

- **observable-words.md**: Catalogue of all supported observable word shapes.
- **json-schema.md**: Complete JSON run-file schema.
- **cli.md**: Command-line interface reference.
- **python-api.md**: Python API reference.
- **bases-and-metrics.md**: Bases (real, bloch, eigen) and metric policies.

## explanation/
Understanding-oriented conceptual content.

- **unifying-idea.md**: Why observable words and the MomentEngine architecture.
- **kpm-overview.md**: Brief introduction to the Kernel Polynomial Method.
- **what-is-new.md**: Architecture changes relative to previous versions.

## Root
- **README.md**: Short overview + links to the four sections above.
- **docs/index.md**: Documentation home page.
