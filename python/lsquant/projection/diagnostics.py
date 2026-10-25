"""Window-dependent exact word errors; no KPM execution or acceptance tuning."""
from .bases import ProjectedBlochBasis
from .spectral import word_error


def compare_windows(source, words, windows, energies, *, cell_area, broadening=0.01,
                    kind="shift", degeneracy_tolerance=1e-6):
    """Report each supported word's oracle error as a function of window.

    Reference matrices and energy grid are held fixed across the sweep.
    Static operators, including the spin current, are consistently projected
    before every word evaluation. Absolute errors remain meaningful when the
    reference word vanishes; relative errors then have value None.
    """
    from .runfile import exact_word_curve
    options = dict(cell_area=cell_area, broadening=broadening,
                   degeneracy_tolerance=degeneracy_tolerance)
    reference = {word: exact_word_curve(source, word, energies, **options)[0] for word in words}
    rows = []
    for window in windows:
        projected = ProjectedBlochBasis(source, dict(kind=kind, window=window))
        errors = {}
        for word in words:
            curve, _ = exact_word_curve(projected, word, energies, **options)
            errors[word] = word_error(reference[word], curve)
        rows.append(dict(window=list(window), kept_bands_per_k=projected.projection.kept.sum(axis=1).tolist(),
                         spectral_bounds=[float(projected.eigenvalues.min()), float(projected.eigenvalues.max())],
                         word_errors=errors))
    return rows
