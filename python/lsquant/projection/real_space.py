"""Cell-gauge Fourier/WS interpolation using the unchanged origin Model.

This is a local hopping representation, not a global dense supercell
projector. Every physical matrix is transformed with the same convention;
S is retained from the original model and is NEVER radius-truncated here.
The discontinuous window projector can have long-range Fourier tails.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
from scipy import sparse

from .origin import origin_module
from .spectral import adjoint, physical_blocks


def cell_gauge(blocks, k_points, centres):
    """M_cell=D M_atom D^+, D_ii=exp(+ik.tau_i), Cartesian k."""
    phase = np.exp(1j*(np.asarray(k_points) @ np.asarray(centres).T))
    return phase[..., :, None]*blocks*phase.conj()[..., None, :]


def atom_gauge(blocks, k_points, centres):
    """M_atom=D^+ M_cell D, inverse of cell_gauge."""
    phase = np.exp(1j*(np.asarray(k_points) @ np.asarray(centres).T))
    return phase.conj()[..., :, None]*blocks*phase[..., None, :]


class LocalMatrix:
    """Origin FFT/WS hoppings with explicit Cartesian/reduced gauge conversion.

    The origin expects i,j,l nested in that order (l fastest). LSQUANT
    uses i fastest; sort by integer mesh coordinates before its reshape.
    Reject shifted/irregular/incomplete meshes rather than aliasing them.
    Origin WS images share boundary hopping weights by their degeneracy.
    """

    def __init__(self, blocks, lattice, centres, supercell, k_points, *, radius=None):
        self.lattice = np.asarray(lattice, float)
        self.centres = np.asarray(centres, float)
        supplied = np.asarray(supercell)
        if supplied.shape != (3,) or np.any(supplied < 1) or np.any(supplied != supplied.astype(int)):
            raise ValueError("supercell must contain three positive integers")
        self.supercell = supplied.astype(int)
        self.radius = np.inf if radius is None else float(radius)
        if self.radius < 0 or np.isnan(self.radius):
            raise ValueError("hopping radius must be nonnegative")
        points = np.asarray(k_points, float)
        data = physical_blocks(blocks)
        if points.shape != (int(np.prod(self.supercell)), 3) or len(data) != len(points):
            raise ValueError("FFT needs every point of the declared mesh")
        reduced = points @ self.lattice.T/(2*np.pi)
        coordinates = reduced*self.supercell
        integers = np.rint(coordinates).astype(int) % self.supercell
        if np.max(np.abs(coordinates-np.rint(coordinates))) > 1e-8:
            raise ValueError("FFT requires a Gamma-centred regular mesh")
        flat = np.ravel_multi_index(integers.T, tuple(self.supercell))
        if len(np.unique(flat)) != len(data):
            raise ValueError("FFT mesh contains duplicate points")
        order = np.argsort(flat)
        # Model is reused unchanged, including the FFT sign and WS tie rule.
        self.origin = origin_module().Model(SimpleNamespace(cell=self.lattice), self.centres,
                                             self.supercell, cell_gauge(data, points, self.centres)[order],
                                             reduced[order])
        self.training_points = points
        self.training_blocks = data

    def blocks(self, k_points, *, radius=None):
        """Evaluate Hermitian atom-gauge matrices at arbitrary Cartesian k."""
        points = np.asarray(k_points, float)
        single = points.ndim == 1
        points = np.atleast_2d(points)
        reduced = points @ self.lattice.T/(2*np.pi)
        cutoff = self.radius if radius is None else radius
        cells = np.array([self.origin.Hk(k, cutoff) for k in reduced])
        blocks = atom_gauge(cells, points, self.centres)
        return blocks[0] if single else blocks

    def receipt(self):
        """Training-grid reconstruction and radius-only errors, Frobenius norms.

        An off-grid error against the per-k exact projection is separate:
        training-grid interpolation cannot certify between-grid behaviour.
        """
        full = self.blocks(self.training_points, radius=np.inf)
        cut = self.blocks(self.training_points)
        difference = cut-full
        return dict(radius_Angstrom=(None if np.isinf(self.radius) else self.radius),
                    fft_mesh=self.supercell.tolist(),
                    interpolation_max_frobenius=float(np.linalg.norm(full-self.training_blocks, axis=(-2, -1)).max()),
                    truncation_max_frobenius=float(np.linalg.norm(difference, axis=(-2, -1)).max()),
                    retained_hoppings=int(np.count_nonzero(self.origin.d <= self.radius+1e-6)),
                    untruncated_hoppings=int(len(self.origin.t)),
                    certificate="training mesh; independent off-grid comparison reported separately")

    def supercell_csr(self, supercell=None):
        """Periodic cell-gauge CSR, i fastest, orbital index inside each cell.

        This exposes local physical matrices to the 6c supercell/disorder
        route without assuming a diagonal disorder operator or changing S.
        WS images wrapping to the same entry are summed by COO -> CSR.
        """
        supplied = np.asarray(self.supercell if supercell is None else supercell)
        if supplied.shape != (3,) or np.any(supplied < 1) or np.any(supplied != supplied.astype(int)):
            raise ValueError("supercell must contain three positive integers")
        size = supplied.astype(int)
        count, width = int(np.prod(size)), self.origin.J
        selected = self.origin.d <= self.radius+1e-6
        R, i, j, t = (array[selected] for array in
                       (self.origin.R, self.origin.i, self.origin.j, self.origin.t))
        rows, columns, values = [], [], []
        for cell in range(count):
            coordinate = np.array(np.unravel_index(cell, tuple(size), order="F"))
            targets = (coordinate+R) % size
            target_cells = np.ravel_multi_index(targets.T, tuple(size), order="F")
            rows.append(cell*width+i)
            columns.append(target_cells*width+j)
            values.append(t)
        matrix = sparse.coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                                   shape=(count*width, count*width)).tocsr()
        matrix.sum_duplicates()
        matrix.sort_indices()
        return matrix
