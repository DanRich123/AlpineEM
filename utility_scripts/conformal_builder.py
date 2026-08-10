import numpy as np
import os

# note that Fortran will be reading the binary created here and it it sensitive to data types. E.g. int32 vs int64 vs float32. It expects a specific formatting.

class ConformalGeometry:
    """
    Builds a fine-grid material map, conformally averages it down to a
    coarse grid, and writes out the files an FDTD solver expects.
    """

    def __init__(self, x_size, y_size, z_size, conform_num, output_dir=None):
        """
        Parameters
        ----------
        x_size, y_size, z_size : int
            Fine-grid dimensions. Each must be evenly divisible by conform_num.
        conform_num : int
            Number of fine cells per coarse cell along each axis (the
            conformal averaging integer, N).
        output_dir : str, optional
            Directory to work/save in. Defaults to the directory this
            script lives in (mirrors the original script's os.chdir call).
        """
        self.x_size = x_size
        self.y_size = y_size
        self.z_size = z_size
        self.N = conform_num

        assert x_size % self.N == 0 and y_size % self.N == 0 and z_size % self.N == 0, \
            f"Grid size ({x_size},{y_size},{z_size}) must be divisible by conform_num ({conform_num})"

        self.coarse_x = x_size // self.N
        self.coarse_y = y_size // self.N
        self.coarse_z = z_size // self.N

        if output_dir is None:
            output_dir = os.path.dirname(os.path.abspath(__file__))
        self.output_dir = output_dir

        # fine grid: material ID per fine cell (0 = ignored / background)
        self.data = np.zeros((x_size, y_size, z_size), dtype=np.float32)

        # coarse grid: material ID per coarse cell, filled in by run()
        self.course_grid = np.zeros((self.coarse_x, self.coarse_y, self.coarse_z), dtype=np.float32)

        # [ [id, epx, epy, epz], ... ] - id 0 is reserved/ignored by the solver
        self.materials_id_info = []

        # populated by run()
        self.coarse_epx = None
        self.coarse_epy = None
        self.coarse_epz = None
        self.coarse_sigx = None
        self.coarse_sigy = None
        self.coarse_sigz = None

    def add_material(self, mat_id, epx, epy, epz, sigx, sigy, sigz):
        """Register a material (don't use id 0 - it's the ignored background)."""
        if mat_id == 0:
            raise ValueError("Material id 0 is reserved for background/ignored cells.")
        self.materials_id_info.append([mat_id, epx, epy, epz, sigx, sigy, sigz])
        return mat_id

    def set_region(self, slices, mat_id):
        """
        Assign a material id to a region of the fine grid.
        """
        self.data[slices] = mat_id

    def _expand_materials_to_blocks(self):
        """Expand per-material eps values onto the full fine grid."""
        blocks_epx = np.ones((self.x_size, self.y_size, self.z_size), dtype=np.float32)
        blocks_epy = np.ones((self.x_size, self.y_size, self.z_size), dtype=np.float32)
        blocks_epz = np.ones((self.x_size, self.y_size, self.z_size), dtype=np.float32)
        blocks_sigx = np.zeros((self.x_size, self.y_size, self.z_size), dtype=np.float32)
        blocks_sigy = np.zeros((self.x_size, self.y_size, self.z_size), dtype=np.float32)
        blocks_sigz = np.zeros((self.x_size, self.y_size, self.z_size), dtype=np.float32)

        for mat_id, epx_val, epy_val, epz_val, sigx_val, sigy_val, sigz_val in self.materials_id_info:
            mask_id = (self.data == mat_id)
            blocks_epx[mask_id] = epx_val
            blocks_epy[mask_id] = epy_val
            blocks_epz[mask_id] = epz_val
            blocks_sigx[mask_id] = sigx_val
            blocks_sigy[mask_id] = sigy_val
            blocks_sigz[mask_id] = sigz_val

        return blocks_epx, blocks_epy, blocks_epz, blocks_sigx, blocks_sigy, blocks_sigz

    def _average_blocks(self, blocks_epx, blocks_epy, blocks_epz, blocks_sigx, blocks_sigy, blocks_sigz):
        """Harmonic mean along each field's own axis, arithmetic along the other two."""
        N = self.N
        Cx, Cy, Cz = self.coarse_x, self.coarse_y, self.coarse_z

        # Reshape (Cx*N, Cy*N, Cz*N) -> (Cx, N, Cy, N, Cz, N)
        # axis pairs: (0,1)=x, (2,3)=y, (4,5)=z ; index 1 in each pair is the within-block axis
        epx6 = blocks_epx.reshape(Cx, N, Cy, N, Cz, N)
        epy6 = blocks_epy.reshape(Cx, N, Cy, N, Cz, N)
        epz6 = blocks_epz.reshape(Cx, N, Cy, N, Cz, N)

        # epx: harmonic along x (axis 1), arithmetic along y,z (axes 3,5 after reduction -> 2,4)
        coarse_epx = 1.0 / np.mean(1.0 / epx6, axis=1)      # -> (Cx, Cy, N, Cz, N)
        coarse_epx = np.mean(coarse_epx, axis=(2, 4))         # -> (Cx, Cy, Cz)

        # epy: harmonic along y (axis 3), arithmetic along x,z (axes 1,5 -> after reduction 1,4)
        coarse_epy = 1.0 / np.mean(1.0 / epy6, axis=3)       # -> (Cx, N, Cy, Cz, N)
        coarse_epy = np.mean(coarse_epy, axis=(1, 4))          # -> (Cx, Cy, Cz)

        # epz: harmonic along z (axis 5), arithmetic along x,y (axes 1,3)
        coarse_epz = 1.0 / np.mean(1.0 / epz6, axis=5)       # -> (Cx, N, Cy, N, Cz)
        coarse_epz = np.mean(coarse_epz, axis=(1, 3))          # -> (Cx, Cy, Cz)

        # same for sigma but also with some catches for divisions by zero
        # Reshape (Cx*N, Cy*N, Cz*N) -> (Cx, N, Cy, N, Cz, N)
        sigx6 = blocks_sigx.reshape(Cx, N, Cy, N, Cz, N)
        sigy6 = blocks_sigy.reshape(Cx, N, Cy, N, Cz, N)
        sigz6 = blocks_sigz.reshape(Cx, N, Cy, N, Cz, N)

        with np.errstate(divide='ignore', invalid='ignore'):
            coarse_sigx = 1.0 / np.mean(1.0 / sigx6, axis=1)
            coarse_sigy = 1.0 / np.mean(1.0 / sigy6, axis=3)
            coarse_sigz = 1.0 / np.mean(1.0 / sigz6, axis=5)

        coarse_sigx = np.mean(coarse_sigx, axis=(2, 4))
        coarse_sigy = np.mean(coarse_sigy, axis=(1, 4))
        coarse_sigz = np.mean(coarse_sigz, axis=(1, 3))

        coarse_sigx = np.nan_to_num(coarse_sigx, nan=0.0, posinf=0.0, neginf=0.0)
        coarse_sigy = np.nan_to_num(coarse_sigy, nan=0.0, posinf=0.0, neginf=0.0)
        coarse_sigz = np.nan_to_num(coarse_sigz, nan=0.0, posinf=0.0, neginf=0.0)

        return (coarse_epx.astype(np.float32, copy=False),
                coarse_epy.astype(np.float32, copy=False),
                coarse_epz.astype(np.float32, copy=False),
                coarse_sigx.astype(np.float32, copy=False),
                coarse_sigy.astype(np.float32, copy=False),
                coarse_sigz.astype(np.float32, copy=False)
                )

    def _map_to_ids(self, coarse_epx, coarse_epy, coarse_epz, coarse_sigx, coarse_sigy, coarse_sigz):
        """Map averaged (epx, epy, epz, sigx, sigy, sigz) back to material IDs, creating new
        'blended' materials as needed."""
        epssig_to_id = {
            (round(float(mat[1]), 5), round(float(mat[2]), 5), round(float(mat[3]), 5), round(float(mat[4]), 5), round(float(mat[5]), 5), round(float(mat[6]), 5)): int(mat[0])
            for mat in self.materials_id_info
        }
        max_id = max([int(mat[0]) for mat in self.materials_id_info] + [0])
        next_id = max_id + 1

        course_grid = np.zeros((self.coarse_x, self.coarse_y, self.coarse_z), dtype=np.float32)

        for i in range(self.coarse_x):
            for j in range(self.coarse_y):
                for k in range(self.coarse_z):
                    key = (
                        round(float(coarse_epx[i, j, k]), 5),
                        round(float(coarse_epy[i, j, k]), 5),
                        round(float(coarse_epz[i, j, k]), 5),
                        round(float(coarse_sigx[i, j, k]), 5),
                        round(float(coarse_sigy[i, j, k]), 5),
                        round(float(coarse_sigz[i, j, k]), 5),
                    )

                    if key == (1.0, 1.0, 1.0, 0.0, 0.0, 0.0):
                        course_grid[i, j, k] = 0
                        continue

                    if key not in epssig_to_id:
                        epssig_to_id[key] = next_id
                        self.materials_id_info.append([next_id, key[0], key[1], key[2], key[3], key[4], key[5]])
                        next_id += 1
                    course_grid[i, j, k] = epssig_to_id[key]

        return course_grid

    """ Fast version that I didn't write myself and haven't validated very well
        def _map_to_ids(self, coarse_epx, coarse_epy, coarse_epz, coarse_sigx, coarse_sigy, coarse_sigz):
        # Map averaged (epx, epy, epz, sigx, sigy, sigz) back to material IDs, creating new
        # 'blended' materials as needed.
        epssig_to_id = {
            (round(float(mat[1]), 5), round(float(mat[2]), 5), round(float(mat[3]), 5),
            round(float(mat[4]), 5), round(float(mat[5]), 5), round(float(mat[6]), 5)): int(mat[0])
            for mat in self.materials_id_info
        }
        max_id = max([int(mat[0]) for mat in self.materials_id_info] + [0])
        next_id = max_id + 1

        shape = coarse_epx.shape  # (coarse_x, coarse_y, coarse_z)

        # Stack and round once, vectorized, in float64 to match float(...) + round(...) precision
        stacked = np.stack(
            [coarse_epx, coarse_epy, coarse_epz, coarse_sigx, coarse_sigy, coarse_sigz],
            axis=-1
        ).astype(np.float64)
        stacked = np.round(stacked, 5)
        flat = stacked.reshape(-1, 6)  # (n_cells, 6)

        # Find unique rows, preserving first-occurrence order (matches original loop's
        # i -> j -> k assignment order for new IDs, since that's C-order/row-major)
        uniq, first_idx, inverse = np.unique(flat, axis=0, return_index=True, return_inverse=True)
        order = np.argsort(first_idx)
        uniq_ordered = uniq[order]

        # Remap inverse indices so they refer to position in uniq_ordered
        rank_of_uniq = np.empty_like(order)
        rank_of_uniq[order] = np.arange(len(order))
        inverse_ordered = rank_of_uniq[inverse]

        # Now only loop over the UNIQUE material combos (typically far fewer than n_cells)
        id_for_rank = np.empty(len(uniq_ordered), dtype=np.float32)
        vacuum_key = (1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        for rank, row in enumerate(uniq_ordered):
            key = tuple(row.tolist())  # 6 python floats, already rounded

            if key == vacuum_key:
                id_for_rank[rank] = 0
                continue

            mat_id = epssig_to_id.get(key)
            if mat_id is None:
                mat_id = next_id
                epssig_to_id[key] = mat_id
                self.materials_id_info.append([mat_id, *key])
                next_id += 1
            id_for_rank[rank] = mat_id

        course_grid = id_for_rank[inverse_ordered].reshape(shape)
        return course_grid
    """

    def run(self):
        """Run the full pipeline: expand materials -> average -> map to IDs."""
        blocks_epx, blocks_epy, blocks_epz, blocks_sigx, blocks_sigy, blocks_sigz= self._expand_materials_to_blocks()
        self.coarse_epx, self.coarse_epy, self.coarse_epz, self.coarse_sigx, self.coarse_sigy, self.coarse_sigz = self._average_blocks(
            blocks_epx, blocks_epy, blocks_epz, blocks_sigx, blocks_sigy, blocks_sigz
        )
        self.course_grid = self._map_to_ids(self.coarse_epx, self.coarse_epy, self.coarse_epz, self.coarse_sigx, self.coarse_sigy, self.coarse_sigz)
        return self.course_grid

    def save(self, bin_name='optional_geom_bulk.bin', mat_name='materials_id.npy'):
        """Write the coarse grid (.bin, Fortran order) and materials table (.npy)."""
        bin_path = os.path.join(self.output_dir, bin_name)
        mat_path = os.path.join(self.output_dir, mat_name)

        self.course_grid.flatten(order='F').tofile(bin_path)
        np.save(mat_path, np.array(self.materials_id_info, dtype=np.float32))
        return bin_path, mat_path

    def plot_slice(self, slice_type , loc):
        """Quick visualization of coarse_grid[:, :, :] (requires matplotlib) at specific location."""
        from matplotlib import pyplot as plt
        X = np.linspace(1, self.coarse_y, self.coarse_y)
        Y = np.linspace(1, self.coarse_z, self.coarse_z)
        X, Y = np.meshgrid(X, Y)
        if slice_type=='x':
            mesh = plt.pcolormesh(X, Y, self.course_grid[loc, :, :])
        if slice_type=='y':
            mesh = plt.pcolormesh(X, Y, self.course_grid[:, loc, :])
        if slice_type=='z':
            mesh = plt.pcolormesh(X, Y, self.course_grid[:, :, loc])
        plt.colorbar(mesh)
        plt.show()
        return mesh