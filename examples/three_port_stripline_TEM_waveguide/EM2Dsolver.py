import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


class EM2DSolver:
    def __init__(self, nx, ny, dx=1.0, dy=1.0):
        self.nx, self.ny = nx, ny
        self.Nx, self.Ny = nx + 1, ny + 1
        self.dx, self.dy = dx, dy
        self.kappa = np.ones((nx, ny))
        self.fixed = {}

    def set_material(self, i_range, j_range, value):
        self.kappa[i_range[0]:i_range[1], j_range[0]:j_range[1]] = value

    def set_boundary(self, i, j, value):
        """Set a single node (i, j) as a Dirichlet boundary with the given value."""
        self.fixed[(i, j)] = value

    def _get_edge_kappa(self, i, j, axis):
        """
        Return the harmonic mean of the (up to 2) cells sharing the edge.

        axis='x': vertical edge between node-columns i and i+1, at node-row j.
                  The two cells that share this edge are:
                    kappa[i, j]   (cell to the right of the edge, above in j)
                    kappa[i, j-1] (cell to the right of the edge, below in j)
        axis='y': horizontal edge between node-rows j and j+1, at node-column i.
                  The two cells that share this edge are:
                    kappa[i,   j] (cell above the edge, to the right in i)
                    kappa[i-1, j] (cell above the edge, to the left  in i)
        """
        cells = []
        if axis == 'x':
            # Edge sits between node columns i and i+1 (stencil goes left: i-1 → i)
            # Cells sharing this vertical edge: same column-index i, rows j and j-1
            if j < self.ny:
                cells.append(self.kappa[i, j])
            if j > 0:
                cells.append(self.kappa[i, j - 1])
        else:  # axis == 'y'
            # Edge sits between node rows j and j+1 (stencil goes down: j-1 → j)
            # Cells sharing this horizontal edge: same row-index j, columns i and i-1
            if i < self.nx:
                cells.append(self.kappa[i, j])
            if i > 0:
                cells.append(self.kappa[i - 1, j])
        vals = np.array(cells)
        return len(vals) / np.sum(1.0 / vals) if len(vals) > 0 else 1.0

    def solve(self):
        N = self.Nx * self.Ny
        idx = lambda i, j: i * self.Ny + j
        A = lil_matrix((N, N))
        b = np.zeros(N)

        for i in range(self.Nx):
            for j in range(self.Ny):
                k = idx(i, j)
                if (i, j) in self.fixed:
                    A[k, k] = 1.0
                    b[k] = self.fixed[(i, j)]
                    continue

                diag = 0.0

                if i > 0:
                    w = self._get_edge_kappa(i - 1, j, 'x') / (self.dx ** 2)
                    A[k, idx(i - 1, j)] = w
                    diag -= w
                if i < self.Nx - 1:
                    w = self._get_edge_kappa(i, j, 'x') / (self.dx ** 2)
                    A[k, idx(i + 1, j)] = w
                    diag -= w
                if j > 0:
                    w = self._get_edge_kappa(i, j - 1, 'y') / (self.dy ** 2)
                    A[k, idx(i, j - 1)] = w
                    diag -= w
                if j < self.Ny - 1:
                    w = self._get_edge_kappa(i, j, 'y') / (self.dy ** 2)
                    A[k, idx(i, j + 1)] = w
                    diag -= w

                # Guard against isolated nodes (no neighbours, no BC)
                A[k, k] = diag if diag != 0.0 else 1.0

        self.potential = spsolve(A.tocsr(), b).reshape((self.Nx, self.Ny))
        return self.potential

    def get_E_fields(self):
        """
        Calculates Ex and Ey, accounting for the local relative permittivity (kappa) 
        at each Yee cell edge.
        """
        if not hasattr(self, 'potential'):
            raise ValueError("You must call solve() before getting the current.")

        V = self.potential
        Ex = -(V[1:, :] - V[:-1, :]) / self.dx   # shape (nx, Ny)
        Ey = -(V[:, 1:] - V[:, :-1]) / self.dy   # shape (Nx, ny)
        return Ex, Ey

    def get_H_fields(self):
        """
        Calculates Hx and Hy, accounting for the local inverse relative permeability (kappa) 
        at each Yee cell edge.
        """
        if not hasattr(self, 'potential'):
            raise ValueError("You must call solve() before getting the current.")
            
        Az = self.potential
        
        # Hx shape is (Nx, ny) -> lives on horizontal cell edges
        Hx = np.zeros((self.Nx, self.ny))
        for i in range(self.Nx):
            for j in range(self.ny):
                # get local 1/mu via harmonic mean of the cells sharing this horizontal edge
                inv_mu = self._get_edge_kappa(i, j, 'y') 
                Hx[i, j] = -inv_mu * (Az[i, j + 1] - Az[i, j]) / self.dy

        # Hy shape is (nx, Ny) -> lives on vertical cell edges
        Hy = np.zeros((self.nx, self.Ny))
        for i in range(self.nx):
            for j in range(self.Ny):
                # get local 1/mu via harmonic mean of the cells sharing this vertical edge
                inv_mu = self._get_edge_kappa(i, j, 'x')
                Hy[i, j] = inv_mu * (Az[i + 1, j] - Az[i, j]) / self.dx

        return Hx, Hy

    def get_net_current(self):
        """
        Calculates the total net current (in Amperes) wrapping around if using H field simulation
        all fixed boundary nodes by integrating the H-field flux leaving them.
        """
        if not hasattr(self, 'potential'):
            raise ValueError("You must call solve() before getting the current.")
            
        Az = self.potential
        I_total = 0.0

        # Loop through every node in the entire grid
        for (i, j), value in self.fixed.items():
            # Only calculate for nodes where Az is set to a non-zero value 
            # (assuming your inner conductor is non-zero and outer shield is 0.0)
            if value == 0.0:
                continue

            # Check all 4 neighboring nodes. If a neighbor is NOT fixed, 
            # it means we are on the boundary edge facing the air gap.
            
            # Left neighbor
            if i > 0 and (i - 1, j) not in self.fixed:
                inv_mu = self._get_edge_kappa(i, j, 'x')
                H_tangent = inv_mu*(Az[i, j] - Az[i - 1, j]) / self.dx 
                I_total += H_tangent * self.dy
                
            # Right neighbor
            if i < self.Nx - 1 and (i + 1, j) not in self.fixed:
                inv_mu = self._get_edge_kappa(i, j, 'x')
                H_tangent = inv_mu*(Az[i, j] - Az[i + 1, j]) / self.dx 
                I_total += H_tangent * self.dy
                
            # Bottom neighbor
            if j > 0 and (i, j - 1) not in self.fixed:
                inv_mu = self._get_edge_kappa(i, j, 'y') 
                H_tangent = inv_mu*(Az[i, j] - Az[i, j - 1]) / self.dy 
                I_total += H_tangent * self.dx
                
            # Top neighbor
            if j < self.Ny - 1 and (i, j + 1) not in self.fixed:
                inv_mu = self._get_edge_kappa(i, j, 'y') 
                H_tangent = inv_mu*(Az[i, j] - Az[i, j + 1]) / self.dy 
                I_total += H_tangent * self.dx

        return I_total