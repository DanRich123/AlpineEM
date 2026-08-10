import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from EM2Dsolver import EM2DSolver

# The magnetostatic and electrostatic solver is the same equation.
# The material is the relative permittivity or inverse relative permeability depending on which solver you intend.
# Similarly the source is voltage in volts or Az (vector potential) in Wb/m depending on which solver you use.
# These can be done together or separately.
# For example, if done together, you can use static modes to get TEM profiles from this.
# You can get E and H from the same simulation if they have equal relative permittivity and
# relative permeability, or if neither spatially varies outside the conductors.
# Always plot up the geometry to verify it matches the FDTD layout in ParaView —
# easiest by using metal walls for checking contrast easily.
# boundaries are nodes while the materials are cells - like in FDTD
# E,H fields are returned from the solver in their yee cell positions but then turned into nodes here for plotting down below

# --- Setup ---
nx, ny = 20, 20
dx, dy = 0.5E-3, 0.5E-3
sim = EM2DSolver(nx, ny, dx, dy)

# --- Configuration ---
cx, cy = nx // 2, ny // 2
inner_hw = 1 # Inner square radius (2x2 center block)
outer_hw = 5 # Outer square radius (10x10 outer block)

# --- Set Outer Shield & Exterior (Nodes) ---
for i in range(sim.Nx):
    for j in range(sim.Ny):
        dx_node = abs(i - cx)
        dy_node = abs(j - cy)
        
        # If we are at the sheath boundary or outside of it
        if dx_node >= outer_hw or dy_node >= outer_hw:
            sim.set_boundary(i, j, value=0.0)

# --- Set Inner Conductor (Nodes) ---
for i in range(cx - inner_hw, cx + inner_hw + 1):
    for j in range(cy - inner_hw, cy + inner_hw + 1):
        sim.set_boundary(i, j, value=1.0)

# --- Set Dielectric (Cells) ---
# Filling the entire "gap" area with a specific relative permittivity
# We can just fill the whole area; the boundary nodes will override 
# the fields inside the conductors anyway.
sim.set_material([0, nx], [0, ny], value=3.0) # e.g., Random dielectric value I assigned it

# --- Solve ---
V = sim.solve()
Ex, Ey = sim.get_E_fields()

# Previous settings still in place unless overwritten
# Only thing that needs changing is to overwrite all relative permittivity with inverse relative permeability values

# --- Set Inverse relative permeability (Cells) ---
# Filling the entire "gap" area with a specific INVERSE relative permeability
# We can just fill the whole area; the boundary nodes will override 
# the fields inside the conductors anyway.
sim.set_material([0, nx], [0, ny], value=1.0/(1.0)) # e.g., vacuum

# --- Solve ---
Az = sim.solve()
Hx, Hy = sim.get_H_fields()
I_total = sim.get_net_current()
# print(f"Total simulated current: {I_total:.6f} Amperes")

# --- Interpolation --- to nodes for plotting and so on (E,H are currently in yee cell positions)
def to_nodes_robust(fx, fy, Nx, Ny):
    f_x_node = np.zeros((Nx, Ny))
    if fx.shape == (Nx - 1, Ny):
        f_x_node[1:-1, :] = (fx[:-1, :] + fx[1:, :]) / 2
        f_x_node[0, :]    = fx[0, :]
        f_x_node[-1, :]   = fx[-1, :]
    else:
        f_x_node[:, 1:-1] = (fx[:, :-1] + fx[:, 1:]) / 2
        f_x_node[:, 0]    = fx[:, 0]
        f_x_node[:, -1]   = fx[:, -1]

    f_y_node = np.zeros((Nx, Ny))
    if fy.shape == (Nx, Ny - 1):
        f_y_node[:, 1:-1] = (fy[:, :-1] + fy[:, 1:]) / 2
        f_y_node[:, 0]    = fy[:, 0]
        f_y_node[:, -1]   = fy[:, -1]
    else:
        f_y_node[1:-1, :] = (fy[:-1, :] + fy[1:, :]) / 2
        f_y_node[0, :]    = fy[0, :]
        f_y_node[-1, :]   = fy[-1, :]

    return f_x_node, f_y_node

Ex_n, Ey_n = to_nodes_robust(Ex, Ey, sim.Nx, sim.Ny)
Hx_n, Hy_n = to_nodes_robust(Hx, Hy, sim.Nx, sim.Ny)

# --- Visualization ---
# Average adjacent nodes to get cell-center values (N nodes → N-1 cells)
V_cells  = (V[:-1, :-1] + V[1:, :-1] + V[:-1, 1:] + V[1:, 1:]) / 4
Az_cells = (Az[:-1, :-1] + Az[1:, :-1] + Az[:-1, 1:] + Az[1:, 1:]) / 4

# Cell edges (N+1 points spanning the domain)
x_edges = np.linspace(0, nx * dx, sim.Nx)
y_edges = np.linspace(0, ny * dy, sim.Ny)
X_edges, Y_edges = np.meshgrid(x_edges, y_edges, indexing='ij')

# Electric
plt.figure()
test = plt.pcolormesh(X_edges/dx, Y_edges/dy, V_cells, cmap='nipy_spectral', shading='flat', alpha=0.8)
plt.colorbar(test, label='Voltage (V)')
plt.streamplot(X_edges.T/dx, Y_edges.T/dy, Ex_n.T, Ey_n.T, color='black', linewidth=1, density=1)
plt.title("Electric: $V$ and $\\vec{E}$")
plt.xlabel('Cells in X')
plt.ylabel('Cells in Y')
plt.show()
plt.savefig('electric fields')

# Magnetic
plt.figure()
test = plt.pcolormesh(X_edges/dx, Y_edges/dy, Az_cells, cmap='nipy_spectral', shading='flat', alpha=0.8)
plt.colorbar(test, label='Az Amplitude (Wb/m)')
plt.streamplot(X_edges.T/dx, Y_edges.T/dy, Hx_n.T, Hy_n.T, color='black', linewidth=1, density=1)
plt.title("Magnetic: $Az$ and $\\vec{H}$")
plt.xlabel('Cells in X')
plt.ylabel('Cells in Y')
plt.show()
plt.savefig('magnetic fields')

# --- Output for FDTD ---
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
data=np.zeros((4,nx+1,ny+1), dtype=np.float32)
data[0,:-1,:]=Ex
data[1,:,:-1]=Ey
data[2,:,:-1]=Hx/I_total
data[3,:-1,:]=Hy/I_total
data.flatten(order='F').tofile('gridded_feed.bin')

# --- Characteristic Impedance IF constant mu and ep for TEM mode ---
# make sure you choose the right values here - this section doesn't take material inputs directly from the solver
ep_0= 8.8541878188E-12
mu_0=1.25663706127E-6
eps_r = 3.0
mu_r  = 1.0
#cell centering first
Ex_cc = (Ex[:, :-1] + Ex[:, 1:]) / 2   
Ey_cc = (Ey[:-1, :] + Ey[1:, :]) / 2   
#Hx_cc = (Hx[:-1, :] + Hx[1:, :]) / 2
#Hy_cc = (Hy[:, :-1] + Hy[:, 1:]) / 2
E_sq = Ex_cc**2 + Ey_cc**2
#H_sq = Hx_cc**2 + Hy_cc**2
W_E = 0.5 * ep_0 * eps_r * np.sum(E_sq) * dx * dy
#W_M = 0.5 * mu_0 * mu_r  * np.sum(H_sq) * dx * dy
# C per unit length from electric energy, V0 = 1V boundary
V0 = 1.0
C = 2 * W_E / V0**2
L = mu_0 * mu_r * ep_0 * eps_r / C
Z0 = np.sqrt(L / C)
vp = 1.0 / np.sqrt(L * C)
print(f"C  per unit length: {C:.6e} F/m")
print(f"L  per unit length: {L:.6e} H/m")
print(f"Z0 from L,C:        {Z0:.4f} Ω")
print(f"vp:                 {vp:.6e} m/s")
