# Geometry Visualization Utilities

Pre-processing visualization tools for generating, inspecting, and styling 3-D FDTD computational domains, material assignments, and boundary setups in ParaView.

---

## Files Overview

| Script / File | Description | Input / Dependencies | Output |
|---|---|---|---|
| `fdtd_geometry_maker.py` | Reads solver configuration (`inputs.txt`) and binary geometry files (`*_geometry.bin`), constructing multiblock VTK rectilinear grids (`.vtr`) wrapped in a master file (`main geometry.vtm`). Converts Yee-cell volumes, thin sheet layers ($x, y, z$), port interfaces, PML/bounding regions, and IGP interfaces into structured VTK blocks. | `inputs.txt`, `*_geometry.bin`, `pyevtk` | `main geometry.vtm` and `unit cell design/*.vtr` files |
| `fdtd_macro.py` | ParaView Python macro that automates rendering of `.vtm` geometry files. Automatically detects active material IDs, applies categorical color maps, configures opacities (e.g., transparent PML/IGP regions), sets up wireframe outlines, and applies standardized annotations. | Active dataset in ParaView, `paraview.simple`, `numpy` | Formatted 3-D rendering, legend, and threshold views inside ParaView |

---

## Technical Notes

* **VTK Grid Mapping:** `fdtd_geometry_maker.py` places thin 2-D sheets ($x, y, z$ directions) and lumped port direction indicators directly on Yee-cell faces using custom sub-cell coordinate arrays (`create_thin_grid()`) to prevent alignment overlaps with volume elements.
* **Material ID Scheme:** Discrete integer ranges map specific simulation domains:
  * `-5`: Internal Ground Plane (IGP)
  * `-4`: Domain Bounding Box
  * `-3`: Perfectly Matched Layers (PML)
  * `-2`: Lumped Port Direction Indicators
  * `-1`: Lumped Ports
  * `≥ 0`: User-Defined Volume & Sheet Material IDs
* **ParaView Pipeline Automation:** `fdtd_macro.py` handles dynamic dataset inspection across multi-block structures, creates threshold filters to hide unassigned background grid cells (ID `-6`), and overrides standard transfer functions with discrete categorical annotations.
