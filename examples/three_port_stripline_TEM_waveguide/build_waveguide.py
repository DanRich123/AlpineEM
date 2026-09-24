"""
Voxel-model builder for the 3-port TEM stripline waveguide described in:
Richardson, Dee, Yaeger, Marsh, Westafer, "A New Method for Determining
Antenna Gain via Transmission Line Based Near Field Measurements in a
Waveguide," PIER C, Vol. 155, 61-66, 2025 (Fig. 2).

Produces a 3D NumPy array (dtype=uint8) where 1 = metal, 0 = air.

COORDINATE SYSTEM
------------------
x : along the waveguide's long axis, 0 .. TOTAL_LEN (229 mm)
    x = 0            -> Port 3 end face
    x = TOTAL_LEN     -> Port 1 end face
y : width direction, centered on y = 0 (waveguide's lateral symmetry axis)
z : height direction, centered on z = 0 (septum / center-conductor midplane)

GEOMETRY, FROM THE PAPER (Fig. 2) + USER-SUPPLIED VALUES
-----------------------------------------------------------
Outer conductor (top & bottom ground plates), in-plane (x-y) shape:
    - two 64 mm tapers (width 50 mm -> 200 mm) flanking a
    - 101 mm center straight section (width 200 mm)
    total length = 64 + 101 + 64 = 229 mm

Septum (center conductor), in-plane (x-y) shape:
    - two 87 mm tapers (width 0 mm -> 100 mm, i.e. a point at each port)
      flanking a
    - 55 mm center straight section (width 100 mm)
    total length = 87 + 55 + 87 = 229 mm
    (Note: septum taper breakpoints do NOT line up with the outer
    conductor's taper breakpoints -- that's consistent with the paper,
    each shape is dimensioned independently in Fig. 2.)

Vertical (z) taper -- ***NOT stated in the paper's extracted text***,
supplied by user from inspection of the hardware:
    - septum-to-plate gap (half-height) = 50 mm in the center section
    - tapers down to ~10 mm at each port
    - assumed to follow the SAME x-breakpoints as the outer conductor's
      width taper (64 / 101 / 64 mm), since the plates are the ground
      structure that physically bounds the gap.
    ==> This produces the faceted, hexagonal 3D shape visible in the
        paper's photo (Fig. 2), rather than flat parallel plates.

Metal thickness: 1 voxel (0.5 mm) for top plate, bottom plate, and septum.

Simplifications (per user):
    - All 3 coax connectors/pins (Ports 1, 2, 3) are OMITTED.
    - Port 1 and Port 3 openings (the x=0 and x=TOTAL_LEN end faces) are
      SEALED with a solid metal cap instead of being left open for a
      coax feed.
    - Port 2 (the monopole, normally a pin through the bottom plate) is
      also simply omitted -- the bottom plate is left solid there, no
      hole.
    - Side walls are OPEN (air) -- this is a stripline, not a fully
      enclosed box, matching "highly transmissive" waveguide described
      in the paper.

Everything here is easy to re-tune -- edit the PARAMETERS block below.
"""

import numpy as np

# ============================== PARAMETERS ==============================
VOXEL = 0.5  # mm, cube size (and metal sheet thickness -> 1 voxel)

# --- Outer conductor (top/bottom plates), longitudinal (x) breakpoints ---
TAPER_LEN = 64.0        # mm
CENTER_LEN = 101.0      # mm
TOTAL_LEN = 2 * TAPER_LEN + CENTER_LEN  # 229 mm

WIDTH_CENTER = 200.0    # mm, outer plate width, center section
WIDTH_PORT = 50.0       # mm, outer plate width, at port opening

# --- Vertical taper (septum-to-plate half-gap), user-supplied ---
GAP_CENTER = 50.0       # mm, half-gap at center section (CONFIRMED by user)
GAP_PORT = 10.0         # mm, half-gap at port ends (user ESTIMATE ~10 mm)

# --- Septum (center conductor), longitudinal (x) breakpoints ---
SEPT_TAPER_LEN = 87.0   # mm
SEPT_CENTER_LEN = 55.0  # mm
SEPT_TOTAL_LEN = 2 * SEPT_TAPER_LEN + SEPT_CENTER_LEN  # 229 mm, == TOTAL_LEN

SEPT_WIDTH_CENTER = 100.0  # mm
SEPT_WIDTH_TIP = 0.0       # mm (tapers to a point/line at the port; coax omitted)

# --- Air padding around the structure (keeps open side walls visible) ---
PAD_Y = 20.0   # mm, extra air on each side in y
PAD_Z = 10.0   # mm, extra air above top plate / below bottom plate

# ==========================================================================


def outer_half_width(x):
    """Outer conductor half-width (mm) at longitudinal position x (mm)."""
    if x <= TAPER_LEN:
        w = np.interp(x, [0, TAPER_LEN], [WIDTH_PORT, WIDTH_CENTER])
    elif x <= TAPER_LEN + CENTER_LEN:
        w = WIDTH_CENTER
    else:
        w = np.interp(x, [TAPER_LEN + CENTER_LEN, TOTAL_LEN],
                      [WIDTH_CENTER, WIDTH_PORT])
    return w / 2.0


def outer_half_gap(x):
    """Septum-to-plate half-gap (mm) at longitudinal position x (mm)."""
    if x <= TAPER_LEN:
        g = np.interp(x, [0, TAPER_LEN], [GAP_PORT, GAP_CENTER])
    elif x <= TAPER_LEN + CENTER_LEN:
        g = GAP_CENTER
    else:
        g = np.interp(x, [TAPER_LEN + CENTER_LEN, TOTAL_LEN],
                      [GAP_CENTER, GAP_PORT])
    return g


def septum_half_width(x):
    """Septum half-width (mm) at longitudinal position x (mm)."""
    if x <= SEPT_TAPER_LEN:
        w = np.interp(x, [0, SEPT_TAPER_LEN], [SEPT_WIDTH_TIP, SEPT_WIDTH_CENTER])
    elif x <= SEPT_TAPER_LEN + SEPT_CENTER_LEN:
        w = SEPT_WIDTH_CENTER
    else:
        w = np.interp(x, [SEPT_TAPER_LEN + SEPT_CENTER_LEN, SEPT_TOTAL_LEN],
                      [SEPT_WIDTH_CENTER, SEPT_WIDTH_TIP])
    return w / 2.0


def build_waveguide():
    nx = int(round(TOTAL_LEN / VOXEL))
    max_half_width = WIDTH_CENTER / 2.0 + PAD_Y
    max_half_gap = GAP_CENTER + VOXEL + PAD_Z  # gap + plate thickness + pad

    ny = int(round(2 * max_half_width / VOXEL))
    nz = int(round(2 * max_half_gap / VOXEL))

    vol = np.zeros((nx, ny, nz), dtype=np.uint8)

    y_centers = (np.arange(ny) + 0.5) * VOXEL - max_half_width
    z_centers = (np.arange(nz) + 0.5) * VOXEL - max_half_gap

    plate_vox = max(1, int(round(VOXEL / VOXEL)))  # 1 voxel thick

    for ix in range(nx):
        x = (ix + 0.5) * VOXEL

        hw_outer = outer_half_width(x)
        hg = outer_half_gap(x)
        hw_sept = septum_half_width(x)

        y_mask_outer = np.abs(y_centers) <= hw_outer
        y_mask_sept = np.abs(y_centers) <= hw_sept

        # top & bottom plate: at z = +/- hg (inner face touching the gap),
        # 1 voxel thick, extending outward (away from septum)
        z_top_mask = (z_centers >= hg) & (z_centers < hg + plate_vox * VOXEL)
        z_bot_mask = (z_centers <= -hg) & (z_centers > -hg - plate_vox * VOXEL)

        if y_mask_outer.any():
            vol[ix, y_mask_outer[:, None] & z_top_mask[None, :]] = 1
            vol[ix, y_mask_outer[:, None] & z_bot_mask[None, :]] = 1

        # septum: centered at z = 0, 1 voxel thick
        z_sept_mask = np.abs(z_centers) <= (plate_vox * VOXEL / 2.0)
        if y_mask_sept.any():
            vol[ix, y_mask_sept[:, None] & z_sept_mask[None, :]] = 1

    # --- Seal the two port end faces (coax omitted) ---
    # Fill the entire outer-conductor cross-section (from bottom plate's
    # outer face to top plate's outer face, full outer width) at the
    # very first and very last x-slice.
    for ix in (0, nx - 1):
        x = (ix + 0.5) * VOXEL
        hw_outer = outer_half_width(x)
        hg = outer_half_gap(x)
        y_mask_outer = np.abs(y_centers) <= hw_outer
        z_mask_full = np.abs(z_centers) <= (hg + plate_vox * VOXEL)
        vol[ix, y_mask_outer[:, None] & z_mask_full[None, :]] = 1

    return vol, (VOXEL, VOXEL, VOXEL), (y_centers, z_centers)


if __name__ == "__main__":
    vol, voxel_size, (y_centers, z_centers) = build_waveguide()
    print("Array shape (nx, ny, nz):", vol.shape)
    print("Voxel size (mm):", voxel_size)
    print("Metal voxel count:", int(vol.sum()))
    print("Array size (MB):", vol.nbytes / 1e6)
    np.save("/home/claude/waveguide_geometry.npy", vol)
    print("Saved to waveguide_geometry.npy")
