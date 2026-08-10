from conformal_builder import ConformalGeometry

# this script allows the user to build optional geometries for fdtd by drawing material IDs straight into the arrays themselves
# it prepares two files: the .bin file expected by the master.py and fdtd executable, and the materials info (.npy) the user submitted here
# the materials must be registered in the master.py file, the user can use the .npy materials file for easy in registering them in master.py, if desired
# the user first establishes if the grid needs to conform to a smaller size when done
# then the user will specific the grid size, add materials, draw the geometry, and finally create the files needed
# the conformal builder currently doesn't support plasma or other pole types, or permeability yet - they follow the same averaging scheme but the class can't handle more inputs yet

# conformal averaging integer value - reduces cells (e.g. 2x2x2 -> 1x1x1 means conform_num = 2) - select conformal_num = 1 for no averaging
conform_num = 2

# fine grid = FDTD grid * conformal size, in all directions
x_size = 200
y_size = 200
z_size = 200

# initilaize
geom = ConformalGeometry(
    x_size=x_size,
    y_size=y_size,
    z_size=z_size,
    conform_num=conform_num,
    output_dir=None
)

# register materials: (mat_id, epx, epy, epz) - don't use id 0, it's the background filtered in the fdtd executable
geom.add_material(mat_id=2, epx=3, epy=3, epz=3, sigx=0, sigy=0, sigz=0)
geom.add_material(mat_id=5, epx=4, epy=4, epz=4, sigx=0, sigy=0, sigz=0)
geom.add_material(mat_id=3, epx=7, epy=7, epz=7, sigx=0, sigy=0, sigz=0)

# draw custom geometry
geom.set_region((slice(20, 24), slice(None, -4), slice(None)), mat_id=5)
geom.set_region((slice(24, 36), slice(None), slice(None)), mat_id=2)

# other geometries and material files from other programs (E.g. plasma proerties from tecplot interpolation) could be imported and added to conform, if desired

# run the main algorithm that will conform the geometry into the FDTD expected grid size
# if new materials are needed due to conformal averaging, they will be created and registerd in the .npy materials file
geom.run()

# write optional_geom_bulk.bin and materials_id_opfile.npy
bin_path, mat_path = geom.save(bin_name='optional_geom_bulk.bin', mat_name='materials_id_opfile.npy')
print(f"Wrote coarse grid to: {bin_path}")
print(f"Wrote materials table to: {mat_path}")
print(f"Final materials table:\n{geom.materials_id_info}")

# optional: visualize a slice (requires matplotlib installation, but not import here)
geom.plot_slice(slice_type='z', loc=6)