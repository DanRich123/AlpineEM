This example performs FDTD to calculate the monostatic RCS of a near PEC (metal) sphere of radius 12.5mm using AlpineEM. It will output time domain data for post processing and the geometry binary output file for both the object and clear cases.

1) Compile the FDTD software using one of the three options based on intended resources:
	--single threaded CPU
	--multi threaded CPUs (openMP)
	--GPU (open ACC)

2) Setup the master.py file to create a text file of inputs and execute the executable compiled in 1).
	--You must enter the solver name into the master.py script (ACC for GPU version is set to default).
	--Several files will be generated that can be used for post processing and geometry viewing.

3) Setup the master_clear.py file to perform a similar function to 2) but now without the sphere present in the simulation.
	--You must enter the solver name into the master.py script (ACC for GPU version is set to default).
	--Several files will be generated that can be used for post processing and geometry viewing.

4) Run the post_process.py file to generate the RCS data in .csv file. There are hypothetical batch scripts for Slurm that can be modified and used.

5) Lastly, if the user wants to view the geometry, they can run the fdtd_geometry_maker.py file to generate Paraview files. These can be loaded directly into Paraview. There is a macro (fdtd_macro.py) that can be loaded into Paraview for ease it setting up the filters for viewing more easily since the user might not be familiar with Paraview and all of it's options.
	--there will be a single Paraview file and a folder with associated files. Only open the single file, it will reference the folder of files.
	--run the macro from within Paraview by loading it via the Macros tab.

Using my computer and settings I got:
open ACC GPU version took 30 seconds per simulation.
open MP multi-thread CPU version took 4 minutes per simulation.
single threaded default version took 8 minutes per simulation.

This can vary based on many factors so use caution when comparing with other software directly.
