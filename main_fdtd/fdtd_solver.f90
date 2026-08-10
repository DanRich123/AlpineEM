Program fdtd

    !LICENSE FILE IS INCLUDED IN THE PARENT FOLDER OF THIS FILE

    !There are two custom compiler flags - one for spice and one for kmax
    !use_kmax_version
    !use_spice_version
    !This allow for 1 .f90 file for combining my what was previously 4 versions together into 1 program
    !There are likely a number of parameters or arrays that will go unused between versions but declared and/or allocated
    !I will do my best to use directives to reduce memory and improve speed for unused items as I notice them

    !added for spice
#ifdef use_spice_version
    use circuit_mod
    use ngspice_interface_mod
#endif

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup all constants and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !These are non-imported constants, some of which will be updated by import variables
    real, parameter :: c=299792458 !speed of light in m/s
    real, parameter :: mu_0=1.25663706E-6 !permeability of free space
    real, parameter :: ep_0=8.85418782E-12 !permittivity of free space
    real :: spread=0.0 !spread factor of pulse
    real :: t_spread=0.0 !width of pulse
    real :: del_t=0.0 !setup delta t
    real :: min_steps=0.0 !min steps needed
    integer :: i=0,ii=0,j=0,jj=0,k=0,kk=0,counter=0,rr=0 !used for loops
    real :: jj_real=0.0 !used for looping conditions in far field when needing a real variable and not an integer
    real :: del_x=0,del_y=0,del_z=0.0 !step sizes
    integer :: nxPML_1=11,nxPML_2=11,nyPML_1=11,nyPML_2=11,nzPML_1=11,nzPML_2=11 !PML lengths (11 surfaces makes 10 cubes)
    integer, parameter :: m = 3, ma = 1 !Used in pml equations
    real :: sig_x_max=0,sig_y_max=0,sig_z_max=0,alpha_x_max=0
    real :: alpha_y_max=0,alpha_z_max=0,kappa_x_max=0,kappa_y_max=0,kappa_z_max=0.0
    real :: Da=0,Db=0.0 !These are for permability of free space - evetually will be moved with Gax and so on below when made into arrays
    integer :: xlow=0,xhigh=0,ylow=0,yhigh=0,zlow=0,zhigh=0 !used for TF/SF formulation
    integer :: buffer=3 !used for TF/SF formulation and boundary in general
    real :: WHx=0,WHy=0,WHz=0,WEx=0,WEy=0,WEz=0 !weights for plane waves
    real :: WHx_mirror=0,WHy_mirror=0,WHz_mirror=0,WEx_mirror=0,WEy_mirror=0,WEz_mirror=0 !weights for mirrored plane waves
    integer :: y_delay=0,x_delay=0,z_delay=0 !delays for plane waves
    integer :: y_delay_mirror=0,x_delay_mirror=0,z_delay_mirror=0 !delays for mirrored plane waves
    real, parameter :: pi=3.14159265358979323846 !constant pi
    real :: sheet_ep_average_x=0.0 !used for sheet updates
    real :: sheet_sig_average_x=0.0 !used for sheet updates
    real :: sheet_ep_average_y=0.0 !used for sheet updates
    real :: sheet_sig_average_y=0.0 !used for sheet updates
    real :: sheet_ep_average_z=0.0 !used for sheet updates
    real :: sheet_sig_average_z=0.0 !used for sheet updates
    real :: sheet_thickness=0.0 !used for sheet updates
    integer :: port_array_size=0 !determined by inputs to make it long enough for no unwanted reflections
    integer :: pbc_x=0,pbc_y=0,pbc_z=0 !used for shift in unit cell versions if used
    real :: ang1=0, ang2=0 !used for far field items in loops
    integer :: vid_size1=0 !used for video sizing (field outputs)
    integer :: vid_size2=0 !used for video sizing (field outputs)
    integer :: ic=0,jc=0,kc=0 !used for phase centering far fields, among other uses
    integer :: len_far_field_arrays=0 !used for far field array length
    real :: r_for_time_relay=0.0 !far field distance for time delay only
    real :: time_f_var=0 !far field time variable
    character(len=50) :: fdtd_input_file_name !name of file to read into this program
    integer :: counter_x_sheet=0 !these 3 are used for determining if sheets are present for ga,gb combinations
    integer :: counter_y_sheet=0
    integer :: counter_z_sheet=0
    real :: xlow_wall=1.0,ylow_wall=1.0,zlow_wall=1.0,xhigh_wall=1.0,yhigh_wall=1.0 &
    ,zhigh_wall=1.0 !defaults to 1.0 and made zero if relevant for igp condition
    real :: Mx_mirror=1.0,My_mirror=1.0,Mz_mirror=1.0,Jx_mirror=1.0,Jy_mirror=1.0 &
    ,Jz_mirror=1.0 ! defaulted to 1 and altered based on mirror conditioning in setup
    real :: x_mirror_offset=0.0, y_mirror_offset=0.0, z_mirror_offset=0.0
    real :: i_mirror=0.0,j_mirror=0.0,k_mirror=0.0 ! all 6 used for mirroring far field calculations
    real :: use_x_mirror=0.0, use_y_mirror=0.0, use_z_mirror=0.0 !used for all mirroring in far field calculations
    integer :: ff_xlow=0, ff_xhigh=0, ff_ylow=0, ff_yhigh=0, ff_zlow=0, ff_zhigh=0 !far field loop bounds, will setup in setup section
    integer :: data_out_phi=1,data_out_theta=1 !used for far field output array checks
    integer :: count_unique_sheets_x=0 !these 3 used for main fdtd looping for speed improvements
    integer :: count_unique_sheets_y=0
    integer :: count_unique_sheets_z=0
    integer :: plasma_counter=0 !this and several below added for plasmas specifically as names suggest
    real :: temp_plasma_x=0, temp_plasma_y=0, temp_plasma_z=0
    integer :: plasma_min_xfields_zpos,plasma_max_xfields_zpos,plasma_min_xfields_ypos,plasma_max_xfields_ypos &
    ,plasma_min_xfields_xpos,plasma_max_xfields_xpos
    integer :: plasma_min_yfields_zpos,plasma_max_yfields_zpos,plasma_min_yfields_ypos,plasma_max_yfields_ypos &
    ,plasma_min_yfields_xpos,plasma_max_yfields_xpos
    integer :: plasma_min_zfields_zpos,plasma_max_zfields_zpos,plasma_min_zfields_ypos,plasma_max_zfields_ypos &
    ,plasma_min_zfields_xpos,plasma_max_zfields_xpos
    integer :: clock_time_start !these 3 for timing of simulation
    integer :: clock_time_end
    integer :: clock_rate
    real :: step_vc_out !added for antenna source clear cases

    !added for kmax
#ifdef use_kmax_version
    complex :: imag=(0.0,1.0)
    real :: f_adj=0
    complex :: temp_refl_TE=0,temp_refl_TM=0,temp_trans_TE=0,temp_trans_TM=0
    complex :: temp_plasma_x_c=(0.0,0.0), temp_plasma_y_c=(0.0,0.0), temp_plasma_z_c=(0.0,0.0)
    integer :: k_pl_start_E=0
    integer :: k_pl_start_H=0
#endif
#ifndef use_kmax_version
    real :: imag=0.0
    real :: temp_plasma_x_c=0.0, temp_plasma_y_c=0.0, temp_plasma_z_c=0.0
#endif

    !added for spice
#ifdef use_spice_version
    type(circuit_t) :: circuit
    character(len=50) :: netlist
    integer :: error_cnt !this is ununsed unless I am testing
    character(len=60) :: pair_str
    character(len=60) :: command_str
    integer :: spice_var
    real :: finalTime
#endif

    !vars that can be filled in directly from the read in file if needed and some related to them
    integer :: x_size=0
    integer :: y_size=0
    integer :: z_size=0
    real :: step_size_x=0
    real :: step_size_y=0
    real :: step_size_z=0
    integer :: time_steps=0
    real :: f_center=0
    integer :: pulse_type=0
    real :: theta,phi,pol=0
    character(10) :: video_on_char
    integer :: video_on=0
    character(3) :: char_slice
    integer :: slice=0
    integer :: slice_location=0
    real :: imped_port=0
    real :: speed_port=0
    character(30) :: filename !max 30 characters is set for default
    character(30) :: material_type
    integer :: num_total_materials=0
    integer :: num_materials=0
    integer :: materials_id=0
    integer :: num_objects=0
    integer :: num_blocks=0
    integer :: num_spheres=0
    integer :: num_cylinders=0
    character(3) :: char_cylinders
    integer :: num_sheet_materials=0
    integer :: sheet_materials_id=0
    integer :: num_sheets=0
    character(3) :: sheet_position_index
    character(3) :: char_ports
    integer :: num_sheets_x=0
    integer :: num_sheets_y=0
    integer :: num_sheets_z=0
    integer :: num_ports=0
    integer :: num_far_field_angles=0
    real :: plane_wave_amp=0.0
    real :: plane_wave_time_delay=0.0
    real :: antenna_amp=0.0
    real :: antenna_time_delay=0.0
    integer :: excitation_port_number=0 !start at 0 so can't be equal to index unless updated at input section to >=1
    real :: imped_value !imported and used to determine sigma based on fixed sheet thickness
    character(len=15) :: excitation_type !slightly larger than needed but trimmed either way
    character(len=5) :: char_is_mirror !this and next 6 are for mirror plane waves if needed for igp, slightly larger than needed but trimmed either way
    integer :: is_mirror=0
    character(len=2) :: char_mirror_type
    integer :: mirror_type=0
    integer :: mirror_height=0
    real :: theta_mirror=0
    real :: phi_mirror=0
    integer :: is_plasma=0 !added for plasmas specifically
    integer :: num_poles=0 !added for J sources in general
    integer :: num_poles_trial=0 !added for setting max poles submitted by the user
    character(len=7) :: mat_prop_type_char !read in for extra material properties - J sources in general
    character(len=50) :: optional_bulk_geom_filename !this and next 2 added for optional geometry inputs
    character(len=5) :: use_optional_bulk_geom_file
    character(len=7) :: port_type !use for basic vs gridded port discrimination
    integer :: sim_xlow,sim_xhigh,sim_ylow,sim_yhigh,sim_zlow,sim_zhigh !integer for boundary type

    !used in kmax version
#ifdef use_kmax_version
    integer :: mode_type=0
    real :: k_count_y=0.0
    real :: k_count_z=0.0
    real :: k_count_x=0.0
    real :: k_num_y_exception=0.0
    real :: k_num_z_exception=0.0
    real :: k_num_x_exception=0.0
    character(len=5) :: mode_type_char
    character(len=2) :: k_direction_char
    integer :: k_direction=0
#endif

    !added for spice
#ifdef use_spice_version
    character(3) :: char_spice_ports
#endif
    integer :: num_spice_ports=0
    real :: spice_time_reduction_factor


    !these are arrays related to read in items
    real, allocatable :: tmp_real(:,:,:) !used to resize arrays after read-in is finalized after each sub-section
    integer, allocatable :: tmp_int(:,:,:) !used to resize arrays after read-in is finalized after each sub-section
    integer, allocatable :: sheets_x(:,:,:)
    integer, allocatable :: sheets_y(:,:,:)
    integer, allocatable :: sheets_z(:,:,:)
    integer, allocatable :: blocks(:,:,:)
    integer, allocatable :: spheres(:,:,:)
    integer, allocatable :: cylinders(:,:,:)
    real, allocatable :: sheet_properties(:,:,:)
    real, allocatable :: materials_properties(:,:,:)
    real, allocatable :: ports (:,:,:)
    integer, dimension(3) :: p_loc
    integer, dimension(3) :: p_loc_d
    character(15), allocatable :: object_type(:)
    real, allocatable :: far_field_angles(:,:)
    character(30), allocatable :: gridded_feed_names(:)
    real, allocatable :: gridded_feed_all(:,:,:,:)
    integer, allocatable :: temp_mat_id(:)

    !added for spice
#ifdef use_spice_version
    real, allocatable :: ports_spice(:,:,:)
    character(len=15), allocatable :: names_of_spice_ports(:)
    character(len=15), allocatable :: names_of_spice_ports_currents(:)
    character(len=15), allocatable :: C_name(:)
    type(string_t), allocatable :: names(:)
    real, allocatable :: Cap(:)
    character(30), allocatable :: gridded_feed_names_spice(:)
    real, allocatable :: gridded_feed_all_spice(:,:,:,:)
#endif

    !Field and auxillary arrays
#ifndef use_kmax_version
    real, allocatable :: Ex(:,:,:)
    real, allocatable :: Ey(:,:,:)
    real, allocatable :: Ez(:,:,:)
    real, allocatable :: Hx(:,:,:)
    real, allocatable :: Hy(:,:,:)
    real, allocatable :: Hz(:,:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: Ex(:,:,:)
    complex, allocatable :: Ey(:,:,:)
    complex, allocatable :: Ez(:,:,:)
    complex, allocatable :: Hx(:,:,:)
    complex, allocatable :: Hy(:,:,:)
    complex, allocatable :: Hz(:,:,:)
#endif
    real, allocatable :: Gax(:,:,:)
    real, allocatable :: Gbx(:,:,:)
    real, allocatable :: Gay(:,:,:)
    real, allocatable :: Gby(:,:,:)
    real, allocatable :: Gaz(:,:,:)
    real, allocatable :: Gbz(:,:,:)

    !added for any material auxillary differential equations
#ifndef use_kmax_version
    real, allocatable :: J_source_x(:,:,:,:)
    real, allocatable :: J_source_y(:,:,:,:)
    real, allocatable :: J_source_z(:,:,:,:)
    real, allocatable :: Ex_oldt(:,:,:)
    real, allocatable :: Ey_oldt(:,:,:)
    real, allocatable :: Ez_oldt(:,:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: J_source_x(:,:,:,:)
    complex, allocatable :: J_source_y(:,:,:,:)
    complex, allocatable :: J_source_z(:,:,:,:)
    complex, allocatable :: Ex_oldt(:,:,:)
    complex, allocatable :: Ey_oldt(:,:,:)
    complex, allocatable :: Ez_oldt(:,:,:)
#endif
    !added specifically for plasmas
    real, allocatable :: J_plasma_ax(:,:,:,:)
    real, allocatable :: J_plasma_bx(:,:,:,:)
    real, allocatable :: J_plasma_ay(:,:,:,:)
    real, allocatable :: J_plasma_by(:,:,:,:)
    real, allocatable :: J_plasma_az(:,:,:,:)
    real, allocatable :: J_plasma_bz(:,:,:,:)
    real, allocatable :: plasma_freq_cell_x(:,:,:,:)
    real, allocatable :: plasma_freq_x(:,:,:,:)
    real, allocatable :: plasma_loss_cell_x(:,:,:,:)
    real, allocatable :: plasma_loss_x(:,:,:,:)
    real, allocatable :: plasma_freq_cell_y(:,:,:,:)
    real, allocatable :: plasma_freq_y(:,:,:,:)
    real, allocatable :: plasma_loss_cell_y(:,:,:,:)
    real, allocatable :: plasma_loss_y(:,:,:,:)
    real, allocatable :: plasma_freq_cell_z(:,:,:,:)
    real, allocatable :: plasma_freq_z(:,:,:,:)
    real, allocatable :: plasma_loss_cell_z(:,:,:,:)
    real, allocatable :: plasma_loss_z(:,:,:,:)
    
    !voltage and current aux arrays
#ifndef use_kmax_version
    real, allocatable :: Voltage(:,:)
    real, allocatable :: Current(:,:)
    real, allocatable :: Voltage_out(:,:)
    real, allocatable :: Current_out(:,:)
    real, allocatable :: V_inc(:)
    real, allocatable :: C_inc(:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: Voltage(:,:)
    complex, allocatable :: Current(:,:)
    complex, allocatable :: Voltage_out(:,:)
    complex, allocatable :: Current_out(:,:)
    complex, allocatable :: V_inc(:)
    complex, allocatable :: C_inc(:)
#endif

    !added for spice
    !voltages are real because they are recieved from spice - real and imag split before spice calculations in kmax version
#ifdef use_spice_version
    real, allocatable :: Spice_Voltage(:)
    real, allocatable :: Spice_Voltage_out(:,:)
#ifndef use_kmax_version
    real, allocatable :: Spice_Current(:)
    real, allocatable :: Spice_Current_out(:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: Spice_Current(:)
    complex, allocatable :: Spice_Current_out(:,:)
#endif
#endif

    !video arrays
#ifndef use_kmax_version
    real, allocatable :: Ex_video(:,:,:)
    real, allocatable :: Hx_video(:,:,:)
    real, allocatable :: Ey_video(:,:,:)
    real, allocatable :: Hy_video(:,:,:)
    real, allocatable :: Ez_video(:,:,:)
    real, allocatable :: Hz_video(:,:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: Ex_video(:,:,:)
    complex, allocatable :: Hx_video(:,:,:)
    complex, allocatable :: Ey_video(:,:,:)
    complex, allocatable :: Hy_video(:,:,:)
    complex, allocatable :: Ez_video(:,:,:)
    complex, allocatable :: Hz_video(:,:,:)
#endif

    !output field arrays for post processing - unit cell specific - for sparameters
#ifndef use_kmax_version
    real, allocatable :: incident(:)
    real, allocatable :: E_reflected(:)
    real, allocatable :: E_transmitted(:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: incident(:)
    complex, allocatable :: E_reflected_TE(:)
    complex, allocatable :: E_reflected_TM(:)
    complex, allocatable :: E_transmitted_TE(:)
    complex, allocatable :: E_transmitted_TM(:)
#endif

    !For post processing - far field
    !these are sources from fields
#ifndef use_kmax_version
    real, allocatable :: My_xlow(:,:)
    real, allocatable :: Mz_xlow(:,:)
    real, allocatable :: Jy_xlow(:,:)
    real, allocatable :: Jz_xlow(:,:)
    real, allocatable :: My_xlow_oldt(:,:)
    real, allocatable :: Mz_xlow_oldt(:,:)
    real, allocatable :: Jy_xlow_oldt(:,:)
    real, allocatable :: Jz_xlow_oldt(:,:)
    real, allocatable :: My_xhigh(:,:)
    real, allocatable :: Mz_xhigh(:,:)
    real, allocatable :: Jy_xhigh(:,:)
    real, allocatable :: Jz_xhigh(:,:)
    real, allocatable :: My_xhigh_oldt(:,:)
    real, allocatable :: Mz_xhigh_oldt(:,:)
    real, allocatable :: Jy_xhigh_oldt(:,:)
    real, allocatable :: Jz_xhigh_oldt(:,:)
    real, allocatable :: Mx_ylow(:,:)
    real, allocatable :: Mz_ylow(:,:)
    real, allocatable :: Jx_ylow(:,:)
    real, allocatable :: Jz_ylow(:,:)
    real, allocatable :: Mx_ylow_oldt(:,:)
    real, allocatable :: Mz_ylow_oldt(:,:)
    real, allocatable :: Jx_ylow_oldt(:,:)
    real, allocatable :: Jz_ylow_oldt(:,:)
    real, allocatable :: Mx_yhigh(:,:)
    real, allocatable :: Mz_yhigh(:,:)
    real, allocatable :: Jx_yhigh(:,:)
    real, allocatable :: Jz_yhigh(:,:)
    real, allocatable :: Mx_yhigh_oldt(:,:)
    real, allocatable :: Mz_yhigh_oldt(:,:)
    real, allocatable :: Jx_yhigh_oldt(:,:)
    real, allocatable :: Jz_yhigh_oldt(:,:)
    real, allocatable :: Mx_zlow(:,:)
    real, allocatable :: My_zlow(:,:)
    real, allocatable :: Jx_zlow(:,:)
    real, allocatable :: Jy_zlow(:,:)
    real, allocatable :: Mx_zlow_oldt(:,:)
    real, allocatable :: My_zlow_oldt(:,:)
    real, allocatable :: Jx_zlow_oldt(:,:)
    real, allocatable :: Jy_zlow_oldt(:,:)
    real, allocatable :: Mx_zhigh(:,:)
    real, allocatable :: My_zhigh(:,:)
    real, allocatable :: Jx_zhigh(:,:)
    real, allocatable :: Jy_zhigh(:,:)
    real, allocatable :: Mx_zhigh_oldt(:,:)
    real, allocatable :: My_zhigh_oldt(:,:)
    real, allocatable :: Jx_zhigh_oldt(:,:)
    real, allocatable :: Jy_zhigh_oldt(:,:)
    !W and U is how we add fields together in time
    real, allocatable :: Wx(:,:)
    real, allocatable :: Wy(:,:)
    real, allocatable :: Wz(:,:)
    real, allocatable :: Ux(:,:)
    real, allocatable :: Uy(:,:)
    real, allocatable :: Uz(:,:)
    real, allocatable :: W_theta(:,:)
    real, allocatable :: W_phi(:,:)
    real, allocatable :: U_theta(:,:)
    real, allocatable :: U_phi(:,:)
    !Then we combine U and W in theta,phi format to E outputs
    !one of these will be cross pol term
    real, allocatable :: E_theta(:,:)
    real, allocatable :: E_phi(:,:)
    real, allocatable :: E_theta_out(:,:)
    real, allocatable :: E_phi_out(:,:)
    !these are for far field phase centering
    real, allocatable :: Ex_ff_pc(:)
    real, allocatable :: Ey_ff_pc(:)
    real, allocatable :: Ez_ff_pc(:)
    real, allocatable :: Hx_ff_pc(:)
    real, allocatable :: Hy_ff_pc(:)
    real, allocatable :: Hz_ff_pc(:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: My_xlow(:,:)
    complex, allocatable :: Mz_xlow(:,:)
    complex, allocatable :: Jy_xlow(:,:)
    complex, allocatable :: Jz_xlow(:,:)
    complex, allocatable :: My_xlow_oldt(:,:)
    complex, allocatable :: Mz_xlow_oldt(:,:)
    complex, allocatable :: Jy_xlow_oldt(:,:)
    complex, allocatable :: Jz_xlow_oldt(:,:)
    complex, allocatable :: My_xhigh(:,:)
    complex, allocatable :: Mz_xhigh(:,:)
    complex, allocatable :: Jy_xhigh(:,:)
    complex, allocatable :: Jz_xhigh(:,:)
    complex, allocatable :: My_xhigh_oldt(:,:)
    complex, allocatable :: Mz_xhigh_oldt(:,:)
    complex, allocatable :: Jy_xhigh_oldt(:,:)
    complex, allocatable :: Jz_xhigh_oldt(:,:)
    complex, allocatable :: Mx_ylow(:,:)
    complex, allocatable :: Mz_ylow(:,:)
    complex, allocatable :: Jx_ylow(:,:)
    complex, allocatable :: Jz_ylow(:,:)
    complex, allocatable :: Mx_ylow_oldt(:,:)
    complex, allocatable :: Mz_ylow_oldt(:,:)
    complex, allocatable :: Jx_ylow_oldt(:,:)
    complex, allocatable :: Jz_ylow_oldt(:,:)
    complex, allocatable :: Mx_yhigh(:,:)
    complex, allocatable :: Mz_yhigh(:,:)
    complex, allocatable :: Jx_yhigh(:,:)
    complex, allocatable :: Jz_yhigh(:,:)
    complex, allocatable :: Mx_yhigh_oldt(:,:)
    complex, allocatable :: Mz_yhigh_oldt(:,:)
    complex, allocatable :: Jx_yhigh_oldt(:,:)
    complex, allocatable :: Jz_yhigh_oldt(:,:)
    complex, allocatable :: Mx_zlow(:,:)
    complex, allocatable :: My_zlow(:,:)
    complex, allocatable :: Jx_zlow(:,:)
    complex, allocatable :: Jy_zlow(:,:)
    complex, allocatable :: Mx_zlow_oldt(:,:)
    complex, allocatable :: My_zlow_oldt(:,:)
    complex, allocatable :: Jx_zlow_oldt(:,:)
    complex, allocatable :: Jy_zlow_oldt(:,:)
    complex, allocatable :: Mx_zhigh(:,:)
    complex, allocatable :: My_zhigh(:,:)
    complex, allocatable :: Jx_zhigh(:,:)
    complex, allocatable :: Jy_zhigh(:,:)
    complex, allocatable :: Mx_zhigh_oldt(:,:)
    complex, allocatable :: My_zhigh_oldt(:,:)
    complex, allocatable :: Jx_zhigh_oldt(:,:)
    complex, allocatable :: Jy_zhigh_oldt(:,:)
    !W and U is how we add fields together in time
    complex, allocatable :: Wx(:,:)
    complex, allocatable :: Wy(:,:)
    complex, allocatable :: Wz(:,:)
    complex, allocatable :: Ux(:,:)
    complex, allocatable :: Uy(:,:)
    complex, allocatable :: Uz(:,:)
    complex, allocatable :: W_theta(:,:)
    complex, allocatable :: W_phi(:,:)
    complex, allocatable :: U_theta(:,:)
    complex, allocatable :: U_phi(:,:)
    !Then we combine U and W in theta,phi format to E outputs
    !one of these will be cross pol term
    complex, allocatable :: E_theta(:,:)
    complex, allocatable :: E_phi(:,:)
    complex, allocatable :: E_theta_out(:,:)
    complex, allocatable :: E_phi_out(:,:)
    !these are for far field phase centering
    complex, allocatable :: Ex_ff_pc(:)
    complex, allocatable :: Ey_ff_pc(:)
    complex, allocatable :: Ez_ff_pc(:)
    complex, allocatable :: Hx_ff_pc(:)
    complex, allocatable :: Hy_ff_pc(:)
    complex, allocatable :: Hz_ff_pc(:)
#endif
    integer, allocatable :: data_out_time(:)

    !special field array for sheets normal
#ifndef use_kmax_version
    real, allocatable :: Ex_special(:,:,:)
    real, allocatable :: Ey_special(:,:,:)
    real, allocatable :: Ez_special(:,:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: Ex_special(:,:,:)
    complex, allocatable :: Ey_special(:,:,:)
    complex, allocatable :: Ez_special(:,:,:)
#endif
    integer, allocatable :: x_sheet_list(:)
    integer, allocatable :: y_sheet_list(:)
    integer, allocatable :: z_sheet_list(:)

    !These are material arrays that will fill the axuillary arrays (G,D)
    !material and grid point arrays
    !first material cells
    real, allocatable :: relative_ep_x_cell(:,:,:)
    real, allocatable :: sigma_x_cell(:,:,:)
    real, allocatable :: relative_ep_y_cell(:,:,:)
    real, allocatable :: sigma_y_cell(:,:,:)
    real, allocatable :: relative_ep_z_cell(:,:,:)
    real, allocatable :: sigma_z_cell(:,:,:)
    real, allocatable :: sheet_ep_x_cell_x(:,:,:)
    real, allocatable :: sheet_ep_y_cell_x(:,:,:)
    real, allocatable :: sheet_ep_z_cell_x(:,:,:)
    real, allocatable :: sheet_sig_x_cell_x(:,:,:)
    real, allocatable :: sheet_sig_y_cell_x(:,:,:)
    real, allocatable :: sheet_sig_z_cell_x(:,:,:)
    real, allocatable :: sheet_ep_x_cell_y(:,:,:)
    real, allocatable :: sheet_ep_y_cell_y(:,:,:)
    real, allocatable :: sheet_ep_z_cell_y(:,:,:)
    real, allocatable :: sheet_sig_x_cell_y(:,:,:)
    real, allocatable :: sheet_sig_y_cell_y(:,:,:)
    real, allocatable :: sheet_sig_z_cell_y(:,:,:)
    real, allocatable :: sheet_ep_x_cell_z(:,:,:)
    real, allocatable :: sheet_ep_y_cell_z(:,:,:)
    real, allocatable :: sheet_ep_z_cell_z(:,:,:)
    real, allocatable :: sheet_sig_x_cell_z(:,:,:)
    real, allocatable :: sheet_sig_y_cell_z(:,:,:)
    real, allocatable :: sheet_sig_z_cell_z(:,:,:)
    !now grid points
    real, allocatable :: relative_ep_x(:,:,:)
    real, allocatable :: sigma_x(:,:,:)
    real, allocatable :: relative_ep_y(:,:,:)
    real, allocatable :: sigma_y(:,:,:)
    real, allocatable :: relative_ep_z(:,:,:)
    real, allocatable :: sigma_z(:,:,:)
    real, allocatable :: sheet_ep_x_x(:,:,:)
    real, allocatable :: sheet_ep_y_x(:,:,:)
    real, allocatable :: sheet_ep_z_x(:,:,:)
    real, allocatable :: sheet_sig_x_x(:,:,:)
    real, allocatable :: sheet_sig_y_x(:,:,:)
    real, allocatable :: sheet_sig_z_x(:,:,:)
    real, allocatable :: sheet_ep_x_y(:,:,:)
    real, allocatable :: sheet_ep_y_y(:,:,:)
    real, allocatable :: sheet_ep_z_y(:,:,:)
    real, allocatable :: sheet_sig_x_y(:,:,:)
    real, allocatable :: sheet_sig_y_y(:,:,:)
    real, allocatable :: sheet_sig_z_y(:,:,:)
    real, allocatable :: sheet_ep_x_z(:,:,:)
    real, allocatable :: sheet_ep_y_z(:,:,:)
    real, allocatable :: sheet_ep_z_z(:,:,:)
    real, allocatable :: sheet_sig_x_z(:,:,:)
    real, allocatable :: sheet_sig_y_z(:,:,:)
    real, allocatable :: sheet_sig_z_z(:,:,:)

    !denominator arrays - used by regular space and PML
    real, allocatable :: den_ex(:)
    real, allocatable :: den_ey(:)
    real, allocatable :: den_ez(:)
    real, allocatable :: den_hx(:)
    real, allocatable :: den_hy(:)
    real, allocatable :: den_hz(:)

    !PML arrays - not all are used in unit cell version
#ifndef use_kmax_version
    real, allocatable :: psi_Ezx_1(:,:,:)
    real, allocatable :: psi_Hyx_1(:,:,:)
    real, allocatable :: psi_Ezy_1(:,:,:)
    real, allocatable :: psi_Hxy_1(:,:,:)
    real, allocatable :: psi_Exz_1(:,:,:)
    real, allocatable :: psi_Hyz_1(:,:,:)
    real, allocatable :: psi_Eyz_1(:,:,:)
    real, allocatable :: psi_Hxz_1(:,:,:)
    real, allocatable :: psi_Eyx_1(:,:,:)
    real, allocatable :: psi_Hzx_1(:,:,:)
    real, allocatable :: psi_Exy_1(:,:,:)
    real, allocatable :: psi_Hzy_1(:,:,:)
    real, allocatable :: psi_Ezx_2(:,:,:)
    real, allocatable :: psi_Hyx_2(:,:,:)
    real, allocatable :: psi_Ezy_2(:,:,:)
    real, allocatable :: psi_Hxy_2(:,:,:)
    real, allocatable :: psi_Exz_2(:,:,:)
    real, allocatable :: psi_Hyz_2(:,:,:)
    real, allocatable :: psi_Eyz_2(:,:,:)
    real, allocatable :: psi_Hxz_2(:,:,:)
    real, allocatable :: psi_Eyx_2(:,:,:)
    real, allocatable :: psi_Hzx_2(:,:,:)
    real, allocatable :: psi_Exy_2(:,:,:)
    real, allocatable :: psi_Hzy_2(:,:,:)
#endif
#ifdef use_kmax_version
    complex, allocatable :: psi_Ezx_1(:,:,:)
    complex, allocatable :: psi_Hyx_1(:,:,:)
    complex, allocatable :: psi_Ezy_1(:,:,:)
    complex, allocatable :: psi_Hxy_1(:,:,:)
    complex, allocatable :: psi_Exz_1(:,:,:)
    complex, allocatable :: psi_Hyz_1(:,:,:)
    complex, allocatable :: psi_Eyz_1(:,:,:)
    complex, allocatable :: psi_Hxz_1(:,:,:)
    complex, allocatable :: psi_Eyx_1(:,:,:)
    complex, allocatable :: psi_Hzx_1(:,:,:)
    complex, allocatable :: psi_Exy_1(:,:,:)
    complex, allocatable :: psi_Hzy_1(:,:,:)
    complex, allocatable :: psi_Ezx_2(:,:,:)
    complex, allocatable :: psi_Hyx_2(:,:,:)
    complex, allocatable :: psi_Ezy_2(:,:,:)
    complex, allocatable :: psi_Hxy_2(:,:,:)
    complex, allocatable :: psi_Exz_2(:,:,:)
    complex, allocatable :: psi_Hyz_2(:,:,:)
    complex, allocatable :: psi_Eyz_2(:,:,:)
    complex, allocatable :: psi_Hxz_2(:,:,:)
    complex, allocatable :: psi_Eyx_2(:,:,:)
    complex, allocatable :: psi_Hzx_2(:,:,:)
    complex, allocatable :: psi_Exy_2(:,:,:)
    complex, allocatable :: psi_Hzy_2(:,:,:)
#endif
    real, allocatable :: be_x_1(:)
    real, allocatable :: ce_x_1(:)
    real, allocatable :: alphae_x_PML_1(:)
    real, allocatable :: sige_x_PML_1(:)
    real, allocatable :: kappae_x_PML_1(:)
    real, allocatable :: bh_x_1(:)
    real, allocatable :: ch_x_1(:)
    real, allocatable :: alphah_x_PML_1(:)
    real, allocatable :: sigh_x_PML_1(:)
    real, allocatable :: kappah_x_PML_1(:)
    real, allocatable :: be_y_1(:)
    real, allocatable :: ce_y_1(:)
    real, allocatable :: alphae_y_PML_1(:)
    real, allocatable :: sige_y_PML_1(:)
    real, allocatable :: kappae_y_PML_1(:)
    real, allocatable :: bh_y_1(:)
    real, allocatable :: ch_y_1(:)
    real, allocatable :: alphah_y_PML_1(:)
    real, allocatable :: sigh_y_PML_1(:)
    real, allocatable :: kappah_y_PML_1(:)
    real, allocatable :: be_z_1(:)
    real, allocatable :: ce_z_1(:)
    real, allocatable :: alphae_z_PML_1(:)
    real, allocatable :: sige_z_PML_1(:)
    real, allocatable :: kappae_z_PML_1(:)
    real, allocatable :: bh_z_1(:)
    real, allocatable :: ch_z_1(:)
    real, allocatable :: alphah_z_PML_1(:)
    real, allocatable :: sigh_z_PML_1(:)
    real, allocatable :: kappah_z_PML_1(:)
    real, allocatable :: be_x_2(:)
    real, allocatable :: ce_x_2(:)
    real, allocatable :: alphae_x_PML_2(:)
    real, allocatable :: sige_x_PML_2(:)
    real, allocatable :: kappae_x_PML_2(:)
    real, allocatable :: bh_x_2(:)
    real, allocatable :: ch_x_2(:)
    real, allocatable :: alphah_x_PML_2(:)
    real, allocatable :: sigh_x_PML_2(:)
    real, allocatable :: kappah_x_PML_2(:)
    real, allocatable :: be_y_2(:)
    real, allocatable :: ce_y_2(:)
    real, allocatable :: alphae_y_PML_2(:)
    real, allocatable :: sige_y_PML_2(:)
    real, allocatable :: kappae_y_PML_2(:)
    real, allocatable :: bh_y_2(:)
    real, allocatable :: ch_y_2(:)
    real, allocatable :: alphah_y_PML_2(:)
    real, allocatable :: sigh_y_PML_2(:)
    real, allocatable :: kappah_y_PML_2(:)
    real, allocatable :: be_z_2(:)
    real, allocatable :: ce_z_2(:)
    real, allocatable :: alphae_z_PML_2(:)
    real, allocatable :: sige_z_PML_2(:)
    real, allocatable :: kappae_z_PML_2(:)
    real, allocatable :: bh_z_2(:)
    real, allocatable :: ch_z_2(:)
    real, allocatable :: alphah_z_PML_2(:)
    real, allocatable :: sigh_z_PML_2(:)
    real, allocatable :: kappah_z_PML_2(:)

    !for optional read in and the output geometry file
    real, allocatable :: optional_bulk_geom(:,:,:)
    real, allocatable :: output_geometry(:,:,:,:)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Inputs read in section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_start,clock_rate)
    write(*, '("Reading in the inputs file and setting up the simulation...")')

    call get_command_argument(1, fdtd_input_file_name)
    open(1,file=trim(fdtd_input_file_name),status="old",action="read")
        read(1,*) x_size,y_size,z_size
        read(1,*) step_size_x,step_size_y,step_size_z
        read(1,*) time_steps
        read(1,*) spice_time_reduction_factor
        read(1,*) f_center
        read(1,*) pulse_type

        read(1,*) sim_xlow,sim_xhigh,sim_ylow,sim_yhigh,sim_zlow,sim_zhigh
        if (sim_xlow==1 .and. sim_xhigh==1) then
            pbc_x=1
        end if
        if (sim_ylow==1 .and. sim_yhigh==1) then
            pbc_y=1
        end if
        if (sim_zlow==1 .and. sim_zhigh==1) then
            pbc_z=1
        end if

        read(1,*) char_is_mirror
        if (adjustl(trim(char_is_mirror))=='IGP') then
            is_mirror=1
            read(1,*) char_mirror_type
            select case(char_mirror_type)
            case('-x')
                mirror_type=1
            case('+x')
                mirror_type=2
            case('-y')
                mirror_type=3
            case('+y')
                mirror_type=4
            case('-z')
                mirror_type=5
            case('+z')
                mirror_type=6
            end select
            read(1,*) mirror_height
        end if

        !note that in many places throughout it is assumed only one of the excitation types is used at one time
        !this acts as a filter and improved readability in coding by adding the amps and time delays in place
#ifndef use_kmax_version
        read(1,'(a)') excitation_type
        select case (adjustl(trim(excitation_type)))
        case ('plane wave')
            read(1,*) plane_wave_amp, plane_wave_time_delay
            read(1,*) theta,phi,pol
            if (mirror_type==1) then
                theta_mirror=theta
                phi_mirror=180.0-phi
            end if
            if (mirror_type==2) then
                theta_mirror=theta
                phi_mirror=180.0-phi
            end if
            if (mirror_type==3) then
                theta_mirror=theta
                phi_mirror=360.0-phi
            end if
            if (mirror_type==4) then
                theta_mirror=theta
                phi_mirror=360.0-phi
            end if
            if (mirror_type==5) then
                theta_mirror=180.0-theta
                phi_mirror=phi
            end if
            if (mirror_type==6) then
                theta_mirror=180.0-theta
                phi_mirror=phi
            end if
        case ('antenna')
            read(1,*) antenna_amp, antenna_time_delay
            read(1,*) excitation_port_number
        end select

        !print warnings or errors and stop the program in some cases
        if (pbc_x+pbc_y+pbc_z==3) then
            write(*,*) "WARNING: ONLY PBC WALLS DETECTED - ONLY ANTENNA SOURCES CAN BE GENERATED - ONLY OUTPUT FIELDS ARE USEFUL"
        end if
        if (pbc_x+pbc_y+pbc_z==2 .and. is_mirror==1) then
            STOP "ERROR: USER ADDED AN IGP CONDITION TO A UNIT CELL W/ 2 SETS OF PBC - THIS IS NOT NEEDED AND INVALID"
        end if
        if (pbc_x==1 .and. plane_wave_amp>0) then
            if (mod(theta,180.0)/=0.0 .and. mod(phi,180.0)/=90.0) then
                STOP "ERROR: INCIDENT PLANE WAVE W/ OBLIQUE ANGLE W/ RESPECT TO PBC DETECTED - USE KMAX"
            end if
        end if
        if (pbc_y==1 .and. plane_wave_amp>0) then
            if (mod(theta,180.0)/=0.0 .and. mod(phi,180.0)/=0.0) then
                STOP "ERROR: INCIDENT PLANE WAVE W/ OBLIQUE ANGLE W/ RESPECT TO PBC DETECTED - USE KMAX"
            end if
        end if
        if (pbc_z==1 .and. plane_wave_amp>0) then
            if (mod(theta,180.0)/=90.0) then
                STOP "ERROR: INCIDENT PLANE WAVE W/ OBLIQUE ANGLE W/ RESPECT TO PBC DETECTED - USE KMAX"
            end if
        end if
#endif

#ifdef use_kmax_version
        read(1,'(a)') excitation_type
        select case (adjustl(trim(excitation_type)))
        case ('plane wave')
            read(1,*) plane_wave_amp, plane_wave_time_delay
            read(1,*) mode_type_char
            read(1,*) k_direction_char
            if (trim(k_direction_char)=='-') then
                k_direction=0
            else if (trim(k_direction_char)=='+') then
                k_direction=1
            else
                STOP "ERROR: INVALID PROPOGATION DIRECTION DETECTED"
            end if

            !these are used to determine if H or E is excited.
            !also used to determined which field compoenent is used in summation for post processing.
            if (trim(mode_type_char)=='TE') then
                mode_type=0
            end if
            if (trim(mode_type_char)=='TM') then
                mode_type=1
            end if
        case ('antenna')
            read(1,*) antenna_amp, antenna_time_delay
            read(1,*) excitation_port_number
        end select

        if (pbc_y+pbc_z==2) then
            read(1,*) k_count_y, k_count_z
        else if (pbc_x+pbc_z==2) then
            read(1,*) k_count_x, k_count_z
        else if (pbc_x+pbc_y==2) then
            read(1,*) k_count_x, k_count_y
        end if

        !print warnings or errors and stop the program in some cases
        if (pbc_x+pbc_y+pbc_z==3) then
            write(*,*) "WARNING: ONLY PBC WALLS DETECTED - ONLY ANTENNA SOURCES CAN BE GENERATED - ONLY OUTPUT FIELDS ARE USEFUL"
        end if
        if (pbc_x+pbc_y+pbc_z==2 .and. is_mirror==1) then
            STOP "ERROR: USER ADDED AN IGP CONDITION TO A UNIT CELL W/ 2 SETS OF PBC - THIS IS NOT NEEDED AND INVALID"
        end if
        if (pbc_x+pbc_y+pbc_z==1) then
            STOP "ERROR: USER ONLY ADDED 1 PBC W/ KMAX - THIS IS INVALID"
        end if
#endif

        read(1,*) num_far_field_angles
        if (num_far_field_angles>0) then
            allocate(far_field_angles(num_far_field_angles,2))
        else 
            allocate(far_field_angles(1,1))
        end if
        do i=1,num_far_field_angles
            read(1,*) far_field_angles(i,1),far_field_angles(i,2)
            far_field_angles(i,1)=pi/180.0*(far_field_angles(i,1))
            far_field_angles(i,2)=pi/180.0*(far_field_angles(i,2))
        end do

        read(1,*) video_on_char
        if (trim(video_on_char)=='yes') then
            video_on=1
            read(1,*) char_slice
            select case (trim(char_slice))
            case ('x')
                slice=0
            case ('y')
                slice=1
            case ('z')
                slice=2
            end select
            read(1,*) slice_location
        end if
        read(1,*) filename

        read(1,*) num_total_materials
        i=1
        j=1
        if (num_total_materials>0) then
            !max size it could be, will resize below
            allocate(materials_properties(num_total_materials,15,3)) !6 poles allowed currently
            materials_properties(:,:,:) = 0 ! initialized to zero for J source purposes
            allocate(sheet_properties(num_total_materials,4,3))
        end if
        do ii=1, num_total_materials
            read(1,*) material_type
            select case (trim(material_type))
            case('volume')
                num_materials=num_materials+1
                read(1,*) mat_prop_type_char
                read(1,*) materials_properties(i,1,1)
                read(1,*) materials_properties(i,2,1),materials_properties(i,2,2),materials_properties(i,2,3)
                read(1,*) materials_properties(i,3,1),materials_properties(i,3,2),materials_properties(i,3,3)
                if (trim(mat_prop_type_char)=='plasma') then
                    is_plasma=1
                    read(1,*) num_poles_trial
                    do plasma_counter=1, num_poles_trial
                        read(1,*) materials_properties(i,2*plasma_counter+2,1),materials_properties(i,2*plasma_counter+2,2),materials_properties(i,2*plasma_counter+2,3)
                        read(1,*) materials_properties(i,2*plasma_counter+3,1),materials_properties(i,2*plasma_counter+3,2),materials_properties(i,2*plasma_counter+3,3)
                    end do
                    !keep largest number of poles user submits
                    if (num_poles_trial>num_poles) then
                        num_poles=num_poles_trial
                    end if
                end if
                i=i+1
            case ('sheet')
                num_sheet_materials=num_sheet_materials+1
                !old method imported thickness and ep,sig anisotropic for each and was a 'true' way to do it.
                !but this meant having to fix a thickness that would numerically work while then choosing a sigma,ep that matched manufacterer reporting impedance.
                !there was also a constraint for managing the thickness of each layer correctly. I had defaulted to requiring a fixed thickness per elevation.
                !new method here simplifies but restricts the user inputs - should work better in long run I think.
                !I left the bones the same below so that it can be modified easily for future use if we desire a permittivity modification by the user.
                read(1,*) sheet_properties(j,1,1)
                read(1,*) imped_value
                sheet_thickness=min(step_size_x,step_size_y,step_size_z)/5000.0
                sheet_properties(j,3,1)=1.0
                sheet_properties(j,3,2)=1.0
                sheet_properties(j,3,3)=1.0
                if (imped_value==0) then
                    !should account for min 1E-12 thickness or ~1E-9 cell size with copper plating
                    sheet_properties(j,4,1)=1E20
                    sheet_properties(j,4,2)=1E20
                    sheet_properties(j,4,3)=1E20 
                end if
                if (imped_value/=0) then
                    sheet_properties(j,4,1)=1.0/(imped_value*sheet_thickness)
                    sheet_properties(j,4,2)=1.0/(imped_value*sheet_thickness)
                    sheet_properties(j,4,3)=1.0/(imped_value*sheet_thickness)
                end if
                j=j+1
            end select
        end do
        if (num_materials > 0 .and. num_materials < num_total_materials) then
            allocate(tmp_real(num_materials, 15, 3))
            tmp_real = materials_properties(1:num_materials, :, :)
            call move_alloc(tmp_real, materials_properties)
        else if (num_materials == 0 .and. allocated(materials_properties)) then
            deallocate(materials_properties)
        end if
        if (num_sheet_materials < num_total_materials .and. num_sheet_materials>0) then
            allocate(tmp_real(num_sheet_materials, 4, 3))
            tmp_real = sheet_properties(1:num_sheet_materials, :, :)
            call move_alloc(tmp_real, sheet_properties)
        else if (num_sheet_materials == 0 .and. allocated(sheet_properties)) then
            deallocate(sheet_properties)
        end if

        read(1,*) use_optional_bulk_geom_file
        if (trim(use_optional_bulk_geom_file)=='yes') then
            read(1,*) optional_bulk_geom_filename
        end if

        read(1,*) num_objects
        if (num_objects>0) then
            !max size it could be, will resize below
            allocate(object_type(num_objects))
            allocate(spheres(num_objects,3,3))
            allocate(blocks(num_objects,3,3))
            allocate(cylinders(num_objects,4,3))
        end if
        i=1
        j=1
        k=1
        do ii=1, num_objects
            read(1,*) object_type(ii)
            select case (trim(object_type(ii)))
            case('block')
                num_blocks=num_blocks+1
                read(1,*) blocks(i,1,1)
                read(1,*) blocks(i,2,1),blocks(i,2,2),blocks(i,2,3)
                read(1,*) blocks(i,3,1),blocks(i,3,2),blocks(i,3,3)
                i=i+1
            case('sphere')
                num_spheres=num_spheres+1
                read(1,*) spheres(j,1,1)
                read(1,*) spheres(j,2,1),spheres(j,2,2),spheres(j,2,3)
                read(1,*) spheres(j,3,1)
                j=j+1
            case('cylinder')
                num_cylinders=num_cylinders+1
                read(1,*) cylinders(k,1,1)
                read(1,*) char_cylinders
                if (trim(char_cylinders)=='x') then
                    cylinders(k,2,1)=0
                else if (trim(char_cylinders)=='y') then
                    cylinders(k,2,1)=1
                else if (trim(char_cylinders)=='z') then
                    cylinders(k,2,1)=2
                end if
                read(1,*) cylinders(k,3,1),cylinders(k,3,2),cylinders(k,3,3)
                read(1,*) cylinders(k,4,1),cylinders(k,4,2)
                k=k+1
            end select
        end do
        if (num_blocks > 0 .and. num_blocks < num_objects) then
            allocate(tmp_int(num_blocks, 3, 3))
            tmp_int = blocks(1:num_blocks, :, :)
            call move_alloc(tmp_int, blocks)
        else if (num_blocks == 0 .and. allocated(blocks)) then
            deallocate(blocks)
        end if
        if (num_spheres > 0 .and. num_spheres < num_objects) then
            allocate(tmp_int(num_spheres, 3, 3))
            tmp_int = spheres(1:num_spheres, :, :)
            call move_alloc(tmp_int, spheres)
        else if (num_spheres == 0 .and. allocated(spheres)) then
            deallocate(spheres)
        end if
        if (num_cylinders > 0 .and. num_cylinders < num_objects) then
            allocate(tmp_int(num_cylinders, 4, 3))
            tmp_int = cylinders(1:num_cylinders, :, :)
            call move_alloc(tmp_int, cylinders)
        else if (num_cylinders == 0 .and. allocated(cylinders)) then
            deallocate(cylinders)
        end if

        read(1,*) num_sheets
        if (num_sheets>0) then
            allocate(sheets_x(num_sheets,4,2))
            allocate(sheets_y(num_sheets,4,2))
            allocate(sheets_z(num_sheets,4,2))
        end if
        i=1
        j=1
        k=1
        do ii=1, num_sheets
            read(1,*) sheet_position_index
            select case (trim(sheet_position_index))
            case ('x')
                num_sheets_x=num_sheets_x+1
                read(1,*) sheets_x(i,1,1)
                read(1,*) sheets_x(i,2,1)
                read(1,*) sheets_x(i,3,1),sheets_x(i,3,2)
                read(1,*) sheets_x(i,4,1),sheets_x(i,4,2)
                i=i+1
            case ('y')
                num_sheets_y=num_sheets_y+1
                read(1,*) sheets_y(j,1,1)
                read(1,*) sheets_y(j,2,1)
                read(1,*) sheets_y(j,3,1),sheets_y(j,3,2)
                read(1,*) sheets_y(j,4,1),sheets_y(j,4,2)
                j=j+1
            case ('z')
                num_sheets_z=num_sheets_z+1
                read(1,*) sheets_z(k,1,1)
                read(1,*) sheets_z(k,2,1)
                read(1,*) sheets_z(k,3,1),sheets_z(k,3,2)
                read(1,*) sheets_z(k,4,1),sheets_z(k,4,2)
                k=k+1
            end select
        end do
        if (num_sheets_x > 0 .and. num_sheets_x < num_sheets) then
            allocate(tmp_int(num_sheets_x, 4, 2))
            tmp_int = sheets_x(1:num_sheets_x, :, :)
            call move_alloc(tmp_int, sheets_x)
        else if (num_sheets_x == 0 .and. allocated(sheets_x)) then
            deallocate(sheets_x)
        end if
        if (num_sheets_y > 0 .and. num_sheets_y < num_sheets) then
            allocate(tmp_int(num_sheets_y, 4, 2))
            tmp_int = sheets_y(1:num_sheets_y, :, :)
            call move_alloc(tmp_int, sheets_y)
        else if (num_sheets_y == 0 .and. allocated(sheets_y)) then
            deallocate(sheets_y)
        end if
        if (num_sheets_z > 0 .and. num_sheets_z < num_sheets) then
            allocate(tmp_int(num_sheets_z, 4, 2))
            tmp_int = sheets_z(1:num_sheets_z, :, :)
            call move_alloc(tmp_int, sheets_z)
        else if (num_sheets_z == 0 .and. allocated(sheets_z)) then
            deallocate(sheets_z)
        end if

        read(1,*) num_ports
        if (num_ports>0) then
            allocate(ports(num_ports,5,4))
            allocate(gridded_feed_names(num_ports))
        else
            allocate(ports(1,1,1))
        end if
        do i=1, num_ports
            read(1,*) port_type
            if (trim(port_type)=='basic') then
                read(1,*) char_ports
                select case (trim(char_ports))
                case ('x')
                    ports(i,1,1)=0
                case ('y')
                    ports(i,1,1)=1
                case ('z')
                    ports(i,1,1)=2
                end select
                read(1,*) imped_port
                speed_port=c
                !L is 2 and C is 3 - 1 and 4 are unused from old notation in case we ever add back R and G
                ports(i,2,2)=imped_port/speed_port
                ports(i,2,3)=1.0/(imped_port*speed_port)
                read(1,*) ports(i,3,1),ports(i,3,2),ports(i,3,3)
                if (ports(i,3,1)==1 .or. ports(i,3,2)==1 .or. ports(i,3,3)==1) then
                    STOP "ERROR: FDTD SOFTWARE DOES NOT CURRENTLY SUPPORT PLACING AN INTERNAL LUMPED PORT AT A FIRST CUBE"
                end if
                read(1,*) ports(i,4,1),ports(i,4,2),ports(i,4,3) 
                select case (trim(char_ports))
                case ('x')
                    ports(i,4,2)=ports(i,4,2)+1
                    ports(i,4,3)=ports(i,4,3)+1
                case ('y')
                    ports(i,4,1)=ports(i,4,1)+1
                    ports(i,4,3)=ports(i,4,3)+1
                case ('z')
                    ports(i,4,1)=ports(i,4,1)+1
                    ports(i,4,2)=ports(i,4,2)+1
                end select
                ports(i,5,1)=0
            end if
            if (trim(port_type)=='gridded') then
                read(1,*) char_ports
                select case (trim(char_ports))
                case ('-x')
                    ports(i,1,1)=0
                    ports(i,1,2)=0
                case ('+x')
                    ports(i,1,1)=0
                    ports(i,1,2)=1
                case ('-y')
                    ports(i,1,1)=1
                    ports(i,1,2)=2
                case ('+y')
                    ports(i,1,1)=1
                    ports(i,1,2)=3
                case ('-z')
                    ports(i,1,1)=2
                    ports(i,1,2)=4
                case ('+z')
                    ports(i,1,1)=2
                    ports(i,1,2)=5
                end select
                read(1,*) imped_port
                speed_port=c
                !L is 2 and C is 3 - 1 and 4 are unused from old notation in case we ever add back R and G
                ports(i,2,2)=imped_port/speed_port
                ports(i,2,3)=1.0/(imped_port*speed_port)
                read(1,*) ports(i,3,1),ports(i,3,2),ports(i,3,3)
                if (ports(i,3,1)==1 .or. ports(i,3,2)==1 .or. ports(i,3,3)==1) then
                    STOP "ERROR: FDTD SOFTWARE DOES NOT CURRENTLY SUPPORT PLACING AN INTERNAL LUMPED PORT AT A FIRST CUBE"
                end if
                read(1,*) ports(i,4,1),ports(i,4,2)
                read(1,*) gridded_feed_names(i)
                ports(i,5,1)=1
            end if
        end do

        read(1,*) num_spice_ports
#ifdef use_spice_version
        if (num_spice_ports>0) then
            read(1,*) netlist
            allocate(ports_spice(num_spice_ports,5,3))
            allocate(names_of_spice_ports(num_spice_ports))
            allocate(names_of_spice_ports_currents(num_spice_ports))
            allocate(C_name(num_spice_ports))
            allocate(gridded_feed_names_spice(num_spice_ports))
        else
            allocate(ports_spice(1,1,1))
            allocate(names_of_spice_ports_currents(1))
        end if
        do i=1, num_spice_ports
            read(1,*) port_type
            if (trim(port_type)=='basic') then
                read(1,*) char_spice_ports
                select case (trim(char_spice_ports))
                case('x')
                    ports_spice(i,1,1)=0
                case('y')
                    ports_spice(i,1,1)=1
                case('z')
                    ports_spice(i,1,1)=2
                end select
                read(1,*) ports_spice(i,2,1),ports_spice(i,2,2),ports_spice(i,2,3)
                if (ports_spice(i,2,1)==1 .or. ports_spice(i,2,2)==1 .or. ports_spice(i,2,3)==1) then
                    STOP "ERROR: FDTD SOFTWARE DOES NOT CURRENTLY SUPPORT PLACING A SPICE LUMPED PORT AT A FIRST CUBE"
                end if
                read(1,*) ports_spice(i,3,1),ports_spice(i,3,2),ports_spice(i,3,3)
                read(1,*) ports_spice(i,4,1)
                read(1,*) names_of_spice_ports(i)
                read(1,*) names_of_spice_ports_currents(i)
                read(1,*) C_name(i)
                select case (trim(char_spice_ports))
                case ('x')
                    ports_spice(i,3,2)=ports_spice(i,3,2)+1
                    ports_spice(i,3,3)=ports_spice(i,3,3)+1
                case ('y')
                    ports_spice(i,3,1)=ports_spice(i,3,1)+1
                    ports_spice(i,3,3)=ports_spice(i,3,3)+1
                case ('z')
                    ports_spice(i,3,1)=ports_spice(i,3,1)+1
                    ports_spice(i,3,2)=ports_spice(i,3,2)+1
                end select  
                ports_spice(i,5,1)=0
            end if
            if (trim(port_type)=='gridded') then
                read(1,*) char_ports
                select case (trim(char_ports))
                case ('-x')
                    ports_spice(i,1,1)=0
                    ports_spice(i,1,2)=0
                case ('+x')
                    ports_spice(i,1,1)=0
                    ports_spice(i,1,2)=1
                case ('-y')
                    ports_spice(i,1,1)=1
                    ports_spice(i,1,2)=2
                case ('+y')
                    ports_spice(i,1,1)=1
                    ports_spice(i,1,2)=3
                case ('-z')
                    ports_spice(i,1,1)=2
                    ports_spice(i,1,2)=4
                case ('+z')
                    ports_spice(i,1,1)=2
                    ports_spice(i,1,2)=5
                end select
                read(1,*) ports_spice(i,2,1),ports_spice(i,2,2),ports_spice(i,2,3)
                if (ports_spice(i,2,1)==1 .or. ports_spice(i,2,2)==1 .or. ports_spice(i,2,3)==1) then
                    STOP "ERROR: FDTD SOFTWARE DOES NOT CURRENTLY SUPPORT PLACING A SPICE LUMPED PORT AT A FIRST CUBE"
                end if
                read(1,*) ports_spice(i,3,1),ports_spice(i,3,2)
                read(1,*) gridded_feed_names_spice(i)
                read(1,*) ports_spice(i,4,1)
                read(1,*) names_of_spice_ports(i)
                read(1,*) names_of_spice_ports_currents(i)
                ports_spice(i,5,1)=1
            end if
        end do
#endif

#ifdef use_spice_version
    #ifdef use_kmax_version
        if ((mod(num_spice_ports, 2) /= 0)) then
            STOP "ERROR: USER IS MISSING ONE OR MORE SPICE PORTS - THERE SHOULD BE AN EVEN NUMBER"
        end if
    #endif
#endif

        if (allocated(tmp_real)) deallocate(tmp_real)
        if (allocated(tmp_int))  deallocate(tmp_int)

    close(1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Assign some variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Several items here might not necessarily be needed every time but it is a tiny fraction of time and easier to read

    !Establish time increment based on CFL criteria - uses 98% CFL criteria
    !Changes with dimension, PDE order, LTE, etc.
    del_x=step_size_x
    del_y=step_size_y
    del_z=step_size_z
    del_t=0.98/(c*(1.0/del_x**2+1.0/del_y**2+1.0/del_z**2)**0.5)/spice_time_reduction_factor

    !add 1 for loop bc size-1 is the limits
    !add another if pbc because we model an extra layer of cubes
    x_size=x_size+1+pbc_x
    y_size=y_size+1+pbc_y
    z_size=z_size+1+pbc_z

    !far field distance - only used for time delay of far fields so doesn't actually need to be in the far field
    !outputs remove r dependence of magnitude and phase from travel distance
    r_for_time_relay=1.0

    !this is slight overkill of condition - it can technically be a little smaller on lHS of >
    !additionally, the rhs of the = sign could be less too...
    !this updates if 1.0 above is too small
    if ((del_x*x_size+del_y*y_size+del_z*z_size)>r_for_time_relay) then
        r_for_time_relay=(del_x*x_size+del_y*y_size+del_z*z_size)*2
    end if

    !This is slightly larger than it needs to be, I think, but is fine - not a huge time or memory constraint
    !this is because not all array elements are updated - only those that are used. 
    !will revisit this later to improve code maybe
    len_far_field_arrays=int(time_steps+&
    (r_for_time_relay+del_x*x_size+del_y*y_size+del_z*z_size)/(c*del_t))
    
    !size needed for no signal coming back into the simultion + a 100 buffer zone
    !based on speed in a vacuum to get end and back, so should be good indefinitely.
    port_array_size=int(time_steps*del_t*c/min(del_x,del_y,del_z)/2.0+100)

    !if pbc in respective direction then set pml size to zero
    if (pbc_x==1) then
        nxPML_1=0
        nxPML_2=0
    end if
    if (pbc_y==1) then
        nyPML_1=0
        nyPML_2=0
    end if
    if (pbc_z==1) then
        nzPML_1=0
        nzPML_2=0
    end if

    !prep output fields (videos) - shift for pbc if needed
    if (slice==0) then
        vid_size1=y_size-1-pbc_y
        vid_size2=z_size-1-pbc_z
    end if
    if (slice==1) then
        vid_size1=x_size-1-pbc_x
        vid_size2=z_size-1-pbc_z
    end if
    if (slice==2) then
        vid_size1=x_size-1-pbc_x
        vid_size2=y_size-1-pbc_y
    end if

    !convert to radians 
    theta=pi/180.0*(theta)
    phi=pi/180.0*(phi)
    pol=pi/180.0*(pol)

    !mirror conversion 
    theta_mirror=pi/180.0*(theta_mirror)
    phi_mirror=pi/180.0*(phi_mirror)

    !far field data uses coordinate system centered in the middle of the simulation space we care about (and the user knows about)
    ic=int((x_size-1-pbc_x)/2.0)
    jc=int((y_size-1-pbc_y)/2.0)
    kc=int((z_size-1-pbc_z)/2.0)

    !low and high vars are used in a number of places - they define the TF/SF boundary layers
    !the low/high are defaulted to zero so set them if they are needed
    xlow=1
    xhigh=x_size-1
    ylow=1
    yhigh=y_size-1
    zlow=1
    zhigh=z_size-1
    if (pbc_x==0) then
        xlow=nxPML_1+buffer
        xhigh=x_size-nxPML_2-buffer
    end if
    if (pbc_y==0) then
        ylow=nyPML_1+buffer
        yhigh=y_size-nyPML_2-buffer 
    end if
    if (pbc_z==0) then  
        zlow=nzPML_1+buffer
        zhigh=z_size-nzPML_2-buffer
    end if

    !if mirror is used then there are several parameters I use that are now ready to be set before moving on.
    !based on direction of which plane +-x,y,z
    !recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1) then
        select case(mirror_type)
        case(1)
            xlow=mirror_height
            xlow_wall=0
            x_mirror_offset = 2.0 * (xlow)
            ic=xlow
            nxPML_1=0
            use_x_mirror=1.0
        case(2)
            xhigh=mirror_height
            xhigh_wall=0
            x_mirror_offset = 2.0 * (xhigh)
            ic=xhigh
            nxPML_2=0
            use_x_mirror=1.0
        case(3)
            ylow=mirror_height
            ylow_wall=0
            y_mirror_offset = 2.0 * (ylow)
            jc=ylow
            nyPML_1=0
            use_y_mirror=1.0
        case(4)
            yhigh=mirror_height
            yhigh_wall=0
            y_mirror_offset = 2.0 * (yhigh)
            jc=yhigh
            nyPML_2=0
            use_y_mirror=1.0
        case(5)
            zlow=mirror_height
            zlow_wall=0
            z_mirror_offset = 2.0 * (zlow)
            kc=zlow
            nzPML_1=0
            use_z_mirror=1.0
        case(6)
            zhigh=mirror_height
            zhigh_wall=0
            z_mirror_offset = 2.0 * (zhigh)
            kc=zhigh
            nzPML_2=0
            use_z_mirror=1.0
        end select
    end if

    !delays for plane waves
    !varying with incident angles, so this is a lookup table
    if ((theta>=0.0) .and. (theta<=(pi/2.0))) then
        if ((phi>=0.0) .and. (phi<=(pi/2.0))) then
            x_delay=xlow
            y_delay=ylow
            z_delay=zlow
        end if
        if ((phi>(pi/2.0)) .and. (phi<=pi)) then
            x_delay=xhigh
            y_delay=ylow
            z_delay=zlow
        end if
        if ((phi>pi) .and. (phi<=(3.0*pi/2.0))) then
            x_delay=xhigh
            y_delay=yhigh
            z_delay=zlow
        end if
        if ((phi>(3.0*pi/2.0)) .and. (phi<(2.0*pi))) then
            x_delay=xlow
            y_delay=yhigh
            z_delay=zlow
        end if
    end if
    if ((theta>(pi/2.0)) .and. (theta<=pi)) then
        if ((phi>=0) .and. (phi<=(pi/2.0))) then
            x_delay=xlow
            y_delay=ylow
            z_delay=zhigh
        end if
        if ((phi>(pi/2.0)) .and. (phi<=pi)) then
            x_delay=xhigh
            y_delay=ylow
            z_delay=zhigh
        end if
        if ((phi>pi) .and. (phi<=(3.0*pi/2.0))) then
            x_delay=xhigh
            y_delay=yhigh
            z_delay=zhigh
        end if
        if ((phi>(3.0*pi/2.0)) .and. (phi<(2.0*pi))) then
            x_delay=xlow
            y_delay=yhigh
            z_delay=zhigh
        end if
    end if

    !weights for plane waves
    !default is 1 V/m but user can modify
#ifndef use_kmax_version
    WHx=plane_wave_amp*(sin(pol)*sin(phi)+cos(pol)*cos(theta)*cos(phi))*sqrt(ep_0/mu_0)
    WHy=plane_wave_amp*(-1*sin(pol)*cos(phi)+cos(pol)*cos(theta)*sin(phi))*sqrt(ep_0/mu_0)
    WHz=plane_wave_amp*(-1*cos(pol)*sin(theta))*sqrt(ep_0/mu_0)
    WEx=plane_wave_amp*(cos(pol)*sin(phi)-sin(pol)*cos(theta)*cos(phi))
    WEy=plane_wave_amp*(-1*cos(pol)*cos(phi)-sin(pol)*cos(theta)*sin(phi))
    WEz=plane_wave_amp*(sin(pol)*sin(theta))
#endif
#ifdef use_kmax_version
    !factors of 2 preferred because it launches both directions bc no TF/SF used - doesn't exactly 'split' depending on angle
    if (mode_type==0) then !Then it's a TE excitation and we need E amplitudes - they default to zero otherwise
        WEx=2.0*plane_wave_amp
        WEy=2.0*plane_wave_amp
        WEz=2.0*plane_wave_amp
    end if    
    if (mode_type==1) then !Then it's a TM excitation and we need H amplitudes - they default to zero otherwise
        WHx=2.0*sqrt(ep_0/mu_0)*plane_wave_amp
        WHz=2.0*sqrt(ep_0/mu_0)*plane_wave_amp
        WHy=2.0*sqrt(ep_0/mu_0)*plane_wave_amp
    end if
#endif

    !delay_mirrors for mirrored plane waves
    !varying with incident angles, so this is a lookup table
    if ((theta_mirror>=0.0) .and. (theta_mirror<=(pi/2.0))) then
        if ((phi_mirror>=0.0) .and. (phi_mirror<=(pi/2.0))) then
            x_delay_mirror=xlow
            y_delay_mirror=ylow
            z_delay_mirror=zlow
        end if
        if ((phi_mirror>(pi/2.0)) .and. (phi_mirror<=pi)) then
            x_delay_mirror=xhigh
            y_delay_mirror=ylow
            z_delay_mirror=zlow
        end if
        if ((phi_mirror>pi) .and. (phi_mirror<=(3.0*pi/2.0))) then
            x_delay_mirror=xhigh
            y_delay_mirror=yhigh
            z_delay_mirror=zlow
        end if
        if ((phi_mirror>(3.0*pi/2.0)) .and. (phi_mirror<(2.0*pi))) then
            x_delay_mirror=xlow
            y_delay_mirror=yhigh
            z_delay_mirror=zlow
        end if
    end if
    if ((theta_mirror>(pi/2.0)) .and. (theta_mirror<=pi)) then
        if ((phi_mirror>=0) .and. (phi_mirror<=(pi/2.0))) then
            x_delay_mirror=xlow
            y_delay_mirror=ylow
            z_delay_mirror=zhigh
        end if
        if ((phi_mirror>(pi/2.0)) .and. (phi_mirror<=pi)) then
            x_delay_mirror=xhigh
            y_delay_mirror=ylow
            z_delay_mirror=zhigh
        end if
        if ((phi_mirror>pi) .and. (phi_mirror<=(3.0*pi/2.0))) then
            x_delay_mirror=xhigh
            y_delay_mirror=yhigh
            z_delay_mirror=zhigh
        end if
        if ((phi_mirror>(3.0*pi/2.0)) .and. (phi_mirror<(2.0*pi))) then
            x_delay_mirror=xlow
            y_delay_mirror=yhigh
            z_delay_mirror=zhigh
        end if
    end if

    !specific to mirror delays only.
    !we are now ready to define more mirror variables if needed.
    !need to shift the delay so origin makes PEC plane halfway.
    !recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1) then
        select case(mirror_type)
        case(1)
            x_delay_mirror=mirror_height+(mirror_height-xhigh)
        case(2)
            x_delay_mirror=mirror_height+(mirror_height-xlow)
        case(3)
            y_delay_mirror=mirror_height+(mirror_height-yhigh)
        case(4)
            y_delay_mirror=mirror_height+(mirror_height-ylow)
        case(5)
            z_delay_mirror=mirror_height+(mirror_height-zhigh)
        case(6)
            z_delay_mirror=mirror_height+(mirror_height-zlow)
        end select
    end if

    !weights for mirrored plane waves if needed
    !recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (is_mirror==1) then
        if ((mirror_type==1) .or. (mirror_type==2)) then !yz-plane is ground, x is normal

            WEx_mirror = WEx         ! Normal to ground - not flipped for PEC
            WEy_mirror = -1.0 * WEy  ! Tangential - flip for PEC
            WEz_mirror = -1.0 * WEz  ! Tangential - flip for PEC
            WHx_mirror = -1.0 * WHx  ! Normal - flip naturally to cancel
            WHy_mirror = WHy         ! Tangential - don't flip
            WHz_mirror = WHz         ! Tangential - don't flip

            !Opposite from E,H because of cross product
            Jy_mirror=-1.0
            Jz_mirror=-1.0
            Mx_mirror=-1.0
        end if

        if ((mirror_type==3) .or. (mirror_type==4)) then !xz-plane is ground, y is normal

            WEx_mirror = -1.0 * WEx  ! Tangential - flip
            WEy_mirror = WEy         ! Normal - not flipped
            WEz_mirror = -1.0 * WEz  ! Tangential - flip
            WHx_mirror = WHx         ! Tangential - don't flip
            WHy_mirror = -1.0 * WHy  ! Normal - flip naturally
            WHz_mirror = WHz         ! Tangential - don't flip

            !Opposite from E,H because of cross product
            Jx_mirror=-1.0
            Jz_mirror=-1.0
            My_mirror=-1.0
        end if

        if ((mirror_type==5) .or. (mirror_type==6)) then !xy-plane is ground, z is normal

            WEx_mirror = -1.0 * WEx  ! Tangential - flip
            WEy_mirror = -1.0 * WEy  ! Tangential - flip
            WEz_mirror = WEz         ! Normal - not flipped
            WHx_mirror = WHx         ! Tangential - don't flip
            WHy_mirror = WHy         ! Tangential - don't flip
            WHz_mirror = -1.0 * WHz  ! Normal - flip naturally

            !Opposite from E,H because of cross product
            Jx_mirror=-1.0
            Jy_mirror=-1.0
            Mz_mirror=-1.0
        end if
    end if

    !Set far field looping bounds
    ff_xlow=xlow-2
    ff_xhigh=xhigh+1
    ff_ylow=ylow-2
    ff_yhigh=yhigh+1
    ff_zlow=zlow-2
    ff_zhigh=zhigh+1

    !if mirror then we are changing how we handle one of the walls
    !this ensures the 4 perpendicular walls will only loop until the pec plane, not on or through it.
    !recall mirror types are 1,2,3,4,5,6 for -x,+x,-y,+y,-z,+z
    if (mirror_type==1) then
        ff_xlow=xlow+1
    end if
    if (mirror_type==2) then
        ff_xhigh=xhigh-2
    end if
    if (mirror_type==3) then
        ff_ylow=ylow+1
    end if
    if (mirror_type==4) then
        ff_yhigh=yhigh-2
    end if
    if (mirror_type==5) then
        ff_zlow=zlow+1
    end if
    if (mirror_type==6) then
        ff_zhigh=zhigh-2
    end if

    !if pbc condition applies in a direction, then overwrite the bounds
    !these are used for S parameters sweep as well if unit cell is present
    if (pbc_x==1) then
        ff_xlow=1
        ff_xhigh=x_size-2
    end if
    if (pbc_y==1) then
        ff_ylow=1
        ff_yhigh=y_size-2
    end if
    if (pbc_z==1) then
        ff_zlow=1
        ff_zhigh=z_size-2
    end if
    
    !These are used in pml equations
    !assumed +-x (or y or z) have the same size
    !The log value can increase with better precision
#ifndef use_kmax_version
    sig_x_max = -log(1.0E-6)*((m+1)/(del_x*(mu_0/ep_0*1.0)**0.5))/(2*nxPML_1)
    sig_y_max = -log(1.0E-6)*((m+1)/(del_y*(mu_0/ep_0*1.0)**0.5))/(2*nyPML_1)
    sig_z_max = -log(1.0E-6)*((m+1)/(del_z*(mu_0/ep_0*1.0)**0.5))/(2*nzPML_1)
    alpha_x_max = 0.03
    alpha_y_max = alpha_x_max
    alpha_z_max = alpha_x_max
    kappa_x_max = 7.0
    kappa_y_max = kappa_x_max
    kappa_z_max = kappa_x_max
#endif
#ifdef use_kmax_version
    sig_x_max = -log(1.0E-6)*((m+1)/(del_x*(mu_0/ep_0*1.0)**0.5))/(2*nxPML_1)*1.25
    sig_y_max = -log(1.0E-6)*((m+1)/(del_y*(mu_0/ep_0*1.0)**0.5))/(2*nyPML_1)*1.25
    sig_z_max = -log(1.0E-6)*((m+1)/(del_z*(mu_0/ep_0*1.0)**0.5))/(2*nzPML_1)*1.25
    alpha_x_max = 0.03+0.9*(sqrt(k_count_y**2+k_count_z**2+k_count_x**2))/377.0
    alpha_y_max = alpha_x_max
    alpha_z_max = alpha_x_max
    kappa_x_max = 7.0
    kappa_y_max = kappa_x_max
    kappa_z_max = kappa_x_max
#endif

    !Setup pulse information based on type - these have historically been good
    !can add pw and antenna amp time delay cause one will be zero always
#ifndef use_kmax_version
    if (pulse_type==1) then
        spread=1.0/(2.0*pi*f_center)
        t_spread=5.0*spread+plane_wave_time_delay+antenna_time_delay
    end if
    if (pulse_type==2) then
        spread=1.0/(2.0*pi*f_center)
        t_spread=6.0*spread+plane_wave_time_delay+antenna_time_delay
    end if
#endif
#ifdef use_kmax_version
    !currently independent of pulse type - but we might change this
    spread=4.0/(2.0*pi*f_center)
    t_spread=5.0*spread+plane_wave_time_delay+antenna_time_delay
    f_adj=(sqrt(k_count_y**2+k_count_z**2+k_count_x**2)*c)/(2.0*pi)+f_center
#endif

#ifdef use_kmax_version
        !need some exceptions for plane generation and wave port measurment for zeros and for starting plane setup
        if (pbc_y+pbc_z==2) then
            if (k_direction==1) then
                k_pl_start_H=xlow-1
                k_pl_start_E=xlow
            else if (k_direction==0) then
                k_pl_start_H=xhigh
                k_pl_start_E=xhigh
            end if
            if (k_count_y==0 .and. k_count_z/=0) then
                k_num_y_exception=0.0
                k_num_z_exception=(k_count_z/sqrt(k_count_z**2+k_count_y**2))
            end if
            if (k_count_y/=0 .and. k_count_z==0) then
                k_num_y_exception=(k_count_y/sqrt(k_count_z**2+k_count_y**2))
                k_num_z_exception=0.0
            end if
            if ((k_count_y==0) .and. (k_count_z==0)) then
                !This will create y directed fields for TM or TE mode that is used
                k_num_y_exception=0.0
                k_num_z_exception=1.0
            end if
            if ((k_count_y/=0) .and. (k_count_z/=0)) then
                k_num_y_exception=(k_count_y/sqrt(k_count_z**2+k_count_y**2))
                k_num_z_exception=(k_count_z/sqrt(k_count_z**2+k_count_y**2))
            end if
        end if
        if (pbc_x+pbc_z==2) then
            if (k_direction==1) then
                k_pl_start_H=ylow-1
                k_pl_start_E=ylow
            else if (k_direction==0) then
                k_pl_start_H=yhigh
                k_pl_start_E=yhigh
            end if
            if (k_count_x==0 .and. k_count_z/=0) then
                k_num_x_exception=0.0
                k_num_z_exception=(k_count_z/sqrt(k_count_z**2+k_count_x**2))
            end if
            if (k_count_x/=0 .and. k_count_z==0) then
                k_num_x_exception=(k_count_x/sqrt(k_count_z**2+k_count_x**2))
                k_num_z_exception=0.0
            end if
            if ((k_count_x==0) .and. (k_count_z==0)) then
                !This will create x directed fields for TM or TE mode that is used
                k_num_x_exception=0.0
                k_num_z_exception=1.0
            end if
            if ((k_count_x/=0) .and. (k_count_z/=0)) then
                k_num_x_exception=(k_count_x/sqrt(k_count_z**2+k_count_x**2))
                k_num_z_exception=(k_count_z/sqrt(k_count_z**2+k_count_x**2))
            end if
        end if
        if (pbc_x+pbc_y==2) then
            if (k_direction==1) then
                k_pl_start_H=zlow-1
                k_pl_start_E=zlow
            else if (k_direction==0) then
                k_pl_start_H=zhigh
                k_pl_start_E=zhigh
            end if
            if (k_count_x==0 .and. k_count_y/=0) then
                k_num_x_exception=0.0
                k_num_y_exception=(k_count_y/sqrt(k_count_y**2+k_count_x**2))
            end if
            if (k_count_x/=0 .and. k_count_y==0) then
                k_num_x_exception=(k_count_x/sqrt(k_count_y**2+k_count_x**2))
                k_num_y_exception=0.0
            end if
            if ((k_count_x==0) .and. (k_count_y==0)) then
                !This will create x directed fields for TM or TE mode that is used
                k_num_x_exception=0.0
                k_num_y_exception=1.0
            end if
            if ((k_count_x/=0) .and. (k_count_y/=0)) then
                k_num_x_exception=(k_count_x/sqrt(k_count_y**2+k_count_x**2))
                k_num_y_exception=(k_count_y/sqrt(k_count_y**2+k_count_x**2))
            end if
        end if
#endif

    !Establish minimum steps required - but likely needs way more - prints to user in .dat to help orient them
    !kmax won't print this because many angles excited so it's less useful in that case
    min_steps=2*t_spread/del_t+sqrt(((sin(theta)*cos(phi)*(x_size)*del_x)**2 &
    +(sin(theta)*sin(phi)*(y_size)*del_y)**2+(cos(theta)*(z_size)*del_z)**2))/(c*del_t)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setting up the geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !permittivity and permability arrays will be created from geometry inputs
    !ultimiately, I will be create the G,D, and J arrays (if needed) after this section

    !I allocate dummy arrays when unused - this lets me avoid initilization and acc modifications throughout
    !this is true for many array allocations throughout, even outside this section

    !need to allocate and initialize several material arrays
    
    !material and grid point arrays
    !material cells
    allocate(relative_ep_x_cell(x_size-1,y_size-1,z_size-1)) 
    allocate(relative_ep_y_cell(x_size-1,y_size-1,z_size-1)) 
    allocate(relative_ep_z_cell(x_size-1,y_size-1,z_size-1)) 
    allocate(sigma_x_cell(x_size-1,y_size-1,z_size-1)) 
    allocate(sigma_y_cell(x_size-1,y_size-1,z_size-1)) 
    allocate(sigma_z_cell(x_size-1,y_size-1,z_size-1)) 
    !grid points
    allocate(relative_ep_x(x_size-1,y_size-1,z_size-1)) 
    allocate(relative_ep_y(x_size-1,y_size-1,z_size-1)) 
    allocate(relative_ep_z(x_size-1,y_size-1,z_size-1))
    allocate(sigma_x(x_size-1,y_size-1,z_size-1)) 
    allocate(sigma_y(x_size-1,y_size-1,z_size-1)) 
    allocate(sigma_z(x_size-1,y_size-1,z_size-1)) 
    if (num_sheets_x>0) then
        !material cells
        allocate(sheet_ep_x_cell_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_cell_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_cell_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_cell_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_cell_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_cell_x(x_size-1,y_size-1,z_size-1))
        !grid points
        allocate(sheet_ep_x_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_x(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_x(x_size-1,y_size-1,z_size-1))
    else
        !material cells
        allocate(sheet_ep_x_cell_x(1,1,1))
        allocate(sheet_ep_y_cell_x(1,1,1))
        allocate(sheet_ep_z_cell_x(1,1,1))
        allocate(sheet_sig_x_cell_x(1,1,1))
        allocate(sheet_sig_y_cell_x(1,1,1))
        allocate(sheet_sig_z_cell_x(1,1,1))
        !grid points
        allocate(sheet_ep_x_x(1,1,1))
        allocate(sheet_ep_y_x(1,1,1))
        allocate(sheet_ep_z_x(1,1,1))
        allocate(sheet_sig_x_x(1,1,1))
        allocate(sheet_sig_y_x(1,1,1))
        allocate(sheet_sig_z_x(1,1,1))
    end if
    if (num_sheets_y>0) then
        !material cells
        allocate(sheet_ep_x_cell_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_cell_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_cell_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_cell_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_cell_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_cell_y(x_size-1,y_size-1,z_size-1))
        !grid points
        allocate(sheet_ep_x_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_y(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_y(x_size-1,y_size-1,z_size-1))
    else 
        !material cells
        allocate(sheet_ep_x_cell_y(1,1,1))
        allocate(sheet_ep_y_cell_y(1,1,1))
        allocate(sheet_ep_z_cell_y(1,1,1))
        allocate(sheet_sig_x_cell_y(1,1,1))
        allocate(sheet_sig_y_cell_y(1,1,1))
        allocate(sheet_sig_z_cell_y(1,1,1))
        !grid points
        allocate(sheet_ep_x_y(1,1,1))
        allocate(sheet_ep_y_y(1,1,1))
        allocate(sheet_ep_z_y(1,1,1))
        allocate(sheet_sig_x_y(1,1,1))
        allocate(sheet_sig_y_y(1,1,1))
        allocate(sheet_sig_z_y(1,1,1))
    end if
    if (num_sheets_z>0) then
        !material cells
        allocate(sheet_ep_x_cell_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_cell_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_cell_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_cell_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_cell_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_cell_z(x_size-1,y_size-1,z_size-1))
        !grid points
        allocate(sheet_ep_x_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_y_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_ep_z_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_x_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_y_z(x_size-1,y_size-1,z_size-1))
        allocate(sheet_sig_z_z(x_size-1,y_size-1,z_size-1))
    else
        !material cells
        allocate(sheet_ep_x_cell_z(1,1,1))
        allocate(sheet_ep_y_cell_z(1,1,1))
        allocate(sheet_ep_z_cell_z(1,1,1))
        allocate(sheet_sig_x_cell_z(1,1,1))
        allocate(sheet_sig_y_cell_z(1,1,1))
        allocate(sheet_sig_z_cell_z(1,1,1))
        !grid points
        allocate(sheet_ep_x_z(1,1,1))
        allocate(sheet_ep_y_z(1,1,1))
        allocate(sheet_ep_z_z(1,1,1))
        allocate(sheet_sig_x_z(1,1,1))
        allocate(sheet_sig_y_z(1,1,1))
        allocate(sheet_sig_z_z(1,1,1))
    end if

    !added specifically for plasmas
    if (is_plasma==1) then
        ! X-component allocations
        allocate(plasma_freq_cell_x(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_freq_x(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_cell_x(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_x(num_poles,x_size-1,y_size-1,z_size-1))
        ! Y-component allocations
        allocate(plasma_freq_cell_y(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_freq_y(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_cell_y(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_y(num_poles,x_size-1,y_size-1,z_size-1))
        ! Z-component allocations
        allocate(plasma_freq_cell_z(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_freq_z(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_cell_z(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(plasma_loss_z(num_poles,x_size-1,y_size-1,z_size-1))
    else
        allocate(plasma_freq_cell_x(1,1,1,1), plasma_freq_x(1,1,1,1))
        allocate(plasma_loss_cell_x(1,1,1,1), plasma_loss_x(1,1,1,1))
        allocate(plasma_freq_cell_y(1,1,1,1), plasma_freq_y(1,1,1,1))
        allocate(plasma_loss_cell_y(1,1,1,1), plasma_loss_y(1,1,1,1))
        allocate(plasma_freq_cell_z(1,1,1,1), plasma_freq_z(1,1,1,1))
        allocate(plasma_loss_cell_z(1,1,1,1), plasma_loss_z(1,1,1,1))
    end if

    !for output geometry file
    allocate(output_geometry(4,x_size-1-pbc_x,y_size-1-pbc_y,z_size-1-pbc_z)) 

    !initialize all arrays

    !material cells first
    relative_ep_x_cell(:,:,:)=1.0
    sigma_x_cell(:,:,:)=0.0
    relative_ep_y_cell(:,:,:)=1.0
    sigma_y_cell(:,:,:)=0.0
    relative_ep_z_cell(:,:,:)=1.0
    sigma_z_cell(:,:,:)=0.0
    !x directed sheet material cells
    sheet_ep_x_cell_x(:,:,:)=1.0
    sheet_ep_y_cell_x(:,:,:)=1.0
    sheet_ep_z_cell_x(:,:,:)=1.0
    sheet_sig_x_cell_x(:,:,:)=0.0
    sheet_sig_y_cell_x(:,:,:)=0.0
    sheet_sig_z_cell_x(:,:,:)=0.0
    !y directed sheet material cells
    sheet_ep_x_cell_y(:,:,:)=1.0
    sheet_ep_y_cell_y(:,:,:)=1.0
    sheet_ep_z_cell_y(:,:,:)=1.0
    sheet_sig_x_cell_y(:,:,:)=0.0
    sheet_sig_y_cell_y(:,:,:)=0.0
    sheet_sig_z_cell_y(:,:,:)=0.0
    !z directed sheet material cells
    sheet_ep_x_cell_z(:,:,:)=1.0
    sheet_ep_y_cell_z(:,:,:)=1.0
    sheet_ep_z_cell_z(:,:,:)=1.0
    sheet_sig_x_cell_z(:,:,:)=0.0
    sheet_sig_y_cell_z(:,:,:)=0.0
    sheet_sig_z_cell_z(:,:,:)=0.0

    !now yee grid points
    relative_ep_x(:,:,:)=1.0
    sigma_x(:,:,:)=0.0
    relative_ep_y(:,:,:)=1.0
    sigma_y(:,:,:)=0.0
    relative_ep_z(:,:,:)=1.0
    sigma_z(:,:,:)=0.0
    !x directed sheet grid points
    sheet_ep_x_x(:,:,:)=1.0
    sheet_ep_y_x(:,:,:)=1.0
    sheet_ep_z_x(:,:,:)=1.0
    sheet_sig_x_x(:,:,:)=0.0
    sheet_sig_y_x(:,:,:)=0.0
    sheet_sig_z_x(:,:,:)=0.0
    !y directed sheet grid points
    sheet_ep_x_y(:,:,:)=1.0
    sheet_ep_y_y(:,:,:)=1.0
    sheet_ep_z_y(:,:,:)=1.0
    sheet_sig_x_y(:,:,:)=0.0
    sheet_sig_y_y(:,:,:)=0.0
    sheet_sig_z_y(:,:,:)=0.0
    !z directed sheet grid points
    sheet_ep_x_z(:,:,:)=1.0
    sheet_ep_y_z(:,:,:)=1.0
    sheet_ep_z_z(:,:,:)=1.0
    sheet_sig_x_z(:,:,:)=0.0
    sheet_sig_y_z(:,:,:)=0.0
    sheet_sig_z_z(:,:,:)=0.0

    ! Specifically for plasmas - both material and yee grid points
    ! Initialize X-component arrays
    plasma_freq_cell_x(:,:,:,:) = 0.0
    plasma_freq_x(:,:,:,:) = 0.0
    plasma_loss_cell_x(:,:,:,:) = 0.0
    plasma_loss_x(:,:,:,:) = 0.0
    ! Initialize Y-component arrays
    plasma_freq_cell_y(:,:,:,:) = 0.0
    plasma_freq_y(:,:,:,:) = 0.0
    plasma_loss_cell_y(:,:,:,:) = 0.0
    plasma_loss_y(:,:,:,:) = 0.0
    ! Initialize Z-component arrays
    plasma_freq_cell_z(:,:,:,:) = 0.0
    plasma_freq_z(:,:,:,:) = 0.0
    plasma_loss_cell_z(:,:,:,:) = 0.0
    plasma_loss_z(:,:,:,:) = 0.0

    !initialize output geometry to zeros (material_ID=0)
    output_geometry(:,:,:,:)=0.0

    !need to build each block, sphere, cylinder, and sheet with the given info
    !this first section is for optional read-in of bulk geometry only for now
    !any geometry drawn will overwrite the optional read-in file
    !the drawing parts are almost identical to the normal parts below with a small difference
    
    !if it exists then do all the things (allocate, read, and create geometry)
    if (trim(use_optional_bulk_geom_file)=='yes') then
        allocate(optional_bulk_geom(x_size-1-pbc_x,y_size-1-pbc_y,z_size-1-pbc_z))       
        open(unit=82, file=trim(optional_bulk_geom_filename), access='stream', status='old')
            read(82) optional_bulk_geom
        close(82)
        allocate(temp_mat_id(0:maxval(int(materials_properties(:, 1, 1)))))
        temp_mat_id(:)=0
        do rr=1, num_materials
            temp_mat_id(int(materials_properties(rr, 1, 1))) = rr
        end do
        do kk=1, z_size-1-pbc_z
            do jj=1, y_size-1-pbc_y
                do i=1, x_size-1-pbc_x
                    !first read the material ID number
                    materials_id = temp_mat_id(int(optional_bulk_geom(i,jj,kk)))
                    if (materials_id /= 0) then
                        !the rest is the same as the other normal section copied over
                        relative_ep_x_cell(i,jj,kk)=materials_properties(materials_id,2,1)
                        relative_ep_y_cell(i,jj,kk)=materials_properties(materials_id,2,2)
                        relative_ep_z_cell(i,jj,kk)=materials_properties(materials_id,2,3)
                        sigma_x_cell(i,jj,kk)=materials_properties(materials_id,3,1)
                        sigma_y_cell(i,jj,kk)=materials_properties(materials_id,3,2)
                        sigma_z_cell(i,jj,kk)=materials_properties(materials_id,3,3)
                        if (is_plasma==1) then 
                            !this doesn't fiter by material, just in general. 
                            !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                            !Just saves some time.
                            !Thus even if all poles aren't used, that's okay, they will yield zero here.
                            do plasma_counter=1,num_poles
                                plasma_freq_cell_x(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+2,1)
                                plasma_freq_cell_y(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+2,2)
                                plasma_freq_cell_z(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+2,3)
                                plasma_loss_cell_x(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+3,1)
                                plasma_loss_cell_y(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+3,2)
                                plasma_loss_cell_z(plasma_counter,i,jj,kk)=materials_properties(materials_id,2*plasma_counter+3,3)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        !now add to output geometry file
        output_geometry(1,:,:,:)=optional_bulk_geom
        deallocate(temp_mat_id)
    end if  

    !now the normal read-in section
    ii=1 !blocks
    jj=1 !spheres
    kk=1 !cylinders
    do counter=1, num_objects
        select case(trim(object_type(counter)))

        case('block')
            !add block properties to the grid
            do rr=1, num_materials
                if (blocks(ii,1,1)==materials_properties(rr,1,1)) then
                    materials_id = rr
                end if
            end do
            do k=(blocks(ii,2,3)),(blocks(ii,2,3)+blocks(ii,3,3)-1)
                do j=(blocks(ii,2,2)),(blocks(ii,2,2)+blocks(ii,3,2)-1)
                    do i=blocks(ii,2,1),(blocks(ii,2,1)+blocks(ii,3,1)-1)
                        relative_ep_x_cell(i,j,k)=materials_properties(materials_id,2,1)
                        relative_ep_y_cell(i,j,k)=materials_properties(materials_id,2,2)
                        relative_ep_z_cell(i,j,k)=materials_properties(materials_id,2,3)
                        sigma_x_cell(i,j,k)=materials_properties(materials_id,3,1)
                        sigma_y_cell(i,j,k)=materials_properties(materials_id,3,2)
                        sigma_z_cell(i,j,k)=materials_properties(materials_id,3,3)
                        if (is_plasma==1) then 
                            !this doesn't fiter by material, just in general. 
                            !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                            !Just saves some time.
                            !Thus even if all poles aren't used, that's okay, they will yield zero here.
                            do plasma_counter=1,num_poles
                                plasma_freq_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,1)
                                plasma_freq_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,2)
                                plasma_freq_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,3)
                                plasma_loss_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,1)
                                plasma_loss_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,2)
                                plasma_loss_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,3)
                            end do
                        end if
                        !now add to geometry file
                        output_geometry(1,i,j,k)=blocks(ii,1,1)
                    end do
                end do
            end do
            ii=ii+1

        case('sphere')
            !add sphere properties to the grid
            do rr=1, num_materials
                if (spheres(jj,1,1)==materials_properties(rr,1,1)) then
                    materials_id = rr
                end if
            end do
            do k=(spheres(jj,2,3)-spheres(jj,3,1)),(spheres(jj,2,3)+&
            spheres(jj,3,1)-1)
                do j=(spheres(jj,2,2)-spheres(jj,3,1)),(spheres(jj,2,2)+&
                spheres(jj,3,1)-1)
                    do i=(spheres(jj,2,1)-spheres(jj,3,1)),(spheres(jj,2,1)+spheres(jj,3,1)-1)
                        if (((real(i)-real(spheres(jj,2,1))+0.5)**2+(real(j)-real(spheres(jj,2,2))+0.5)**2 &
                        +(real(k)-real(spheres(jj,2,3))+0.5)**2)<=real(spheres(jj,3,1))**2) then
                            relative_ep_x_cell(i,j,k)=materials_properties(materials_id,2,1)
                            relative_ep_y_cell(i,j,k)=materials_properties(materials_id,2,2)
                            relative_ep_z_cell(i,j,k)=materials_properties(materials_id,2,3)
                            sigma_x_cell(i,j,k)=materials_properties(materials_id,3,1)
                            sigma_y_cell(i,j,k)=materials_properties(materials_id,3,2)
                            sigma_z_cell(i,j,k)=materials_properties(materials_id,3,3)
                            if (is_plasma==1) then 
                                !this doesn't fiter by material, just in general. 
                                !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                                !Just saves some time.
                                !Thus even if all poles aren't used, that's okay, they will yield zero here.
                                do plasma_counter=1,num_poles
                                    plasma_freq_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,1)
                                    plasma_freq_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,2)
                                    plasma_freq_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,3)
                                    plasma_loss_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,1)
                                    plasma_loss_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,2)
                                    plasma_loss_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,3)
                                end do
                            end if
                            !now add to geometry file
                            output_geometry(1,i,j,k)=spheres(jj,1,1)
                        end if
                    end do
                end do
            end do
            jj=jj+1

        case('cylinder')
            do rr=1, num_materials
                if (cylinders(kk,1,1)==materials_properties(rr,1,1)) then
                    materials_id = rr
                end if
            end do
            if (cylinders(kk,2,1)==0) then
                do k=(cylinders(kk,3,3)-cylinders(kk,4,2)),(cylinders(kk,3,3)+cylinders(kk,4,2)-1)
                    do j=(cylinders(kk,3,2)-cylinders(kk,4,2)),(cylinders(kk,3,2)+cylinders(kk,4,2)-1)
                        do i=cylinders(kk,3,1),(cylinders(kk,3,1)+cylinders(kk,4,1)-1)
                            if (((real(j)-real(cylinders(kk,3,2))+0.5)**2+(real(k)-real(cylinders(kk,3,3))+0.5)**2)<=real(cylinders(kk,4,2))**2) then
                                relative_ep_x_cell(i,j,k)=materials_properties(materials_id,2,1)
                                relative_ep_y_cell(i,j,k)=materials_properties(materials_id,2,2)
                                relative_ep_z_cell(i,j,k)=materials_properties(materials_id,2,3)
                                sigma_x_cell(i,j,k)=materials_properties(materials_id,3,1)
                                sigma_y_cell(i,j,k)=materials_properties(materials_id,3,2)
                                sigma_z_cell(i,j,k)=materials_properties(materials_id,3,3)
                                if (is_plasma==1) then 
                                    !this doesn't fiter by material, just in general. 
                                    !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                                    !Just saves some time.
                                    !Thus even if all poles aren't used, that's okay, they will yield zero here.
                                    do plasma_counter=1,num_poles
                                        plasma_freq_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,1)
                                        plasma_freq_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,2)
                                        plasma_freq_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,3)
                                        plasma_loss_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,1)
                                        plasma_loss_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,2)
                                        plasma_loss_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,3)
                                    end do
                                end if
                                !now add to geometry file
                                output_geometry(1,i,j,k)=cylinders(kk,1,1)
                            end if
                        end do
                    end do
                end do
            end if
            if (cylinders(kk,2,1)==1) then
                do k=(cylinders(kk,3,3)-cylinders(kk,4,2)),(cylinders(kk,3,3)+cylinders(kk,4,2)-1)
                    do j=(cylinders(kk,3,2)),(cylinders(kk,3,2)+cylinders(kk,4,1)-1)
                        do i=(cylinders(kk,3,1)-cylinders(kk,4,2)),(cylinders(kk,3,1)+cylinders(kk,4,2)-1)
                            if (((real(i)-real(cylinders(kk,3,1))+0.5)**2+(real(k)-real(cylinders(kk,3,3))+0.5)**2)<=real(cylinders(kk,4,2))**2) then
                                relative_ep_x_cell(i,j,k)=materials_properties(materials_id,2,1)
                                relative_ep_y_cell(i,j,k)=materials_properties(materials_id,2,2)
                                relative_ep_z_cell(i,j,k)=materials_properties(materials_id,2,3)
                                sigma_x_cell(i,j,k)=materials_properties(materials_id,3,1)
                                sigma_y_cell(i,j,k)=materials_properties(materials_id,3,2)
                                sigma_z_cell(i,j,k)=materials_properties(materials_id,3,3)
                                if (is_plasma==1) then 
                                    !this doesn't fiter by material, just in general. 
                                    !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                                    !Just saves some time.
                                    !Thus even if all poles aren't used, that's okay, they will yield zero here.
                                    do plasma_counter=1,num_poles
                                        plasma_freq_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,1)
                                        plasma_freq_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,2)
                                        plasma_freq_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,3)
                                        plasma_loss_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,1)
                                        plasma_loss_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,2)
                                        plasma_loss_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,3)
                                    end do
                                end if
                                !now add to geometry file
                                output_geometry(1,i,j,k)=cylinders(kk,1,1)
                            end if
                        end do
                    end do
                end do
            end if
            if (cylinders(kk,2,1)==2) then
                do k=(cylinders(kk,3,3)),(cylinders(kk,3,3)+cylinders(kk,4,1)-1)
                    do j=(cylinders(kk,3,2)-cylinders(kk,4,2)),(cylinders(kk,3,2)+cylinders(kk,4,2)-1)
                        do i=(cylinders(kk,3,1)-cylinders(kk,4,2)),(cylinders(kk,3,1)+cylinders(kk,4,2)-1)
                            if (((real(i)-real(cylinders(kk,3,1))+0.5)**2+(real(j)-real(cylinders(kk,3,2))+0.5)**2)<=real(cylinders(kk,4,2))**2) then
                                relative_ep_x_cell(i,j,k)=materials_properties(materials_id,2,1)
                                relative_ep_y_cell(i,j,k)=materials_properties(materials_id,2,2)
                                relative_ep_z_cell(i,j,k)=materials_properties(materials_id,2,3)
                                sigma_x_cell(i,j,k)=materials_properties(materials_id,3,1)
                                sigma_y_cell(i,j,k)=materials_properties(materials_id,3,2)
                                sigma_z_cell(i,j,k)=materials_properties(materials_id,3,3)
                                if (is_plasma==1) then 
                                    !this doesn't fiter by material, just in general. 
                                    !All materials will default these properties to zero anyway if not plasma and filled in accordingly. 
                                    !Just saves some time.
                                    !Thus even if all poles aren't used, that's okay, they will yield zero here.
                                    do plasma_counter=1,num_poles
                                        plasma_freq_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,1)
                                        plasma_freq_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,2)
                                        plasma_freq_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+2,3)
                                        plasma_loss_cell_x(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,1)
                                        plasma_loss_cell_y(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,2)
                                        plasma_loss_cell_z(plasma_counter,i,j,k)=materials_properties(materials_id,2*plasma_counter+3,3)
                                    end do
                                end if
                                !now add to geometry file
                                output_geometry(1,i,j,k)=cylinders(kk,1,1)
                            end if
                        end do
                    end do
                end do
            end if
            kk=kk+1

        end select
    end do

    do ii=1, num_sheets_x
        !loop over sheet size and use material properties
        do rr=1, num_sheet_materials
            if (sheets_x(ii,1,1)==sheet_properties(rr,1,1)) then
                sheet_materials_id = rr
            end if
        end do
        i=sheets_x(ii,2,1) !sheet location in x
        do k=(sheets_x(ii,3,2)),(sheets_x(ii,3,2)+sheets_x(ii,4,2)-1)
            do j=(sheets_x(ii,3,1)),(sheets_x(ii,3,1)+sheets_x(ii,4,1)-1)
                sheet_ep_x_cell_x(i,j,k)=sheet_properties(sheet_materials_id,3,1)
                sheet_ep_y_cell_x(i,j,k)=sheet_properties(sheet_materials_id,3,2)
                sheet_ep_z_cell_x(i,j,k)=sheet_properties(sheet_materials_id,3,3)
                sheet_sig_x_cell_x(i,j,k)=sheet_properties(sheet_materials_id,4,1)
                sheet_sig_y_cell_x(i,j,k)=sheet_properties(sheet_materials_id,4,2)
                sheet_sig_z_cell_x(i,j,k)=sheet_properties(sheet_materials_id,4,3)
                !now add to geometry file
                output_geometry(2,i,j,k)=sheets_x(ii,1,1)
            end do
        end do
    end do
    do ii=1, num_sheets_y
        !loop over sheet size and use material properties
        do rr=1, num_sheet_materials
            if (sheets_y(ii,1,1)==sheet_properties(rr,1,1)) then
                sheet_materials_id = rr
            end if
        end do
        j=sheets_y(ii,2,1) !sheet location in y
        do k=(sheets_y(ii,3,2)),(sheets_y(ii,3,2)+sheets_y(ii,4,2)-1)
            do i=(sheets_y(ii,3,1)),(sheets_y(ii,3,1)+sheets_y(ii,4,1)-1)
                sheet_ep_x_cell_y(i,j,k)=sheet_properties(sheet_materials_id,3,1)
                sheet_ep_y_cell_y(i,j,k)=sheet_properties(sheet_materials_id,3,2)
                sheet_ep_z_cell_y(i,j,k)=sheet_properties(sheet_materials_id,3,3)
                sheet_sig_x_cell_y(i,j,k)=sheet_properties(sheet_materials_id,4,1)
                sheet_sig_y_cell_y(i,j,k)=sheet_properties(sheet_materials_id,4,2)
                sheet_sig_z_cell_y(i,j,k)=sheet_properties(sheet_materials_id,4,3)
                !now add to geometry file
                output_geometry(3,i,j,k)=sheets_y(ii,1,1)
            end do
        end do
    end do
    do ii=1, num_sheets_z
        !loop over sheet size and use material properties
        do rr=1, num_sheet_materials
            if (sheets_z(ii,1,1)==sheet_properties(rr,1,1)) then
                sheet_materials_id = rr
            end if
        end do
        k=sheets_z(ii,2,1) !sheet location in z
        do j=(sheets_z(ii,3,2)),(sheets_z(ii,3,2)+sheets_z(ii,4,2)-1)          
            do i=(sheets_z(ii,3,1)),(sheets_z(ii,3,1)+sheets_z(ii,4,1)-1)
                sheet_ep_x_cell_z(i,j,k)=sheet_properties(sheet_materials_id,3,1)
                sheet_ep_y_cell_z(i,j,k)=sheet_properties(sheet_materials_id,3,2)
                sheet_ep_z_cell_z(i,j,k)=sheet_properties(sheet_materials_id,3,3)
                sheet_sig_x_cell_z(i,j,k)=sheet_properties(sheet_materials_id,4,1)
                sheet_sig_y_cell_z(i,j,k)=sheet_properties(sheet_materials_id,4,2)
                sheet_sig_z_cell_z(i,j,k)=sheet_properties(sheet_materials_id,4,3)
                !now add to geometry file
                output_geometry(4,i,j,k)=sheets_z(ii,1,1)
            end do
        end do
    end do

    ! save geometry before continuing on
    open(unit=83, file=filename(1:LEN_TRIM(filename)-4)//"_"//"geometry.bin", access='stream', status='replace')
        write(83) output_geometry
    close(83)

    !before yee cell creation from material cells - need to create a repeated materials cell layer into size-1 layers if pbc
    if (pbc_x==1) then
        do k=1,z_size-1
            do j=1,y_size-1
                !blocks
                relative_ep_x_cell(x_size-1,j,k)=relative_ep_x_cell(1,j,k)
                relative_ep_y_cell(x_size-1,j,k)=relative_ep_y_cell(1,j,k)
                relative_ep_z_cell(x_size-1,j,k)=relative_ep_z_cell(1,j,k)
                sigma_x_cell(x_size-1,j,k)=sigma_x_cell(1,j,k)
                sigma_y_cell(x_size-1,j,k)=sigma_y_cell(1,j,k)
                sigma_z_cell(x_size-1,j,k)=sigma_z_cell(1,j,k)                
                !sheets x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_cell_x(x_size-1,j,k)=sheet_ep_x_cell_x(1,j,k)
                    sheet_ep_y_cell_x(x_size-1,j,k)=sheet_ep_y_cell_x(1,j,k)
                    sheet_ep_z_cell_x(x_size-1,j,k)=sheet_ep_z_cell_x(1,j,k)
                    sheet_sig_x_cell_x(x_size-1,j,k)=sheet_sig_x_cell_x(1,j,k)
                    sheet_sig_y_cell_x(x_size-1,j,k)=sheet_sig_y_cell_x(1,j,k)
                    sheet_sig_z_cell_x(x_size-1,j,k)=sheet_sig_z_cell_x(1,j,k)
                end if
                !sheets y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_cell_y(x_size-1,j,k)=sheet_ep_x_cell_y(1,j,k)
                    sheet_ep_y_cell_y(x_size-1,j,k)=sheet_ep_y_cell_y(1,j,k)
                    sheet_ep_z_cell_y(x_size-1,j,k)=sheet_ep_z_cell_y(1,j,k)
                    sheet_sig_x_cell_y(x_size-1,j,k)=sheet_sig_x_cell_y(1,j,k)
                    sheet_sig_y_cell_y(x_size-1,j,k)=sheet_sig_y_cell_y(1,j,k)
                    sheet_sig_z_cell_y(x_size-1,j,k)=sheet_sig_z_cell_y(1,j,k)
                end if
                !sheets z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_cell_z(x_size-1,j,k)=sheet_ep_x_cell_z(1,j,k)
                    sheet_ep_y_cell_z(x_size-1,j,k)=sheet_ep_y_cell_z(1,j,k)
                    sheet_ep_z_cell_z(x_size-1,j,k)=sheet_ep_z_cell_z(1,j,k)
                    sheet_sig_x_cell_z(x_size-1,j,k)=sheet_sig_x_cell_z(1,j,k)
                    sheet_sig_y_cell_z(x_size-1,j,k)=sheet_sig_y_cell_z(1,j,k)
                    sheet_sig_z_cell_z(x_size-1,j,k)=sheet_sig_z_cell_z(1,j,k)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_cell_x(plasma_counter,x_size-1,j,k)=plasma_freq_cell_x(plasma_counter,1,j,k)
                        plasma_loss_cell_x(plasma_counter,x_size-1,j,k)=plasma_loss_cell_x(plasma_counter,1,j,k)
                        plasma_freq_cell_y(plasma_counter,x_size-1,j,k)=plasma_freq_cell_y(plasma_counter,1,j,k)
                        plasma_loss_cell_y(plasma_counter,x_size-1,j,k)=plasma_loss_cell_y(plasma_counter,1,j,k)
                        plasma_freq_cell_z(plasma_counter,x_size-1,j,k)=plasma_freq_cell_z(plasma_counter,1,j,k)
                        plasma_loss_cell_z(plasma_counter,x_size-1,j,k)=plasma_loss_cell_z(plasma_counter,1,j,k)
                    end do
                end if
            end do
        end do
    end if
    if (pbc_y==1) then
        do k=1,z_size-1
            do i=1,x_size-1
                relative_ep_x_cell(i,y_size-1,k)=relative_ep_x_cell(i,1,k)
                relative_ep_y_cell(i,y_size-1,k)=relative_ep_y_cell(i,1,k)
                relative_ep_z_cell(i,y_size-1,k)=relative_ep_z_cell(i,1,k)
                sigma_x_cell(i,y_size-1,k)=sigma_x_cell(i,1,k)
                sigma_y_cell(i,y_size-1,k)=sigma_y_cell(i,1,k)
                sigma_z_cell(i,y_size-1,k)=sigma_z_cell(i,1,k)
                !sheets x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_cell_x(i,y_size-1,k)=sheet_ep_x_cell_x(i,1,k)
                    sheet_ep_y_cell_x(i,y_size-1,k)=sheet_ep_y_cell_x(i,1,k)
                    sheet_ep_z_cell_x(i,y_size-1,k)=sheet_ep_z_cell_x(i,1,k)
                    sheet_sig_x_cell_x(i,y_size-1,k)=sheet_sig_x_cell_x(i,1,k)
                    sheet_sig_y_cell_x(i,y_size-1,k)=sheet_sig_y_cell_x(i,1,k)
                    sheet_sig_z_cell_x(i,y_size-1,k)=sheet_sig_z_cell_x(i,1,k)
                end if
                !sheets y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_cell_y(i,y_size-1,k)=sheet_ep_x_cell_y(i,1,k)
                    sheet_ep_y_cell_y(i,y_size-1,k)=sheet_ep_y_cell_y(i,1,k)
                    sheet_ep_z_cell_y(i,y_size-1,k)=sheet_ep_z_cell_y(i,1,k)
                    sheet_sig_x_cell_y(i,y_size-1,k)=sheet_sig_x_cell_y(i,1,k)
                    sheet_sig_y_cell_y(i,y_size-1,k)=sheet_sig_y_cell_y(i,1,k)
                    sheet_sig_z_cell_y(i,y_size-1,k)=sheet_sig_z_cell_y(i,1,k)
                end if
                !sheets z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_cell_z(i,y_size-1,k)=sheet_ep_x_cell_z(i,1,k)
                    sheet_ep_y_cell_z(i,y_size-1,k)=sheet_ep_y_cell_z(i,1,k)
                    sheet_ep_z_cell_z(i,y_size-1,k)=sheet_ep_z_cell_z(i,1,k)
                    sheet_sig_x_cell_z(i,y_size-1,k)=sheet_sig_x_cell_z(i,1,k)
                    sheet_sig_y_cell_z(i,y_size-1,k)=sheet_sig_y_cell_z(i,1,k)
                    sheet_sig_z_cell_z(i,y_size-1,k)=sheet_sig_z_cell_z(i,1,k)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_cell_x(plasma_counter,i,y_size-1,k)=plasma_freq_cell_x(plasma_counter,i,1,k)
                        plasma_loss_cell_x(plasma_counter,i,y_size-1,k)=plasma_loss_cell_x(plasma_counter,i,1,k)
                        plasma_freq_cell_y(plasma_counter,i,y_size-1,k)=plasma_freq_cell_y(plasma_counter,i,1,k)
                        plasma_loss_cell_y(plasma_counter,i,y_size-1,k)=plasma_loss_cell_y(plasma_counter,i,1,k)
                        plasma_freq_cell_z(plasma_counter,i,y_size-1,k)=plasma_freq_cell_z(plasma_counter,i,1,k)
                        plasma_loss_cell_z(plasma_counter,i,y_size-1,k)=plasma_loss_cell_z(plasma_counter,i,1,k)
                    end do
                end if 
            end do
        end do
    end if
    if (pbc_z==1) then
        do j=1,y_size-1
            do i=1,x_size-1
                !blocks
                relative_ep_x_cell(i,j,z_size-1)=relative_ep_x_cell(i,j,1)
                relative_ep_y_cell(i,j,z_size-1)=relative_ep_y_cell(i,j,1)
                relative_ep_z_cell(i,j,z_size-1)=relative_ep_z_cell(i,j,1)
                sigma_x_cell(i,j,z_size-1)=sigma_x_cell(i,j,1)
                sigma_y_cell(i,j,z_size-1)=sigma_y_cell(i,j,1)
                sigma_z_cell(i,j,z_size-1)=sigma_z_cell(i,j,1)
                !sheets x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_cell_x(i,j,z_size-1)=sheet_ep_x_cell_x(i,j,1)
                    sheet_ep_y_cell_x(i,j,z_size-1)=sheet_ep_y_cell_x(i,j,1)
                    sheet_ep_z_cell_x(i,j,z_size-1)=sheet_ep_z_cell_x(i,j,1)
                    sheet_sig_x_cell_x(i,j,z_size-1)=sheet_sig_x_cell_x(i,j,1)
                    sheet_sig_y_cell_x(i,j,z_size-1)=sheet_sig_y_cell_x(i,j,1)
                    sheet_sig_z_cell_x(i,j,z_size-1)=sheet_sig_z_cell_x(i,j,1)
                end if
                !sheets y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_cell_y(i,j,z_size-1)=sheet_ep_x_cell_y(i,j,1)
                    sheet_ep_y_cell_y(i,j,z_size-1)=sheet_ep_y_cell_y(i,j,1)
                    sheet_ep_z_cell_y(i,j,z_size-1)=sheet_ep_z_cell_y(i,j,1)
                    sheet_sig_x_cell_y(i,j,z_size-1)=sheet_sig_x_cell_y(i,j,1)
                    sheet_sig_y_cell_y(i,j,z_size-1)=sheet_sig_y_cell_y(i,j,1)
                    sheet_sig_z_cell_y(i,j,z_size-1)=sheet_sig_z_cell_y(i,j,1)
                end if
                !sheets z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_cell_z(i,j,z_size-1)=sheet_ep_x_cell_z(i,j,1)
                    sheet_ep_y_cell_z(i,j,z_size-1)=sheet_ep_y_cell_z(i,j,1)
                    sheet_ep_z_cell_z(i,j,z_size-1)=sheet_ep_z_cell_z(i,j,1)
                    sheet_sig_x_cell_z(i,j,z_size-1)=sheet_sig_x_cell_z(i,j,1)
                    sheet_sig_y_cell_z(i,j,z_size-1)=sheet_sig_y_cell_z(i,j,1)
                    sheet_sig_z_cell_z(i,j,z_size-1)=sheet_sig_z_cell_z(i,j,1)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_cell_x(plasma_counter,i,j,z_size-1)=plasma_freq_cell_x(plasma_counter,i,j,1)
                        plasma_loss_cell_x(plasma_counter,i,j,z_size-1)=plasma_loss_cell_x(plasma_counter,i,j,1)
                        plasma_freq_cell_y(plasma_counter,i,j,z_size-1)=plasma_freq_cell_y(plasma_counter,i,j,1)
                        plasma_loss_cell_y(plasma_counter,i,j,z_size-1)=plasma_loss_cell_y(plasma_counter,i,j,1)
                        plasma_freq_cell_z(plasma_counter,i,j,z_size-1)=plasma_freq_cell_z(plasma_counter,i,j,1)
                        plasma_loss_cell_z(plasma_counter,i,j,z_size-1)=plasma_loss_cell_z(plasma_counter,i,j,1)
                    end do
                end if
            end do
        end do
    end if

    !now I need to convert the material cells and sheets to yee grid points

    !bulk materials first
    do k=2, z_size-1
        do j=2, y_size-1
            do i=2,x_size-1
                relative_ep_x(i,j,k)=(relative_ep_x_cell(i,j,k)+relative_ep_x_cell(i,j-1,k)+&
                relative_ep_x_cell(i,j,k-1)+relative_ep_x_cell(i,j-1,k-1))/4.0
                relative_ep_y(i,j,k)=(relative_ep_y_cell(i,j,k)+relative_ep_y_cell(i-1,j,k)+&
                relative_ep_y_cell(i,j,k-1)+relative_ep_y_cell(i-1,j,k-1))/4.0
                relative_ep_z(i,j,k)=(relative_ep_z_cell(i,j,k)+relative_ep_z_cell(i-1,j,k)+&
                relative_ep_z_cell(i,j-1,k)+relative_ep_z_cell(i-1,j-1,k))/4.0
                sigma_x(i,j,k)=(sigma_x_cell(i,j,k)+sigma_x_cell(i,j-1,k)+sigma_x_cell(i,j,k-1)&
                +sigma_x_cell(i,j-1,k-1))/4.0
                sigma_y(i,j,k)=(sigma_y_cell(i,j,k)+sigma_y_cell(i-1,j,k)+sigma_y_cell(i,j,k-1)&
                +sigma_y_cell(i-1,j,k-1))/4.0
                sigma_z(i,j,k)=(sigma_z_cell(i,j,k)+sigma_z_cell(i-1,j,k)+sigma_z_cell(i,j-1,k)&
                +sigma_z_cell(i-1,j-1,k))/4.0
            end do
        end do
    end do
    !now sheets
    do k=2,z_size-1
        do j=2,y_size-1
            do i=2,x_size-1
                !x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_x(i,j,k)=(sheet_ep_x_cell_x(i,j,k)+sheet_ep_x_cell_x(i,j-1,k)+sheet_ep_x_cell_x(i,j,k-1)&
                    +sheet_ep_x_cell_x(i,j-1,k-1))/4.0
                    sheet_ep_y_x(i,j,k)=(sheet_ep_y_cell_x(i,j,k)+sheet_ep_y_cell_x(i,j,k-1)+sheet_ep_y_cell_x(i,j,k)&
                    +sheet_ep_y_cell_x(i,j,k-1))/4.0
                    sheet_ep_z_x(i,j,k)=(sheet_ep_z_cell_x(i,j,k)+sheet_ep_z_cell_x(i,j-1,k)+sheet_ep_z_cell_x(i,j,k)&
                    +sheet_ep_z_cell_x(i,j-1,k))/4.0
                    sheet_sig_x_x(i,j,k)=(sheet_sig_x_cell_x(i,j,k)+sheet_sig_x_cell_x(i,j-1,k)+sheet_sig_x_cell_x(i,j,k-1)+&
                    sheet_sig_x_cell_x(i,j-1,k-1))/4.0
                    sheet_sig_y_x(i,j,k)=(sheet_sig_y_cell_x(i,j,k)+sheet_sig_y_cell_x(i,j,k-1)+sheet_sig_y_cell_x(i,j,k)+&
                    sheet_sig_y_cell_x(i,j,k-1))/4.0
                    sheet_sig_z_x(i,j,k)=(sheet_sig_z_cell_x(i,j,k)+sheet_sig_z_cell_x(i,j-1,k)+sheet_sig_z_cell_x(i,j,k)+&
                    sheet_sig_z_cell_x(i,j-1,k))/4.0
                end if
                !y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_y(i,j,k)=(sheet_ep_x_cell_y(i,j,k)+sheet_ep_x_cell_y(i,j,k)+sheet_ep_x_cell_y(i,j,k-1)&
                    +sheet_ep_x_cell_y(i,j,k-1))/4.0
                    sheet_ep_y_y(i,j,k)=(sheet_ep_y_cell_y(i,j,k)+sheet_ep_y_cell_y(i,j,k-1)+sheet_ep_y_cell_y(i-1,j,k)&
                    +sheet_ep_y_cell_y(i-1,j,k-1))/4.0
                    sheet_ep_z_y(i,j,k)=(sheet_ep_z_cell_y(i,j,k)+sheet_ep_z_cell_y(i-1,j,k)+sheet_ep_z_cell_y(i,j,k)&
                    +sheet_ep_z_cell_y(i-1,j,k))/4.0
                    sheet_sig_x_y(i,j,k)=(sheet_sig_x_cell_y(i,j,k)+sheet_sig_x_cell_y(i,j,k)+sheet_sig_x_cell_y(i,j,k-1)+&
                    sheet_sig_x_cell_y(i,j,k-1))/4.0
                    sheet_sig_y_y(i,j,k)=(sheet_sig_y_cell_y(i,j,k)+sheet_sig_y_cell_y(i,j,k-1)+sheet_sig_y_cell_y(i-1,j,k)+&
                    sheet_sig_y_cell_y(i-1,j,k-1))/4.0
                    sheet_sig_z_y(i,j,k)=(sheet_sig_z_cell_y(i,j,k)+sheet_sig_z_cell_y(i-1,j,k)+sheet_sig_z_cell_y(i,j,k)+&
                    sheet_sig_z_cell_y(i-1,j,k))/4.0
                end if
                !z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_z(i,j,k)=(sheet_ep_x_cell_z(i,j,k)+sheet_ep_x_cell_z(i,j-1,k)+sheet_ep_x_cell_z(i,j,k)&
                    +sheet_ep_x_cell_z(i,j-1,k))/4.0
                    sheet_ep_y_z(i,j,k)=(sheet_ep_y_cell_z(i,j,k)+sheet_ep_y_cell_z(i-1,j,k)+sheet_ep_y_cell_z(i,j,k)&
                    +sheet_ep_y_cell_z(i-1,j,k))/4.0
                    sheet_ep_z_z(i,j,k)=(sheet_ep_z_cell_z(i,j,k)+sheet_ep_z_cell_z(i-1,j,k)+sheet_ep_z_cell_z(i,j-1,k)&
                    +sheet_ep_z_cell_z(i-1,j-1,k))/4.0
                    sheet_sig_x_z(i,j,k)=(sheet_sig_x_cell_z(i,j,k)+sheet_sig_x_cell_z(i,j-1,k)+sheet_sig_x_cell_z(i,j,k)+&
                    sheet_sig_x_cell_z(i,j-1,k))/4.0
                    sheet_sig_y_z(i,j,k)=(sheet_sig_y_cell_z(i,j,k)+sheet_sig_y_cell_z(i-1,j,k)+sheet_sig_y_cell_z(i,j,k)+&
                    sheet_sig_y_cell_z(i-1,j,k))/4.0
                    sheet_sig_z_z(i,j,k)=(sheet_sig_z_cell_z(i,j,k)+sheet_sig_z_cell_z(i-1,j,k)+sheet_sig_z_cell_z(i,j-1,k)+&
                    sheet_sig_z_cell_z(i-1,j-1,k))/4.0
                end if
            end do
        end do
    end do

    if (is_plasma==1) then
        do k=2,z_size-1
            do j=2,y_size-1
                do i=2,x_size-1
                    do plasma_counter=1, num_poles
                        plasma_freq_x(plasma_counter,i,j,k)=(plasma_freq_cell_x(plasma_counter,i,j,k)+plasma_freq_cell_x(plasma_counter,i,j-1,k)+&
                        plasma_freq_cell_x(plasma_counter,i,j,k-1)+plasma_freq_cell_x(plasma_counter,i,j-1,k-1))/4.0
                        plasma_loss_x(plasma_counter,i,j,k)=(plasma_loss_cell_x(plasma_counter,i,j,k)+plasma_loss_cell_x(plasma_counter,i,j-1,k)+&
                        plasma_loss_cell_x(plasma_counter,i,j,k-1)+plasma_loss_cell_x(plasma_counter,i,j-1,k-1))/4.0

                        plasma_freq_y(plasma_counter,i,j,k)=(plasma_freq_cell_y(plasma_counter,i,j,k) + plasma_freq_cell_y(plasma_counter,i-1,j,k) + &
                        plasma_freq_cell_y(plasma_counter,i,j,k-1) + plasma_freq_cell_y(plasma_counter,i-1,j,k-1))/4.0
                        plasma_loss_y(plasma_counter,i,j,k)=(plasma_loss_cell_y(plasma_counter,i,j,k) + plasma_loss_cell_y(plasma_counter,i-1,j,k) + &
                        plasma_loss_cell_y(plasma_counter,i,j,k-1) + plasma_loss_cell_y(plasma_counter,i-1,j,k-1))/4.0

                        plasma_freq_z(plasma_counter,i,j,k)=(plasma_freq_cell_z(plasma_counter,i,j,k) + plasma_freq_cell_z(plasma_counter,i-1,j,k) + &
                        plasma_freq_cell_z(plasma_counter,i,j-1,k) + plasma_freq_cell_z(plasma_counter,i-1,j-1,k))/4.0
                        plasma_loss_z(plasma_counter,i,j,k)=(plasma_loss_cell_z(plasma_counter,i,j,k) + plasma_loss_cell_z(plasma_counter,i-1,j,k) + &
                        plasma_loss_cell_z(plasma_counter,i,j-1,k) + plasma_loss_cell_z(plasma_counter,i-1,j-1,k))/4.0
                    end do
                end do
            end do
        end do
    end if

    !we can't fill in yee grid point 1. This is fine if PML because defaults are vacuum. If pbc, we need to pull from the end.
    if (pbc_x==1) then
        do k=1,z_size-1
            do j=1,y_size-1
                relative_ep_x(1,j,k)=relative_ep_x(x_size-1,j,k)
                relative_ep_y(1,j,k)=relative_ep_y(x_size-1,j,k)
                relative_ep_z(1,j,k)=relative_ep_z(x_size-1,j,k)
                sigma_x(1,j,k)=sigma_x(x_size-1,j,k)
                sigma_y(1,j,k)=sigma_y(x_size-1,j,k)
                sigma_z(1,j,k)=sigma_z(x_size-1,j,k)
                !sheets x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_x(1,j,k)=sheet_ep_x_x(x_size-1,j,k)
                    sheet_ep_y_x(1,j,k)=sheet_ep_y_x(x_size-1,j,k)
                    sheet_ep_z_x(1,j,k)=sheet_ep_z_x(x_size-1,j,k)
                    sheet_sig_x_x(1,j,k)=sheet_sig_x_x(x_size-1,j,k)
                    sheet_sig_y_x(1,j,k)=sheet_sig_y_x(x_size-1,j,k)
                    sheet_sig_z_x(1,j,k)=sheet_sig_z_x(x_size-1,j,k)
                end if
                !sheets y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_y(1,j,k)=sheet_ep_x_y(x_size-1,j,k)
                    sheet_ep_y_y(1,j,k)=sheet_ep_y_y(x_size-1,j,k)
                    sheet_ep_z_y(1,j,k)=sheet_ep_z_y(x_size-1,j,k)
                    sheet_sig_x_y(1,j,k)=sheet_sig_x_y(x_size-1,j,k)
                    sheet_sig_y_y(1,j,k)=sheet_sig_y_y(x_size-1,j,k)
                    sheet_sig_z_y(1,j,k)=sheet_sig_z_y(x_size-1,j,k)
                end if
                !sheets z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_z(1,j,k)=sheet_ep_x_z(x_size-1,j,k)
                    sheet_ep_y_z(1,j,k)=sheet_ep_y_z(x_size-1,j,k)
                    sheet_ep_z_z(1,j,k)=sheet_ep_z_z(x_size-1,j,k)
                    sheet_sig_x_z(1,j,k)=sheet_sig_x_z(x_size-1,j,k)
                    sheet_sig_y_z(1,j,k)=sheet_sig_y_z(x_size-1,j,k)
                    sheet_sig_z_z(1,j,k)=sheet_sig_z_z(x_size-1,j,k)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_x(plasma_counter,1,j,k)=plasma_freq_x(plasma_counter,x_size-1,j,k)
                        plasma_loss_x(plasma_counter,1,j,k)=plasma_loss_x(plasma_counter,x_size-1,j,k)
                        plasma_freq_y(plasma_counter,1,j,k)=plasma_freq_y(plasma_counter,x_size-1,j,k)
                        plasma_loss_y(plasma_counter,1,j,k)=plasma_loss_y(plasma_counter,x_size-1,j,k)
                        plasma_freq_z(plasma_counter,1,j,k)=plasma_freq_z(plasma_counter,x_size-1,j,k)
                        plasma_loss_z(plasma_counter,1,j,k)=plasma_loss_z(plasma_counter,x_size-1,j,k)
                    end do
                end if
            end do
        end do
    end if
    if (pbc_z==1) then
        do j=1,y_size-1
            do i=1,x_size-1
                !blocks
                relative_ep_x(i,j,1)=relative_ep_x(i,j,z_size-1)
                relative_ep_y(i,j,1)=relative_ep_y(i,j,z_size-1)
                relative_ep_z(i,j,1)=relative_ep_z(i,j,z_size-1)
                sigma_x(i,j,1)=sigma_x(i,j,z_size-1)
                sigma_y(i,j,1)=sigma_y(i,j,z_size-1)
                sigma_z(i,j,1)=sigma_z(i,j,z_size-1)
                !sheets normal x
                if (num_sheets_x>0) then
                    sheet_ep_x_x(i,j,1)=sheet_ep_x_x(i,j,z_size-1)
                    sheet_ep_y_x(i,j,1)=sheet_ep_y_x(i,j,z_size-1)
                    sheet_ep_z_x(i,j,1)=sheet_ep_z_x(i,j,z_size-1)
                    sheet_sig_x_x(i,j,1)=sheet_sig_x_x(i,j,z_size-1)
                    sheet_sig_y_x(i,j,1)=sheet_sig_y_x(i,j,z_size-1)
                    sheet_sig_z_x(i,j,1)=sheet_sig_z_x(i,j,z_size-1)
                end if
                !sheets normal y
                if (num_sheets_y>0) then
                    sheet_ep_x_y(i,j,1)=sheet_ep_x_y(i,j,z_size-1)
                    sheet_ep_y_y(i,j,1)=sheet_ep_y_y(i,j,z_size-1)
                    sheet_ep_z_y(i,j,1)=sheet_ep_z_y(i,j,z_size-1)
                    sheet_sig_x_y(i,j,1)=sheet_sig_x_y(i,j,z_size-1)
                    sheet_sig_y_y(i,j,1)=sheet_sig_y_y(i,j,z_size-1)
                    sheet_sig_z_y(i,j,1)=sheet_sig_z_y(i,j,z_size-1)
                end if
                !sheets normal z
                if (num_sheets_z>0) then
                    sheet_ep_x_z(i,j,1)=sheet_ep_x_z(i,j,z_size-1)
                    sheet_ep_y_z(i,j,1)=sheet_ep_y_z(i,j,z_size-1)
                    sheet_ep_z_z(i,j,1)=sheet_ep_z_z(i,j,z_size-1)
                    sheet_sig_x_z(i,j,1)=sheet_sig_x_z(i,j,z_size-1)
                    sheet_sig_y_z(i,j,1)=sheet_sig_y_z(i,j,z_size-1)
                    sheet_sig_z_z(i,j,1)=sheet_sig_z_z(i,j,z_size-1)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_x(plasma_counter,i,j,1)=plasma_freq_x(plasma_counter,i,j,z_size-1)
                        plasma_loss_x(plasma_counter,i,j,1)=plasma_loss_x(plasma_counter,i,j,z_size-1)
                        plasma_freq_y(plasma_counter,i,j,1)=plasma_freq_y(plasma_counter,i,j,z_size-1)
                        plasma_loss_y(plasma_counter,i,j,1)=plasma_loss_y(plasma_counter,i,j,z_size-1)
                        plasma_freq_z(plasma_counter,i,j,1)=plasma_freq_z(plasma_counter,i,j,z_size-1)
                        plasma_loss_z(plasma_counter,i,j,1)=plasma_loss_z(plasma_counter,i,j,z_size-1)
                    end do
                end if
            end do
        end do
    end if
    if (pbc_y==1) then
        do k=1,z_size-1
            do i=1,x_size-1
                relative_ep_x(i,1,k)=relative_ep_x(i,y_size-1,k)
                relative_ep_y(i,1,k)=relative_ep_y(i,y_size-1,k)
                relative_ep_z(i,1,k)=relative_ep_z(i,y_size-1,k)
                sigma_x(i,1,k)=sigma_x(i,y_size-1,k)
                sigma_y(i,1,k)=sigma_y(i,y_size-1,k)
                sigma_z(i,1,k)=sigma_z(i,y_size-1,k)
                !sheets x normal
                if (num_sheets_x>0) then
                    sheet_ep_x_x(i,1,k)=sheet_ep_x_x(i,y_size-1,k)
                    sheet_ep_y_x(i,1,k)=sheet_ep_y_x(i,y_size-1,k)
                    sheet_ep_z_x(i,1,k)=sheet_ep_z_x(i,y_size-1,k)
                    sheet_sig_x_x(i,1,k)=sheet_sig_x_x(i,y_size-1,k)
                    sheet_sig_y_x(i,1,k)=sheet_sig_y_x(i,y_size-1,k)
                    sheet_sig_z_x(i,1,k)=sheet_sig_z_x(i,y_size-1,k)
                end if
                !sheets y normal
                if (num_sheets_y>0) then
                    sheet_ep_x_y(i,1,k)=sheet_ep_x_y(i,y_size-1,k)
                    sheet_ep_y_y(i,1,k)=sheet_ep_y_y(i,y_size-1,k)
                    sheet_ep_z_y(i,1,k)=sheet_ep_z_y(i,y_size-1,k)
                    sheet_sig_x_y(i,1,k)=sheet_sig_x_y(i,y_size-1,k)
                    sheet_sig_y_y(i,1,k)=sheet_sig_y_y(i,y_size-1,k)
                    sheet_sig_z_y(i,1,k)=sheet_sig_z_y(i,y_size-1,k)
                end if
                !sheets z normal
                if (num_sheets_z>0) then
                    sheet_ep_x_z(i,1,k)=sheet_ep_x_z(i,y_size-1,k)
                    sheet_ep_y_z(i,1,k)=sheet_ep_y_z(i,y_size-1,k)
                    sheet_ep_z_z(i,1,k)=sheet_ep_z_z(i,y_size-1,k)
                    sheet_sig_x_z(i,1,k)=sheet_sig_x_z(i,y_size-1,k)
                    sheet_sig_y_z(i,1,k)=sheet_sig_y_z(i,y_size-1,k)
                    sheet_sig_z_z(i,1,k)=sheet_sig_z_z(i,y_size-1,k)
                end if
                if (is_plasma==1) then
                    do plasma_counter=1, num_poles
                        plasma_freq_x(plasma_counter,i,1,k)=plasma_freq_x(plasma_counter,i,y_size-1,k)
                        plasma_loss_x(plasma_counter,i,1,k)=plasma_loss_x(plasma_counter,i,y_size-1,k)
                        plasma_freq_y(plasma_counter,i,1,k)=plasma_freq_y(plasma_counter,i,y_size-1,k)
                        plasma_loss_y(plasma_counter,i,1,k)=plasma_loss_y(plasma_counter,i,y_size-1,k)
                        plasma_freq_z(plasma_counter,i,1,k)=plasma_freq_z(plasma_counter,i,y_size-1,k)
                        plasma_loss_z(plasma_counter,i,1,k)=plasma_loss_z(plasma_counter,i,y_size-1,k)
                    end do
                end if
            end do
        end do
    end if

    !at the end of each contained section, deallocate memory we don't need anymore before allocating more
    ! Standard cell arrays
    deallocate(relative_ep_x_cell, sigma_x_cell, &
            relative_ep_y_cell, sigma_y_cell, &
            relative_ep_z_cell, sigma_z_cell)

    ! X-directed sheet material cells
    deallocate(sheet_ep_x_cell_x, sheet_ep_y_cell_x, sheet_ep_z_cell_x, &
            sheet_sig_x_cell_x, sheet_sig_y_cell_x, sheet_sig_z_cell_x)

    ! Y-directed sheet material cells
    deallocate(sheet_ep_x_cell_y, sheet_ep_y_cell_y, sheet_ep_z_cell_y, &
            sheet_sig_x_cell_y, sheet_sig_y_cell_y, sheet_sig_z_cell_y)

    ! Z-directed sheet material cells
    deallocate(sheet_ep_x_cell_z, sheet_ep_y_cell_z, sheet_ep_z_cell_z, &
            sheet_sig_x_cell_z, sheet_sig_y_cell_z, sheet_sig_z_cell_z)

    ! Plasma cell arrays
    deallocate(plasma_freq_cell_x, plasma_loss_cell_x, &
            plasma_freq_cell_y, plasma_loss_cell_y, &
            plasma_freq_cell_z, plasma_loss_cell_z)

    ! optional geoemtry and output geometry arrays are no longer needed either
    deallocate(output_geometry)
    if (trim(use_optional_bulk_geom_file)=='yes') then
        deallocate(optional_bulk_geom)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Creating Aux Eqs from geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !permability set to vacuum for all x,y,z
    DA = 1.0
    DB = (del_t/(mu_0)) 

    allocate(Gax(x_size-1,y_size-1,z_size-1))
    allocate(Gbx(x_size-1,y_size-1,z_size-1))
    allocate(Gay(x_size-1,y_size-1,z_size-1))
    allocate(Gby(x_size-1,y_size-1,z_size-1))
    allocate(Gaz(x_size-1,y_size-1,z_size-1))
    allocate(Gbz(x_size-1,y_size-1,z_size-1))

    !added specifically for plasmas
    if (is_plasma==1) then
        ! X-component allocations
        allocate(J_plasma_ax(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(J_plasma_bx(num_poles,x_size-1,y_size-1,z_size-1))
        ! Y-component allocations
        allocate(J_plasma_ay(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(J_plasma_by(num_poles,x_size-1,y_size-1,z_size-1))
        ! Z-component allocations
        allocate(J_plasma_az(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(J_plasma_bz(num_poles,x_size-1,y_size-1,z_size-1))
    else
        allocate(J_plasma_ax(1,1,1,1), J_plasma_bx(1,1,1,1))
        allocate(J_plasma_ay(1,1,1,1), J_plasma_by(1,1,1,1))
        allocate(J_plasma_az(1,1,1,1), J_plasma_bz(1,1,1,1))
        J_plasma_ax(1,1,1,1)=1.0
        J_plasma_bx(1,1,1,1)=0.0
        J_plasma_ay(1,1,1,1)=1.0
        J_plasma_by(1,1,1,1)=0.0
        J_plasma_az(1,1,1,1)=1.0
        J_plasma_bz(1,1,1,1)=0.0
    end if

    !G and D and Ja(b) don't need initilization because they are entirely filled before main fdtd runs
    !unless J isn't used and is dummied to 1, then we need to do something about it to avoid access issue
    !the access issue is for the speed up bounds section below in setup

    !now the aux equations for fdtd updates for bulk materials, will get overwritten if sheets present.
    if (num_poles==0) then
        do k = 1,z_size-1
            do j = 1,y_size-1
                do i = 1,x_size-1

                    gax(i,j,k) = (1.0 - sigma_x(i,j,k)*del_t / (2.0*ep_0*relative_ep_x(i,j,k))) / &
                    (1.0 + sigma_x(i,j,k) * del_t / (2.0*ep_0*relative_ep_x(i,j,k)))
                    
                    gbx(i,j,k) = (del_t/(ep_0*relative_ep_x(i,j,k))) / &
                    (1.0 + sigma_x(i,j,k)*del_t / (2.0*ep_0*relative_ep_x(i,j,k)))

                    gay(i,j,k) = (1.0 - sigma_y(i,j,k)*del_t / (2.0*ep_0*relative_ep_y(i,j,k))) / &
                    (1.0 + sigma_y(i,j,k) * del_t / (2.0*ep_0*relative_ep_y(i,j,k)))

                    gby(i,j,k) = (del_t/(ep_0*relative_ep_y(i,j,k))) / &
                    (1.0 + sigma_y(i,j,k)*del_t / (2.0*ep_0*relative_ep_y(i,j,k)))

                    gaz(i,j,k) = (1.0 - sigma_z(i,j,k)*del_t / (2.0*ep_0*relative_ep_z(i,j,k))) / &
                    (1.0 + sigma_z(i,j,k) * del_t / (2.0*ep_0*relative_ep_z(i,j,k)))
                    
                    gbz(i,j,k) = (del_t/(ep_0*relative_ep_z(i,j,k))) / &
                    (1.0 + sigma_z(i,j,k)*del_t / (2.0*ep_0*relative_ep_z(i,j,k)))

                end do   
            end do
        end do
    end if
    if (is_plasma==1) then
        do k = 1,z_size-1
            do j = 1,y_size-1
                do i = 1,x_size-1

                    temp_plasma_x=0.0
                    temp_plasma_y=0.0
                    temp_plasma_z=0.0

                    do plasma_counter=1, num_poles

                        J_plasma_ax(plasma_counter,i,j,k)=(1.0-plasma_loss_x(plasma_counter,i,j,k)*del_t/2.0)/&
                        (1.0+plasma_loss_x(plasma_counter,i,j,k)*del_t/2.0)

                        J_plasma_bx(plasma_counter,i,j,k) = &
                            REAL( (DBLE(plasma_freq_x(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                            (1.0 + plasma_loss_x(plasma_counter,i,j,k) * del_t / 2.0), 4)

                        J_plasma_ay(plasma_counter,i,j,k)=(1.0-plasma_loss_y(plasma_counter,i,j,k)*del_t/2.0)/&
                        (1.0+plasma_loss_y(plasma_counter,i,j,k)*del_t/2.0)

                        J_plasma_by(plasma_counter,i,j,k) = &
                            REAL( (DBLE(plasma_freq_y(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                            (1.0 + plasma_loss_y(plasma_counter,i,j,k) * del_t / 2.0), 4)

                        J_plasma_az(plasma_counter,i,j,k)=(1.0-plasma_loss_z(plasma_counter,i,j,k)*del_t/2.0)/&
                        (1.0+plasma_loss_z(plasma_counter,i,j,k)*del_t/2.0)

                        J_plasma_bz(plasma_counter,i,j,k) = &
                            REAL( (DBLE(plasma_freq_z(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                            (1.0 + plasma_loss_z(plasma_counter,i,j,k) * del_t / 2.0), 4)

                        temp_plasma_x=temp_plasma_x+J_plasma_bx(plasma_counter,i,j,k)
                        temp_plasma_y=temp_plasma_y+J_plasma_by(plasma_counter,i,j,k)
                        temp_plasma_z=temp_plasma_z+J_plasma_bz(plasma_counter,i,j,k)

                    end do

                    gax(i,j,k) = (1.0 - (sigma_x(i,j,k)+temp_plasma_x)*del_t / (2.0*ep_0*relative_ep_x(i,j,k))) / &
                    (1.0 + (sigma_x(i,j,k)+temp_plasma_x) * del_t / (2.0*ep_0*relative_ep_x(i,j,k)))
                    
                    gbx(i,j,k) = (del_t/(ep_0*relative_ep_x(i,j,k))) / &
                    (1.0 + (sigma_x(i,j,k)+temp_plasma_x)*del_t / (2.0*ep_0*relative_ep_x(i,j,k)))

                    gay(i,j,k) = (1.0 - (sigma_y(i,j,k)+temp_plasma_y)*del_t / (2.0*ep_0*relative_ep_y(i,j,k))) / &
                    (1.0 + (sigma_y(i,j,k)+temp_plasma_y) * del_t / (2.0*ep_0*relative_ep_y(i,j,k)))

                    gby(i,j,k) = (del_t/(ep_0*relative_ep_y(i,j,k))) / &
                    (1.0 + (sigma_y(i,j,k)+temp_plasma_y)*del_t / (2.0*ep_0*relative_ep_y(i,j,k)))

                    gaz(i,j,k) = (1.0 - (sigma_z(i,j,k)+temp_plasma_z)*del_t / (2.0*ep_0*relative_ep_z(i,j,k))) / &
                    (1.0 + (sigma_z(i,j,k)+temp_plasma_z) * del_t / (2.0*ep_0*relative_ep_z(i,j,k)))
                    
                    gbz(i,j,k) = (del_t/(ep_0*relative_ep_z(i,j,k))) / &
                    (1.0 + (sigma_z(i,j,k)+temp_plasma_z)*del_t / (2.0*ep_0*relative_ep_z(i,j,k)))

                end do   
            end do
        end do
    end if

    !I will loop though all x,y,z and see if a sheet conductivity (and epsilon in case we modify code later) has been udpated by user from vacuum values.
    !If so, GA,GB will be overwritten with an appropriate average.
    if (num_poles==0) then
        do k=1, z_size-1
            do j=1, y_size-1
                do i=1, x_size-1

                    !default all each time before begining
                    !Not necessary at every i,j,k, only when sheets are present, but this is the easiest way to do it since we don't know ahead of time if x,y, or z normal.
                    !Note a big time constraint.
                    sheet_ep_average_x=relative_ep_x(i,j,k)
                    sheet_sig_average_x=sigma_x(i,j,k)
                    sheet_ep_average_y=relative_ep_y(i,j,k)
                    sheet_sig_average_y=sigma_y(i,j,k)
                    sheet_ep_average_z=relative_ep_z(i,j,k)
                    sheet_sig_average_z=sigma_z(i,j,k)
                    counter_x_sheet=0
                    counter_y_sheet=0
                    counter_z_sheet=0

                    !First check for x normal sheets and add to average appropriately
                    if (num_sheets_x>0) then
                        if ((sheet_sig_y_x(i,j,k)>0) .or. (sheet_ep_y_x(i,j,k)>1)) then
                            sheet_ep_average_y=sheet_ep_average_y+(sheet_thickness/del_x)*(sheet_ep_y_x(i,j,k)-relative_ep_y(i,j,k))
                            sheet_sig_average_y=sheet_sig_average_y+(sheet_thickness/del_x)*(sheet_sig_y_x(i,j,k)-sigma_y(i,j,k))
                            counter_y_sheet=1
                        end if
                        if ((sheet_sig_z_x(i,j,k)>0) .or. (sheet_ep_z_x(i,j,k)>1)) then
                            sheet_ep_average_z=sheet_ep_average_z+(sheet_thickness/del_x)*(sheet_ep_z_x(i,j,k)-relative_ep_z(i,j,k))
                            sheet_sig_average_z=sheet_sig_average_z+(sheet_thickness/del_x)*(sheet_sig_z_x(i,j,k)-sigma_z(i,j,k))
                            counter_z_sheet=1
                        end if
                    end if

                    !Second check for y normal sheets and add to average appropriately
                    if (num_sheets_y>0) then
                        if ((sheet_sig_x_y(i,j,k)>0) .or. (sheet_ep_x_y(i,j,k)>1)) then
                            sheet_ep_average_x=sheet_ep_average_x+(sheet_thickness/del_y)*(sheet_ep_x_y(i,j,k)-relative_ep_x(i,j,k))
                            sheet_sig_average_x=sheet_sig_average_x+(sheet_thickness/del_y)*(sheet_sig_x_y(i,j,k)-sigma_x(i,j,k))
                            counter_x_sheet=1
                        end if
                        if ((sheet_sig_z_y(i,j,k)>0) .or. (sheet_ep_z_y(i,j,k)>1)) then
                            sheet_ep_average_z=sheet_ep_average_z+(sheet_thickness/del_y)*(sheet_ep_z_y(i,j,k)-relative_ep_z(i,j,k))
                            sheet_sig_average_z=sheet_sig_average_z+(sheet_thickness/del_y)*(sheet_sig_z_y(i,j,k)-sigma_z(i,j,k))
                            counter_z_sheet=1
                        end if
                    end if

                    !Lastly check for z normal sheets and add to average appropriately
                    if (num_sheets_z>0) then
                        if ((sheet_sig_x_z(i,j,k)>0) .or. (sheet_ep_x_z(i,j,k)>1)) then
                            sheet_ep_average_x=sheet_ep_average_x+(sheet_thickness/del_z)*(sheet_ep_x_z(i,j,k)-relative_ep_x(i,j,k))
                            sheet_sig_average_x=sheet_sig_average_x+(sheet_thickness/del_z)*(sheet_sig_x_z(i,j,k)-sigma_x(i,j,k))
                            counter_x_sheet=1
                        end if
                        if ((sheet_sig_y_z(i,j,k)>0) .or. (sheet_ep_y_z(i,j,k)>1)) then
                            sheet_ep_average_y=sheet_ep_average_y+(sheet_thickness/del_z)*(sheet_ep_y_z(i,j,k)-relative_ep_y(i,j,k))
                            sheet_sig_average_y=sheet_sig_average_y+(sheet_thickness/del_z)*(sheet_sig_y_z(i,j,k)-sigma_y(i,j,k))
                            counter_y_sheet=1
                        end if
                    end if
                    
                    !now for every i,j,k we know if a sheet is present and what the thickness is and have partially updated the averages
                    !now we need to finish the averages and update GA,GB
                    !if an average was updated above in any capacity, it moved the counter from zero to 1.
                    if ((counter_z_sheet)>0) then
                        !averages for z
                        gaz(i,j,k) = (1.0 - sheet_sig_average_z*del_t / (2.0*ep_0*sheet_ep_average_z)) / &
                        (1.0 + sheet_sig_average_z * del_t / (2.0*ep_0*sheet_ep_average_z))
                        gbz(i,j,k) = (del_t/(ep_0*sheet_ep_average_z)) / &
                        (1.0 + sheet_sig_average_z*del_t / (2.0*ep_0*sheet_ep_average_z))
                    end if
                    if ((counter_y_sheet)>0) then
                        !averages for y
                        gay(i,j,k) = (1.0 - sheet_sig_average_y*del_t / (2.0*ep_0*sheet_ep_average_y)) / &
                        (1.0 + sheet_sig_average_y * del_t / (2.0*ep_0*sheet_ep_average_y))
                        gby(i,j,k) = (del_t/(ep_0*sheet_ep_average_y)) / &
                        (1.0 + sheet_sig_average_y*del_t / (2.0*ep_0*sheet_ep_average_y))
                    end if
                    if ((counter_x_sheet)>0) then 
                        !averges for x
                        gax(i,j,k) = (1.0 - sheet_sig_average_x*del_t / (2.0*ep_0*sheet_ep_average_x)) / &
                        (1.0 + sheet_sig_average_x * del_t / (2.0*ep_0*sheet_ep_average_x))
                        gbx(i,j,k) = (del_t/(ep_0*sheet_ep_average_x)) / &
                        (1.0 + sheet_sig_average_x*del_t / (2.0*ep_0*sheet_ep_average_x))
                    end if

                end do
            end do
        end do
    end if
    if (is_plasma==1) then
        do k=1, z_size-1
            do j=1, y_size-1
                do i=1, x_size-1

                    !default all each time before begining
                    !Not necessary at every i,j,k, only when sheets are present, but this is the easiest way to do it since we don't know ahead of time if x,y, or z normal.
                    !Note a big time constraint.
                    sheet_ep_average_x=relative_ep_x(i,j,k)
                    sheet_sig_average_x=sigma_x(i,j,k)
                    sheet_ep_average_y=relative_ep_y(i,j,k)
                    sheet_sig_average_y=sigma_y(i,j,k)
                    sheet_ep_average_z=relative_ep_z(i,j,k)
                    sheet_sig_average_z=sigma_z(i,j,k)
                    counter_x_sheet=0
                    counter_y_sheet=0
                    counter_z_sheet=0

                    !First check for x normal sheets and add to average appropriately
                    if (num_sheets_x>0) then
                        if ((sheet_sig_y_x(i,j,k)>0) .or. (sheet_ep_y_x(i,j,k)>1)) then
                            sheet_ep_average_y=sheet_ep_average_y+(sheet_thickness/del_x)*(sheet_ep_y_x(i,j,k)-relative_ep_y(i,j,k))
                            sheet_sig_average_y=sheet_sig_average_y+(sheet_thickness/del_x)*(sheet_sig_y_x(i,j,k)-sigma_y(i,j,k))
                            counter_y_sheet=1
                        end if
                        if ((sheet_sig_z_x(i,j,k)>0) .or. (sheet_ep_z_x(i,j,k)>1)) then
                            sheet_ep_average_z=sheet_ep_average_z+(sheet_thickness/del_x)*(sheet_ep_z_x(i,j,k)-relative_ep_z(i,j,k))
                            sheet_sig_average_z=sheet_sig_average_z+(sheet_thickness/del_x)*(sheet_sig_z_x(i,j,k)-sigma_z(i,j,k))
                            counter_z_sheet=1
                        end if
                    end if

                    !Second check for y normal sheets and add to average appropriately
                    if (num_sheets_y>0) then
                        if ((sheet_sig_x_y(i,j,k)>0) .or. (sheet_ep_x_y(i,j,k)>1)) then
                            sheet_ep_average_x=sheet_ep_average_x+(sheet_thickness/del_y)*(sheet_ep_x_y(i,j,k)-relative_ep_x(i,j,k))
                            sheet_sig_average_x=sheet_sig_average_x+(sheet_thickness/del_y)*(sheet_sig_x_y(i,j,k)-sigma_x(i,j,k))
                            counter_x_sheet=1
                        end if
                        if ((sheet_sig_z_y(i,j,k)>0) .or. (sheet_ep_z_y(i,j,k)>1)) then
                            sheet_ep_average_z=sheet_ep_average_z+(sheet_thickness/del_y)*(sheet_ep_z_y(i,j,k)-relative_ep_z(i,j,k))
                            sheet_sig_average_z=sheet_sig_average_z+(sheet_thickness/del_y)*(sheet_sig_z_y(i,j,k)-sigma_z(i,j,k))
                            counter_z_sheet=1
                        end if
                    end if

                    !Lastly check for z normal sheets and add to average appropriately
                    if (num_sheets_z>0) then
                        if ((sheet_sig_x_z(i,j,k)>0) .or. (sheet_ep_x_z(i,j,k)>1)) then
                            sheet_ep_average_x=sheet_ep_average_x+(sheet_thickness/del_z)*(sheet_ep_x_z(i,j,k)-relative_ep_x(i,j,k))
                            sheet_sig_average_x=sheet_sig_average_x+(sheet_thickness/del_z)*(sheet_sig_x_z(i,j,k)-sigma_x(i,j,k))
                            counter_x_sheet=1
                        end if
                        if ((sheet_sig_y_z(i,j,k)>0) .or. (sheet_ep_y_z(i,j,k)>1)) then
                            sheet_ep_average_y=sheet_ep_average_y+(sheet_thickness/del_z)*(sheet_ep_y_z(i,j,k)-relative_ep_y(i,j,k))
                            sheet_sig_average_y=sheet_sig_average_y+(sheet_thickness/del_z)*(sheet_sig_y_z(i,j,k)-sigma_y(i,j,k))
                            counter_y_sheet=1
                        end if
                    end if
                    
                    !now for every i,j,k we know if a sheet is present and what the thickness is and have partially updated the averages
                    !now we need to finish the averages and update GA,GB
                    !if an average was updated above in any capacity, it moved the counter from zero to 1.
                    if ((counter_z_sheet)>0) then
                        temp_plasma_z=0.0
                        do plasma_counter=1, num_poles
                            J_plasma_bz(plasma_counter,i,j,k) = &
                                REAL( (DBLE(plasma_freq_z(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                                (1.0 + plasma_loss_z(plasma_counter,i,j,k) * del_t / 2.0), 4)
                            temp_plasma_z=temp_plasma_z+J_plasma_bz(plasma_counter,i,j,k)
                        end do
                        !averages for z
                        gaz(i,j,k) = (1.0 - (sheet_sig_average_z+temp_plasma_z)*del_t / (2.0*ep_0*sheet_ep_average_z)) / &
                        (1.0 + (sheet_sig_average_z+temp_plasma_z) * del_t / (2.0*ep_0*sheet_ep_average_z))
                        gbz(i,j,k) = (del_t/(ep_0*sheet_ep_average_z)) / &
                        (1.0 + (sheet_sig_average_z+temp_plasma_z)*del_t / (2.0*ep_0*sheet_ep_average_z))
                    end if
                    if ((counter_y_sheet)>0) then
                        temp_plasma_y=0.0
                        do plasma_counter=1, num_poles
                            J_plasma_by(plasma_counter,i,j,k) = &
                                REAL( (DBLE(plasma_freq_y(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                                (1.0 + plasma_loss_y(plasma_counter,i,j,k) * del_t / 2.0), 4)
                            temp_plasma_y=temp_plasma_y+J_plasma_by(plasma_counter,i,j,k)
                        end do
                        !averages for y
                        gay(i,j,k) = (1.0 - (sheet_sig_average_y+temp_plasma_y)*del_t / (2.0*ep_0*sheet_ep_average_y)) / &
                        (1.0 + (sheet_sig_average_y+temp_plasma_y) * del_t / (2.0*ep_0*sheet_ep_average_y))
                        gby(i,j,k) = (del_t/(ep_0*sheet_ep_average_y)) / &
                        (1.0 + (sheet_sig_average_y+temp_plasma_y)*del_t / (2.0*ep_0*sheet_ep_average_y))
                    end if
                    if ((counter_x_sheet)>0) then 
                        temp_plasma_x=0.0
                        do plasma_counter=1, num_poles
                            J_plasma_bx(plasma_counter,i,j,k) = &
                                REAL( (DBLE(plasma_freq_x(plasma_counter,i,j,k))**2 * ep_0 * del_t / 2.0D0) / &
                                (1.0 + plasma_loss_x(plasma_counter,i,j,k) * del_t / 2.0), 4)
                            temp_plasma_x=temp_plasma_x+J_plasma_bx(plasma_counter,i,j,k)
                        end do
                        !averges for x
                        gax(i,j,k) = (1.0 - (sheet_sig_average_x+temp_plasma_x)*del_t / (2.0*ep_0*sheet_ep_average_x)) / &
                        (1.0 + (sheet_sig_average_x+temp_plasma_x) * del_t / (2.0*ep_0*sheet_ep_average_x))
                        gbx(i,j,k) = (del_t/(ep_0*sheet_ep_average_x)) / &
                        (1.0 + (sheet_sig_average_x+temp_plasma_x)*del_t / (2.0*ep_0*sheet_ep_average_x))
                    end if

                end do
            end do
        end do
    end if

    !now do two things: 
    !1. check if plasma was actually even used
    !2. if so, get a bounding for each component of J,E so we don't search the whole space in main fdtd update
    !This is purely for speed purposes.
    if (ANY(J_plasma_ax>1) .or. ANY(J_plasma_bx>0) .or. ANY(J_plasma_ay>1) .or. ANY(J_plasma_by>0) .or. ANY(J_plasma_az>1) .or. ANY(J_plasma_bz>0)) then
        is_plasma=1

        !x fields
        plasma_max_xfields_zpos=2
        plasma_min_xfields_zpos=z_size-1
        plasma_max_xfields_ypos=2
        plasma_min_xfields_ypos=y_size-1
        plasma_max_xfields_xpos=1
        plasma_min_xfields_xpos=x_size-1
        !y fields
        plasma_max_yfields_zpos=2
        plasma_min_yfields_zpos=z_size-1
        plasma_max_yfields_ypos=1
        plasma_min_yfields_ypos=y_size-1
        plasma_max_yfields_xpos=2
        plasma_min_yfields_xpos=x_size-1
        !z fields
        plasma_max_zfields_zpos=1
        plasma_min_zfields_zpos=z_size-1
        plasma_max_zfields_ypos=2
        plasma_min_zfields_ypos=y_size-1
        plasma_max_zfields_xpos=2
        plasma_min_zfields_xpos=x_size-1

        !first x directed fields
        do k=2, z_size-1
            do j=2, y_size-1
                do i=1, x_size-1
                    do ii=1, num_poles
                        if ((J_plasma_ax(ii,i,j,k)>1) .or. (J_plasma_bx(ii,i,j,k)>0)) then
                            if (k<plasma_min_xfields_zpos) then
                                plasma_min_xfields_zpos=k
                            end if
                            if (k>plasma_max_xfields_zpos) then
                                plasma_max_xfields_zpos=k
                            end if
                            if (j<plasma_min_xfields_ypos) then
                                plasma_min_xfields_ypos=j
                            end if
                            if (j>plasma_max_xfields_ypos) then
                                plasma_max_xfields_ypos=j
                            end if
                            if (i<plasma_min_xfields_xpos) then
                                plasma_min_xfields_xpos=i
                            end if
                            if (i>plasma_max_xfields_xpos) then
                                plasma_max_xfields_xpos=i
                            end if
                            !any pole is sufficent to mark the index
                            exit
                        end if
                    end do
                end do
            end do
        end do

        !now y directed fields
        do k=2, z_size-1
            do j=1, y_size-1
                do i=2, x_size-1
                    do ii=1, num_poles
                        if ((J_plasma_ay(ii,i,j,k)>1) .or. (J_plasma_by(ii,i,j,k)>0)) then
                            if (k<plasma_min_yfields_zpos) then
                                plasma_min_yfields_zpos=k
                            end if
                            if (k>plasma_max_yfields_zpos) then
                                plasma_max_yfields_zpos=k
                            end if
                            if (j<plasma_min_yfields_ypos) then
                                plasma_min_yfields_ypos=j
                            end if
                            if (j>plasma_max_yfields_ypos) then
                                plasma_max_yfields_ypos=j
                            end if
                            if (i<plasma_min_yfields_xpos) then
                                plasma_min_yfields_xpos=i
                            end if
                            if (i>plasma_max_yfields_xpos) then
                                plasma_max_yfields_xpos=i
                            end if
                            !any pole is sufficient to mark the index
                            exit
                        end if
                    end do
                end do
            end do
        end do

        !now lastly z directed fields
        do k=1, z_size-1
            do j=2, y_size-1
                do i=2, x_size-1
                    do ii=1, num_poles
                        if ((J_plasma_az(ii,i,j,k)>1) .or. (J_plasma_bz(ii,i,j,k)>0)) then
                            if (k<plasma_min_zfields_zpos) then
                                plasma_min_zfields_zpos=k
                            end if
                            if (k>plasma_max_zfields_zpos) then
                                plasma_max_zfields_zpos=k
                            end if
                            if (j<plasma_min_zfields_ypos) then
                                plasma_min_zfields_ypos=j
                            end if
                            if (j>plasma_max_zfields_ypos) then
                                plasma_max_zfields_ypos=j
                            end if
                            if (i<plasma_min_zfields_xpos) then
                                plasma_min_zfields_xpos=i
                            end if
                            if (i>plasma_max_zfields_xpos) then
                                plasma_max_zfields_xpos=i
                            end if
                            !any pole is sufficient to mark the index
                            exit
                        end if
                    end do
                end do
            end do
        end do
    else
        is_plasma=0
    end if

    !if you want to check the bounds we ended up sweeping, do these lines
    !write(*,*) plasma_max_xfields_zpos,plasma_min_xfields_zpos,plasma_max_xfields_ypos, plasma_min_xfields_ypos,plasma_max_xfields_xpos,plasma_min_xfields_xpos
    !write(*,*) plasma_max_yfields_zpos,plasma_min_yfields_zpos,plasma_max_yfields_ypos, plasma_min_yfields_ypos,plasma_max_yfields_xpos,plasma_min_yfields_xpos
    !write(*,*) plasma_max_zfields_zpos,plasma_min_zfields_zpos,plasma_max_zfields_ypos, plasma_min_zfields_ypos,plasma_max_zfields_xpos,plasma_min_zfields_xpos

    !at the end of each contained section, deallocate memory we don't need anymore before allocating more
    ! Standard arrays - SPICE version is unique in that relative x,y,z don't get deallocated here
#ifdef use_spice_version
    deallocate(sigma_x, &
            sigma_y, &
            sigma_z)
#endif
#ifndef use_spice_version
    deallocate(relative_ep_x, sigma_x, &
            relative_ep_y, sigma_y, &
            relative_ep_z, sigma_z)
#endif

    ! X-directed sheet material - keeping xx
    deallocate(sheet_ep_y_x, sheet_ep_z_x, &
            sheet_sig_y_x, sheet_sig_z_x)

    ! Y-directed sheet material - keeping yy
    deallocate(sheet_ep_x_y, sheet_ep_z_y, &
            sheet_sig_x_y, sheet_sig_z_y)

    ! Z-directed sheet material - keeping zz
    deallocate(sheet_ep_x_z, sheet_ep_y_z, &
            sheet_sig_x_z, sheet_sig_y_z)

    ! Plasma arrays
    deallocate(plasma_freq_x, plasma_loss_x, &
            plasma_freq_y, plasma_loss_y, &
            plasma_freq_z, plasma_loss_z)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup CPML !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !PML arrays - not all are used in unit cell version 
    !but it's a small memory and time add-on - could refine later to reduce memory requirements a tiny, tiny bit

    allocate(psi_Ezx_1(nxPML_1,y_size-1,z_size-1))
    allocate(psi_Ezx_2(nxPML_2,y_size-1,z_size-1))
    allocate(psi_Hyx_1(nxPML_1-1,y_size-1,z_size-1))
    allocate(psi_Hyx_2(nxPML_2-1,y_size-1,z_size-1))
    allocate(psi_Ezy_1(x_size-1,nyPML_1,z_size-1))
    allocate(psi_Ezy_2(x_size-1,nyPML_2,z_size-1))
    allocate(psi_Hxy_1(x_size-1,nyPML_1-1,z_size-1))
    allocate(psi_Hxy_2(x_size-1,nyPML_2-1,z_size-1))
    allocate(psi_Hxz_1(x_size-1,y_size-1,nzPML_1-1))
    allocate(psi_Hxz_2(x_size-1,y_size-1,nzPML_2-1))
    allocate(psi_Hyz_1(x_size-1,y_size-1,nzPML_1-1))
    allocate(psi_Hyz_2(x_size-1,y_size-1,nzPML_2-1))
    allocate(psi_Exz_1(x_size-1,y_size-1,nzPML_1))
    allocate(psi_Exz_2(x_size-1,y_size-1,nzPML_2))
    allocate(psi_Eyz_1(x_size-1,y_size-1,nzPML_1))
    allocate(psi_Eyz_2(x_size-1,y_size-1,nzPML_2))
    allocate(psi_Hzx_1(nxPML_1-1,y_size-1,z_size))
    allocate(psi_Eyx_1(nxPML_1,y_size-1,z_size-1))
    allocate(psi_Hzx_2(nxPML_2-1,y_size-1,z_size-1))
    allocate(psi_Eyx_2(nxPML_2,y_size-1,z_size-1))
    allocate(psi_Hzy_1(x_size-1,nyPML_1-1,z_size-1))
    allocate(psi_Exy_1(x_size-1,nyPML_1,z_size-1))
    allocate(psi_Hzy_2(x_size-1,nyPML_2-1,z_size-1))
    allocate(psi_Exy_2(x_size-1,nyPML_2,z_size-1))
    allocate(be_x_1(nxPML_1))
    allocate(ce_x_1(nxPML_1))
    allocate(alphae_x_PML_1(nxPML_1))
    allocate(sige_x_PML_1(nxPML_1))
    allocate(kappae_x_PML_1(nxPML_1))
    allocate(bh_x_1(nxPML_1-1))
    allocate(ch_x_1(nxPML_1-1))
    allocate(alphah_x_PML_1(nxPML_1-1))
    allocate(sigh_x_PML_1(nxPML_1-1))
    allocate(kappah_x_PML_1(nxPML_1-1))
    allocate(be_x_2(nxPML_2))
    allocate(ce_x_2(nxPML_2))
    allocate(alphae_x_PML_2(nxPML_2))
    allocate(sige_x_PML_2(nxPML_2))
    allocate(kappae_x_PML_2(nxPML_2))
    allocate(bh_x_2(nxPML_2-1))
    allocate(ch_x_2(nxPML_2-1))
    allocate(alphah_x_PML_2(nxPML_2-1))
    allocate(sigh_x_PML_2(nxPML_2-1))
    allocate(kappah_x_PML_2(nxPML_2-1))
    allocate(be_y_1(nyPML_1))
    allocate(ce_y_1(nyPML_1))
    allocate(alphae_y_PML_1(nyPML_1))
    allocate(sige_y_PML_1(nyPML_1))
    allocate(kappae_y_PML_1(nyPML_1))
    allocate(bh_y_1(nyPML_1-1))
    allocate(ch_y_1(nyPML_1-1))
    allocate(alphah_y_PML_1(nyPML_1-1))
    allocate(sigh_y_PML_1(nyPML_1-1))
    allocate(kappah_y_PML_1(nyPML_1-1))
    allocate(be_y_2(nyPML_2))
    allocate(ce_y_2(nyPML_2))
    allocate(alphae_y_PML_2(nyPML_2))
    allocate(sige_y_PML_2(nyPML_2))
    allocate(kappae_y_PML_2(nyPML_2))
    allocate(bh_y_2(nyPML_2-1))
    allocate(ch_y_2(nyPML_2-1))
    allocate(alphah_y_PML_2(nyPML_2-1))
    allocate(sigh_y_PML_2(nyPML_2-1))
    allocate(kappah_y_PML_2(nyPML_2-1))
    allocate(be_z_1(nzPML_1))
    allocate(ce_z_1(nzPML_1))
    allocate(alphae_z_PML_1(nzPML_1))
    allocate(sige_z_PML_1(nzPML_1))
    allocate(kappae_z_PML_1(nzPML_1))
    allocate(bh_z_1(nzPML_1-1))
    allocate(ch_z_1(nzPML_1-1))
    allocate(alphah_z_PML_1(nzPML_1-1))
    allocate(sigh_z_PML_1(nzPML_1-1))
    allocate(kappah_z_PML_1(nzPML_1-1))
    allocate(be_z_2(nzPML_2))
    allocate(ce_z_2(nzPML_2))
    allocate(alphae_z_PML_2(nzPML_2))
    allocate(sige_z_PML_2(nzPML_2))
    allocate(kappae_z_PML_2(nzPML_2))
    allocate(bh_z_2(nzPML_2-1))
    allocate(ch_z_2(nzPML_2-1))
    allocate(alphah_z_PML_2(nzPML_2-1))
    allocate(sigh_z_PML_2(nzPML_2-1))
    allocate(kappah_z_PML_2(nzPML_2-1))

    !Now initalize PML field arrays, won't intialize other arrays as they are fully filled in below
    psi_Exy_1(:,:,:) = 0.0+0.0*imag
    psi_Exy_2(:,:,:) = 0.0+0.0*imag
    psi_Exz_1(:,:,:) = 0.0+0.0*imag
    psi_Exz_2(:,:,:) = 0.0+0.0*imag
    psi_Eyx_1(:,:,:) = 0.0+0.0*imag
    psi_Eyx_2(:,:,:) = 0.0+0.0*imag
    psi_Eyz_1(:,:,:) = 0.0+0.0*imag
    psi_Eyz_2(:,:,:) = 0.0+0.0*imag
    psi_Ezy_1(:,:,:) = 0.0+0.0*imag
    psi_Ezy_2(:,:,:) = 0.0+0.0*imag
    psi_Ezx_1(:,:,:) = 0.0+0.0*imag
    psi_Ezx_2(:,:,:) = 0.0+0.0*imag
    psi_Hxy_1(:,:,:) = 0.0+0.0*imag
    psi_Hxy_2(:,:,:) = 0.0+0.0*imag
    psi_Hxz_1(:,:,:) = 0.0+0.0*imag
    psi_Hxz_2(:,:,:) = 0.0+0.0*imag
    psi_Hyx_1(:,:,:) = 0.0+0.0*imag
    psi_Hyx_2(:,:,:) = 0.0+0.0*imag
    psi_Hyz_1(:,:,:) = 0.0+0.0*imag
    psi_Hyz_2(:,:,:) = 0.0+0.0*imag
    psi_Hzy_1(:,:,:) = 0.0+0.0*imag
    psi_Hzy_2(:,:,:) = 0.0+0.0*imag
    psi_Hzx_1(:,:,:) = 0.0+0.0*imag
    psi_Hzx_2(:,:,:) = 0.0+0.0*imag

    !Setup PML variables used in pml field arrays
    do i = 1,nxPML_1
        sige_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i) / (nxPML_1 - 1.0) )**m
        alphae_x_PML_1(i) = alpha_x_max*((i-1.0)/(nxPML_1-1.0))**ma
        kappae_x_PML_1(i) = 1.0+(kappa_x_max-1.0)* &
        ((nxPML_1 - i) / (nxPML_1 - 1.0))**m
        be_x_1(i) = EXP(-(sige_x_PML_1(i) / kappae_x_PML_1(i) + &
        alphae_x_PML_1(i))*del_t/ep_0)
        if ((sige_x_PML_1(i) == 0.0) .and. &
        (alphae_x_PML_1(i) == 0.0) .and. (i == nxPML_1)) then
            ce_x_1(i) = 0.0
        else
            ce_x_1(i) = sige_x_PML_1(i)*(be_x_1(i)-1.0)/ &
            (sige_x_PML_1(i)+kappae_x_PML_1(i)*alphae_x_PML_1(i)) &
            /kappae_x_PML_1(i)
        end if
    end do

    do i = 1,nxPML_1-1
        sigh_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i - 0.5)/(nxPML_1-1.0))**m
        alphah_x_PML_1(i) = alpha_x_max*((i-0.5)/(nxPML_1-1.0))**ma
        kappah_x_PML_1(i) = 1.0+(kappa_x_max-1.0)* &
        ((nxPML_1 - i - 0.5) / (nxPML_1 - 1.0))**m
        bh_x_1(i) = EXP(-(sigh_x_PML_1(i) / kappah_x_PML_1(i) + &
        alphah_x_PML_1(i))*del_t/ep_0)
        ch_x_1(i) = sigh_x_PML_1(i)*(bh_x_1(i)-1.0)/  &
        (sigh_x_PML_1(i)+kappah_x_PML_1(i)*alphah_x_PML_1(i)) &
        /kappah_x_PML_1(i)
    end do

    do i = 1,nxPML_2
        sige_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i) / (nxPML_2 - 1.0) )**m
        alphae_x_PML_2(i) = alpha_x_max*((i-1.0)/(nxPML_2-1.0))**ma
        kappae_x_PML_2(i) = 1.0+(kappa_x_max-1.0)* &
        ((nxPML_2 - i) / (nxPML_2 - 1.0))**m
        be_x_2(i) = EXP(-(sige_x_PML_2(i) / kappae_x_PML_2(i) +  &
        alphae_x_PML_2(i))*del_t/ep_0)
        if ((sige_x_PML_2(i) == 0.0) .and. &
        (alphae_x_PML_2(i) == 0.0) .and. (i == nxPML_2)) then
            ce_x_2(i) = 0.0
        else
            ce_x_2(i) = sige_x_PML_2(i)*(be_x_2(i)-1.0)/ &
            (sige_x_PML_2(i)+kappae_x_PML_2(i)*alphae_x_PML_2(i)) &
            /kappae_x_PML_2(i)
        end if
    end do

    do i = 1,nxPML_2-1
        sigh_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i - 0.5)/(nxPML_2-1.0))**m
        alphah_x_PML_2(i) = alpha_x_max*((i-0.5)/(nxPML_2-1.0))**ma
        kappah_x_PML_2(i) = 1.0+(kappa_x_max-1.0)* &
        ((nxPML_2 - i - 0.5) / (nxPML_2 - 1.0))**m
        bh_x_2(i) = EXP(-(sigh_x_PML_2(i) / kappah_x_PML_2(i) + &
        alphah_x_PML_2(i))*del_t/ep_0)
        ch_x_2(i) = sigh_x_PML_2(i)*(bh_x_2(i)-1.0)/ &
        (sigh_x_PML_2(i)+kappah_x_PML_2(i)*alphah_x_PML_2(i)) &
        /kappah_x_PML_2(i)
    end do

    do j = 1,nyPML_1
        sige_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j ) / (nyPML_1 - 1.0) )**m
        alphae_y_PML_1(j) = alpha_y_max*((j-1)/(nyPML_1-1.0))**ma
        kappae_y_PML_1(j) = 1.0+(kappa_y_max-1.0)* &
        ((nyPML_1 - j) / (nyPML_1 - 1.0))**m
        be_y_1(j) = EXP(-(sige_y_PML_1(j) / kappae_y_PML_1(j) +   &
        alphae_y_PML_1(j))*del_t/ep_0)
        if ((sige_y_PML_1(j) == 0.0) .and. &
        (alphae_y_PML_1(j) == 0.0) .and. (j == nyPML_1)) then
            ce_y_1(j) = 0.0
        else
            ce_y_1(j) = sige_y_PML_1(j)*(be_y_1(j)-1.0)/ &
            (sige_y_PML_1(j)+kappae_y_PML_1(j)*alphae_y_PML_1(j)) &
            /kappae_y_PML_1(j)
        end if
    end do

    do j = 1,nyPML_1-1
        sigh_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j - 0.5)/(nyPML_1-1.0))**m
        alphah_y_PML_1(j) = alpha_y_max*((j-0.5)/(nyPML_1-1.0))**ma
        kappah_y_PML_1(j) = 1.0+(kappa_y_max-1.0)* &
        ((nyPML_1 - j - 0.5) / (nyPML_1 - 1.0))**m
        bh_y_1(j) = EXP(-(sigh_y_PML_1(j) / kappah_y_PML_1(j) + &
        alphah_y_PML_1(j))*del_t/ep_0)
        ch_y_1(j) = sigh_y_PML_1(j)*(bh_y_1(j)-1.0)/ &
        (sigh_y_PML_1(j)+kappah_y_PML_1(j)*alphah_y_PML_1(j)) &
        /kappah_y_PML_1(j)
    end do

    do j = 1,nyPML_2
        sige_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j ) / (nyPML_2 - 1.0) )**m
        alphae_y_PML_2(j) = alpha_y_max*((j-1)/(nyPML_2-1.0))**ma
        kappae_y_PML_2(j) = 1.0+(kappa_y_max-1.0)* &
        ((nyPML_2 - j) / (nyPML_2 - 1.0))**m
        be_y_2(j) = EXP(-(sige_y_PML_2(j) / kappae_y_PML_2(j) + &
        alphae_y_PML_2(j))*del_t/ep_0)
        if ((sige_y_PML_2(j) == 0.0) .and. &
        (alphae_y_PML_2(j) == 0.0) .and. (j == nyPML_2)) then
            ce_y_2(j) = 0.0
        else
            ce_y_2(j) = sige_y_PML_2(j)*(be_y_2(j)-1.0)/ &
            (sige_y_PML_2(j)+kappae_y_PML_2(j)*alphae_y_PML_2(j)) &
            /kappae_y_PML_2(j)
        end if
    end do

    do j = 1,nyPML_2-1
        sigh_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j - 0.5)/(nyPML_2-1.0))**m
        alphah_y_PML_2(j) = alpha_y_max*((j-0.5)/(nyPML_2-1.0))**ma
        kappah_y_PML_2(j) = 1.0+(kappa_y_max-1.0)* &
        ((nyPML_2 - j - 0.5) / (nyPML_2 - 1.0))**m
        bh_y_2(j) = EXP(-(sigh_y_PML_2(j) / kappah_y_PML_2(j) + &
        alphah_y_PML_2(j))*del_t/ep_0)
        ch_y_2(j) = sigh_y_PML_2(j)*(bh_y_2(j)-1.0)/ &
        (sigh_y_PML_2(j)+kappah_y_PML_2(j)*alphah_y_PML_2(j)) &
        /kappah_y_PML_2(j)
    end do

    do k = 1,nzPML_1
        sige_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k ) / (nzPML_1 - 1.0) )**m
        alphae_z_PML_1(k) = alpha_z_max*((k-1)/(nzPML_1-1.0))**ma
        kappae_z_PML_1(k) = 1.0+(kappa_z_max-1.0)* &
        ((nzPML_1 - k) / (nzPML_1 - 1.0))**m
        be_z_1(k) = EXP(-(sige_z_PML_1(k) / kappae_z_PML_1(k) + &
        alphae_z_PML_1(k))*del_t/ep_0)
        if ((sige_z_PML_1(k) == 0.0) .and. &
        (alphae_z_PML_1(k) == 0.0) .and. (k == nzPML_1)) then
            ce_z_1(k) = 0.0
        else
            ce_z_1(k) = sige_z_PML_1(k)*(be_z_1(k)-1.0)/ &
            (sige_z_PML_1(k)+kappae_z_PML_1(k)*alphae_z_PML_1(k)) &
            /kappae_z_PML_1(k)
        end if
    end do

    do k = 1,nzPML_1-1
        sigh_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k - 0.5)/(nzPML_1-1.0))**m
        alphah_z_PML_1(k) = alpha_z_max*((k-0.5)/(nzPML_1-1.0))**ma
        kappah_z_PML_1(k) = 1.0+(kappa_z_max-1.0)* &
        ((nzPML_1 - k - 0.5) / (nzPML_1 - 1.0))**m
        bh_z_1(k) = EXP(-(sigh_z_PML_1(k) / kappah_z_PML_1(k) + &
        alphah_z_PML_1(k))*del_t/ep_0)
        ch_z_1(k) = sigh_z_PML_1(k)*(bh_z_1(k)-1.0)/ &
        (sigh_z_PML_1(k)+kappah_z_PML_1(k)*alphah_z_PML_1(k)) &
        /kappah_z_PML_1(k)
    end do

    do k = 1,nzPML_2
        sige_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k ) / (nzPML_2 - 1.0) )**m
        alphae_z_PML_2(k) = alpha_z_max*((k-1)/(nzPML_2-1.0))**ma
        kappae_z_PML_2(k) = 1.0+(kappa_z_max-1.0)* &
        ((nzPML_2 - k) / (nzPML_2 - 1.0))**m
        be_z_2(k) = EXP(-(sige_z_PML_2(k) / kappae_z_PML_2(k) + &
        alphae_z_PML_2(k))*del_t/ep_0)
        if ((sige_z_PML_2(k) == 0.0) .and. &
        (alphae_z_PML_2(k) == 0.0) .and. (k == nzPML_2)) then
            ce_z_2(k) = 0.0
        else
            ce_z_2(k) = sige_z_PML_2(k)*(be_z_2(k)-1.0)/ &
            (sige_z_PML_2(k)+kappae_z_PML_2(k)*alphae_z_PML_2(k)) &
            /kappae_z_PML_2(k)
        end if
    end do

    do k = 1,nzPML_2-1
        sigh_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k - 0.5)/(nzPML_2-1.0))**m
        alphah_z_PML_2(k) = alpha_z_max*((k-0.5)/(nzPML_2-1.0))**ma
        kappah_z_PML_2(k) = 1.0+(kappa_z_max-1.0)* &
        ((nzPML_2 - k - 0.5) / (nzPML_2 - 1.0))**m
        bh_z_2(k) = EXP(-(sigh_z_PML_2(k) / kappah_z_PML_2(k) + &
        alphah_z_PML_2(k))*del_t/ep_0)
        ch_z_2(k) = sigh_z_PML_2(k)*(bh_z_2(k)-1.0)/ &
        (sigh_z_PML_2(k)+kappah_z_PML_2(k)*alphah_z_PML_2(k)) &
        /kappah_z_PML_2(k)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup the denominatory arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !denominator arrays - used by regular space and PML

    allocate(den_ex(x_size-1))
    allocate(den_ey(y_size-1))
    allocate(den_ez(z_size-1))
    allocate(den_hx(x_size-1))
    allocate(den_hy(y_size-1))
    allocate(den_hz(z_size-1))

    !fully filled in below so not initialized

    !Now the denominator update equations
    ii = nxPML_2-1
    do i = 1,x_size-1
        if (i <= nxPML_1-1) then
            den_hx(i) = 1.0/(kappah_x_PML_1(i)*del_x)
        elseif (i >= x_size+1-nxPML_2) then
            den_hx(i) = 1.0/(kappah_x_PML_2(ii)*del_x)
            ii = ii-1
        else
            den_hx(i) = 1.0/del_x
        end if
    end do

    jj = nyPML_2-1
    do j = 1,y_size-1
        if (j <= nyPML_1-1) then
            den_hy(j) = 1.0/(kappah_y_PML_1(j)*del_y)
        elseif (j >= y_size+1-nyPML_2) then
            den_hy(j) = 1.0/(kappah_y_PML_2(jj)*del_y)
            jj = jj-1
        else
            den_hy(j) = 1.0/del_y
        end if
    end do

    kk = nzPML_2-1
    do k = 1,z_size-1
        if (k <= nzPML_1-1) then
            den_hz(k) = 1.0/(kappah_z_PML_1(k)*del_z)
        elseif (k >= z_size+1-nzPML_2) then
            den_hz(k) = 1.0/(kappah_z_PML_2(kk)*del_z)
            kk = kk - 1
        else
            den_hz(k) = 1.0/del_z
        end if
    end do

    ii = nxPML_2
    do i = 1,x_size-1
        if (i <= nxPML_1) then
            den_ex(i) = 1.0/(kappae_x_PML_1(i)*del_x)
        elseif (i >= x_size+1-nxPML_2) then
            den_ex(i) = 1.0/(kappae_x_PML_2(ii)*del_x)
            ii = ii-1
        else
            den_ex(i) = 1.0/del_x
        end if
    end do

    jj = nyPML_2
    do j = 1,y_size-1
        if (j <= nyPML_1) then
            den_ey(j) = 1.0/(kappae_y_PML_1(j)*del_y)
        elseif (j >= y_size+1-nyPML_2) then
            den_ey(j) = 1.0/(kappae_y_PML_2(jj)*del_y)
            jj = jj-1
        else
            den_ey(j) = 1.0/del_y
        end if
    end do

    kk = nzPML_2
    do k = 1,z_size-1
        if (k <= nzPML_1) then
            den_ez(k) = 1.0/(kappae_z_PML_1(k)*del_z)
        elseif (k >= z_size+1-nzPML_2) then
            den_ez(k) = 1.0/(kappae_z_PML_2(kk)*del_z)
            kk = kk - 1
        else
            den_ez(k) = 1.0/del_z
        end if
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Speed improvements for sheets in main fdtd section !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !to improve speed of sheet based loops in main fdtd loop, I will only scan through the user entered heights

    !for slight speed purposes, always declared, even if not needed, at full size at the cost of memory (pbc conditions)
    allocate(x_sheet_list(x_size-1))
    allocate(y_sheet_list(y_size-1))
    allocate(z_sheet_list(z_size-1))

    !initialize
    x_sheet_list(:)=0
    y_sheet_list(:)=0
    z_sheet_list(:)=0

    if (num_sheets_x>0) then
        count_unique_sheets_x=1
        do i=1, x_size-1
            if ( any(i == sheets_x(:,2,1)) ) then
                x_sheet_list(count_unique_sheets_x)=i
                count_unique_sheets_x=count_unique_sheets_x+1
            end if
        end do
        count_unique_sheets_x = count_unique_sheets_x - 1
    end if

    if (num_sheets_y>0) then
        count_unique_sheets_y=1
        do j=1, y_size-1
            if ( any(j == sheets_y(:,2,1)) ) then
                y_sheet_list(count_unique_sheets_y)=j
                count_unique_sheets_y=count_unique_sheets_y+1
            end if
        end do
        count_unique_sheets_y = count_unique_sheets_y - 1
    end if

    if (num_sheets_z>0) then
        count_unique_sheets_z=1
        do k=1, z_size-1
            if ( any(k == sheets_z(:,2,1)) ) then
                z_sheet_list(count_unique_sheets_z)=k
                count_unique_sheets_z=count_unique_sheets_z+1
            end if
        end do
        count_unique_sheets_z = count_unique_sheets_z - 1
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Setup gridded feeds if relevant !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    j=0
    do i=1, num_ports
        if (ports(i,5,1)==1) then
            j=j+1
        end if
    end do

    if (j>0) then
        allocate(gridded_feed_all(num_ports,4,MAXVAL(int(ports(:,4,1)))+1,MAXVAL(int(ports(:,4,2)))+1))
    else
        allocate(gridded_feed_all(1,1,1,1))
    end if

    do i=1,num_ports
        if (ports(i,5,1)==1) then

            open(unit=100+i, file=trim(gridded_feed_names(i)), access='stream', status='old')
                read(100+i) gridded_feed_all(i,:,1:int(ports(i,4,1))+1,1:int(ports(i,4,2))+1)
            close(100+i)  

        end if
    end do

#ifdef use_spice_version
    k=0
    do i=1, num_spice_ports
        if (ports_spice(i,5,1)==1) then
            k=k+1
        end if
    end do

    if (k>0) then
        allocate(gridded_feed_all_spice(num_spice_ports,4,MAXVAL(int(ports_spice(:,3,1)))+1,MAXVAL(int(ports_spice(:,3,2)))+1))
    else
        allocate(gridded_feed_all_spice(1,1,1,1))
    end if

    do i=1,num_spice_ports
        if (ports_spice(i,5,1)==1) then

            open(unit=100+i+num_ports, file=trim(gridded_feed_names_spice(i)), access='stream', status='old')
                read(100+i+num_ports) gridded_feed_all_spice(i,:,1:int(ports_spice(i,3,1))+1,1:int(ports_spice(i,3,2))+1)
            close(100+i+num_ports)  

        end if
    end do
#endif

#ifdef use_spice_version
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Added for Spice Setup !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !allocate and intialize all arrays and parameters needed
    !similar to above, we need to allocate regardless of use for acc passage.

    if (num_spice_ports>0) then
        allocate(Spice_Current(num_spice_ports))
        allocate(Spice_Voltage(num_spice_ports))
        allocate(Spice_Voltage_out(num_spice_ports,time_steps))
        allocate(Spice_Current_out(num_spice_ports,time_steps))
        allocate(names(num_spice_ports))
        allocate(Cap(num_spice_ports))
    else
        allocate(Spice_Current(1))
        allocate(Spice_Voltage(1))
        allocate(Spice_Voltage_out(1,1))
        allocate(Spice_Current_out(1,1))
        allocate(names(1))
        allocate(Cap(1))
    end if

    Spice_Current(:)=0.0+0.0*imag
    Spice_Current_out(:,:)=0.0+0.0*imag
    Spice_Voltage(:)=0.0
    Spice_Voltage_out(:,:)=0.0
    Cap(:)=0.0
    do i=1, num_spice_ports
        names(i)%name = TRIM(names_of_spice_ports(i))
    end do

    !now setup some items
    circuit%time = 0.0
    circuit%dt = del_t
    finalTime = del_t*time_steps
    error_cnt = 0

    if (num_spice_ports>0) then
        ! --- Initialize the circuit ---
        call circuit%init(names=names, netlist=netlist)
        call circuit%setStopTimes(finalTime, circuit%dt)
        !now that it's intialized, we need to update/correct the capacitances for each basic port before starting the time loop
        !we will perform an average of them and normalize to gap distance like current is
        !so average ep * (del*del*non_gap_dim*non_gap_dim)/(gap_dim*del)/gap_dim
        do i=1, num_spice_ports
            if (ports_spice(i,5,1)==0) then          
                !these are the starting and stopping index locations
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                p_loc_d(1)=ports_spice(i,3,1)+p_loc(1)-1
                p_loc_d(2)=ports_spice(i,3,2)+p_loc(2)-1
                p_loc_d(3)=ports_spice(i,3,3)+p_loc(3)-1

                !set Cap(i) to zero to initalize it
                Cap(i)=0
                
                select case (int(ports_spice(i,1,1)))
                case (0)
                    !sum the permittivities
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Cap(i)=Cap(i)+relative_ep_x(ii,jj,kk)*ep_0
                            end do
                        end do
                    end do
                    !normalize to get the actual average
                    Cap(i)=Cap(i)/(ports_spice(i,3,1)*ports_spice(i,3,2)*ports_spice(i,3,3))
                    !now calculate using the formula for a parallel plate
                    Cap(i)=Cap(i)*(ports_spice(i,3,2)*del_y*ports_spice(i,3,3)*del_z)/(ports_spice(i,3,1)*del_x)
                    !normalize to gap distance to be consistnet with currents
                    Cap(i)=Cap(i)/ports_spice(i,3,1)

                case (1)
                    !sum the permittivities
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Cap(i)=Cap(i)+relative_ep_y(ii,jj,kk)*ep_0
                            end do
                        end do
                    end do
                    !normalize to get the actual average
                    Cap(i)=Cap(i)/(ports_spice(i,3,1)*ports_spice(i,3,2)*ports_spice(i,3,3))
                    !now calculate using the formula for a parallel plate
                    Cap(i)=Cap(i)*(ports_spice(i,3,1)*del_x*ports_spice(i,3,3)*del_z)/(ports_spice(i,3,2)*del_y)
                    !normalize to gap distance to be consistnet with currents
                    Cap(i)=Cap(i)/ports_spice(i,3,2)

                case (2)
                    !sum the permittivities
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Cap(i)=Cap(i)+relative_ep_z(ii,jj,kk)*ep_0
                            end do
                        end do
                    end do
                    !normalize to get the actual average
                    Cap(i)=Cap(i)/(ports_spice(i,3,1)*ports_spice(i,3,2)*ports_spice(i,3,3))
                    !now calculate using the formula for a parallel plate
                    Cap(i)=Cap(i)*(ports_spice(i,3,1)*del_x*ports_spice(i,3,2)*del_y)/(ports_spice(i,3,3)*del_z)
                    !normalize to gap distance to be consistnet with currents
                    Cap(i)=Cap(i)/ports_spice(i,3,3)

                end select

                command_str = 'alter @'//trim(C_name(i))//'='
                write(pair_str, '(E15.8)') Cap(i)
                command_str = trim(command_str) // trim(pair_str)
                ! --- Execute the command to update the Capacitance values associated with fdtd
                call command(trim(command_str) // char(0))
            end if
        end do
    end if

    !spice version is unique for deallocating relative permittivity - had to wait for Cap updates
    deallocate(relative_ep_x, &
            relative_ep_y, &
            relative_ep_z)
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Allocate and Initialize the rest of the arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !allocate and initalize the rest of the arrays needed
    !these are the actul fields of interest (E,H,V,I,etc.)

    allocate(Ex(x_size-1,y_size,z_size))
    allocate(Ey(x_size,y_size-1,z_size))
    allocate(Ez(x_size,y_size,z_size-1))
    allocate(Hx(x_size,y_size-1,z_size-1))
    allocate(Hy(x_size-1,y_size,z_size-1))
    allocate(Hz(x_size-1,y_size-1,z_size))

    !added for any material auxillary differential equations
    if (num_poles>0) then
        allocate(J_source_x(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(J_source_y(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(J_source_z(num_poles,x_size-1,y_size-1,z_size-1))
        allocate(Ex_oldt(x_size-1,y_size,z_size))
        allocate(Ey_oldt(x_size,y_size-1,z_size))
        allocate(Ez_oldt(x_size,y_size,z_size-1))
    else 
        allocate(J_source_x(1,1,1,1))
        allocate(J_source_y(1,1,1,1))
        allocate(J_source_z(1,1,1,1))
        allocate(Ex_oldt(1,1,1))
        allocate(Ey_oldt(1,1,1))
        allocate(Ez_oldt(1,1,1))
    end if

    !voltage and current aux arrays
    !I don't use the curret_out arrays but it's a small overhead in memory and I might want them later on
    if (num_ports>0) then
        allocate(Voltage(num_ports,port_array_size))
        allocate(Current(num_ports,port_array_size-1))
        allocate(Voltage_out(num_ports,time_steps))
        allocate(Current_out(num_ports,time_steps))
    else 
        allocate(Voltage(1,1))
        allocate(Current(1,1))
        allocate(Voltage_out(1,1))
        allocate(Current_out(1,1))
    end if
    if (num_ports>0 .and. antenna_amp>0) then
        allocate(V_inc(port_array_size*2))
        allocate(C_inc(port_array_size*2-1))
    else
        allocate(V_inc(1))
        allocate(C_inc(1))
    end if

    !video arrays
    if (video_on==1) then
        allocate(Ex_video(vid_size1,vid_size2,time_steps))
        allocate(Hx_video(vid_size1,vid_size2,time_steps))
        allocate(Ey_video(vid_size1,vid_size2,time_steps))
        allocate(Hy_video(vid_size1,vid_size2,time_steps))
        allocate(Ez_video(vid_size1,vid_size2,time_steps))
        allocate(Hz_video(vid_size1,vid_size2,time_steps))
    else 
        allocate(Ex_video(1,1,1))
        allocate(Hx_video(1,1,1))
        allocate(Ey_video(1,1,1))
        allocate(Hy_video(1,1,1))
        allocate(Ez_video(1,1,1))
        allocate(Hz_video(1,1,1))
    end if

    !incident field we save to use in post processing
    allocate(incident(time_steps))
    !output field arrays for post processing later on - these are for pbc s parameters
#ifndef use_kmax_version
    if (pbc_x+pbc_y+pbc_z==2) then
        allocate(E_reflected(time_steps))
        allocate(E_transmitted(time_steps))
    else
        allocate(E_reflected(1))
        allocate(E_transmitted(1))
    end if
#endif
#ifdef use_kmax_version
    allocate(E_reflected_TE(time_steps))
    allocate(E_reflected_TM(time_steps))
    allocate(E_transmitted_TE(time_steps))
    allocate(E_transmitted_TM(time_steps))
#endif

    !for slight speed purposes, special are always declared full size at the cost of memory (pbc conditions)
    allocate(Ex_special(x_size-1,y_size,z_size))
    allocate(Ey_special(x_size,y_size-1,z_size))
    allocate(Ez_special(x_size,y_size,z_size-1))

    !These are for far field post processing
    if (num_far_field_angles>0) then
        if (pbc_x==0) then
            allocate(My_xlow(ff_yhigh,ff_zhigh))
            allocate(Mz_xlow(ff_yhigh,ff_zhigh))
            allocate(Jy_xlow(ff_yhigh,ff_zhigh))
            allocate(Jz_xlow(ff_yhigh,ff_zhigh))
            allocate(My_xlow_oldt(ff_yhigh,ff_zhigh))
            allocate(Mz_xlow_oldt(ff_yhigh,ff_zhigh))
            allocate(Jy_xlow_oldt(ff_yhigh,ff_zhigh))
            allocate(Jz_xlow_oldt(ff_yhigh,ff_zhigh))
            allocate(My_xhigh(ff_yhigh,ff_zhigh))
            allocate(Mz_xhigh(ff_yhigh,ff_zhigh))
            allocate(Jy_xhigh(ff_yhigh,ff_zhigh))
            allocate(Jz_xhigh(ff_yhigh,ff_zhigh))
            allocate(My_xhigh_oldt(ff_yhigh,ff_zhigh))
            allocate(Mz_xhigh_oldt(ff_yhigh,ff_zhigh))
            allocate(Jy_xhigh_oldt(ff_yhigh,ff_zhigh))
            allocate(Jz_xhigh_oldt(ff_yhigh,ff_zhigh))
        else
            allocate(My_xlow(1,1))
            allocate(Mz_xlow(1,1))
            allocate(Jy_xlow(1,1))
            allocate(Jz_xlow(1,1))
            allocate(My_xlow_oldt(1,1))
            allocate(Mz_xlow_oldt(1,1))
            allocate(Jy_xlow_oldt(1,1))
            allocate(Jz_xlow_oldt(1,1))
            allocate(My_xhigh(1,1))
            allocate(Mz_xhigh(1,1))
            allocate(Jy_xhigh(1,1))
            allocate(Jz_xhigh(1,1))
            allocate(My_xhigh_oldt(1,1))
            allocate(Mz_xhigh_oldt(1,1))
            allocate(Jy_xhigh_oldt(1,1))
            allocate(Jz_xhigh_oldt(1,1))
        end if
        if (pbc_y==0) then
            allocate(Mx_ylow(ff_xhigh,ff_zhigh))
            allocate(Mz_ylow(ff_xhigh,ff_zhigh))
            allocate(Jx_ylow(ff_xhigh,ff_zhigh))
            allocate(Jz_ylow(ff_xhigh,ff_zhigh))
            allocate(Mx_ylow_oldt(ff_xhigh,ff_zhigh))
            allocate(Mz_ylow_oldt(ff_xhigh,ff_zhigh))
            allocate(Jx_ylow_oldt(ff_xhigh,ff_zhigh))
            allocate(Jz_ylow_oldt(ff_xhigh,ff_zhigh))
            allocate(Mx_yhigh(ff_xhigh,ff_zhigh))
            allocate(Mz_yhigh(ff_xhigh,ff_zhigh))
            allocate(Jx_yhigh(ff_xhigh,ff_zhigh))
            allocate(Jz_yhigh(ff_xhigh,ff_zhigh))
            allocate(Mx_yhigh_oldt(ff_xhigh,ff_zhigh))
            allocate(Mz_yhigh_oldt(ff_xhigh,ff_zhigh))
            allocate(Jx_yhigh_oldt(ff_xhigh,ff_zhigh))
            allocate(Jz_yhigh_oldt(ff_xhigh,ff_zhigh))
        else
            allocate(Mx_ylow(1,1))
            allocate(Mz_ylow(1,1))
            allocate(Jx_ylow(1,1))
            allocate(Jz_ylow(1,1))
            allocate(Mx_ylow_oldt(1,1))
            allocate(Mz_ylow_oldt(1,1))
            allocate(Jx_ylow_oldt(1,1))
            allocate(Jz_ylow_oldt(1,1))
            allocate(Mx_yhigh(1,1))
            allocate(Mz_yhigh(1,1))
            allocate(Jx_yhigh(1,1))
            allocate(Jz_yhigh(1,1))
            allocate(Mx_yhigh_oldt(1,1))
            allocate(Mz_yhigh_oldt(1,1))
            allocate(Jx_yhigh_oldt(1,1))
            allocate(Jz_yhigh_oldt(1,1))
        end if
        if (pbc_z==0) then
            allocate(Mx_zlow(ff_xhigh,ff_yhigh))
            allocate(My_zlow(ff_xhigh,ff_yhigh))
            allocate(Jx_zlow(ff_xhigh,ff_yhigh))
            allocate(Jy_zlow(ff_xhigh,ff_yhigh))
            allocate(Mx_zlow_oldt(ff_xhigh,ff_yhigh))
            allocate(My_zlow_oldt(ff_xhigh,ff_yhigh))
            allocate(Jx_zlow_oldt(ff_xhigh,ff_yhigh))
            allocate(Jy_zlow_oldt(ff_xhigh,ff_yhigh))
            allocate(Mx_zhigh(ff_xhigh,ff_yhigh))
            allocate(My_zhigh(ff_xhigh,ff_yhigh))
            allocate(Jx_zhigh(ff_xhigh,ff_yhigh))
            allocate(Jy_zhigh(ff_xhigh,ff_yhigh))
            allocate(Mx_zhigh_oldt(ff_xhigh,ff_yhigh))
            allocate(My_zhigh_oldt(ff_xhigh,ff_yhigh))
            allocate(Jx_zhigh_oldt(ff_xhigh,ff_yhigh))
            allocate(Jy_zhigh_oldt(ff_xhigh,ff_yhigh))
        else
            allocate(Mx_zlow(1,1))
            allocate(My_zlow(1,1))
            allocate(Jx_zlow(1,1))
            allocate(Jy_zlow(1,1))
            allocate(Mx_zlow_oldt(1,1))
            allocate(My_zlow_oldt(1,1))
            allocate(Jx_zlow_oldt(1,1))
            allocate(Jy_zlow_oldt(1,1))
            allocate(Mx_zhigh(1,1))
            allocate(My_zhigh(1,1))
            allocate(Jx_zhigh(1,1))
            allocate(Jy_zhigh(1,1))
            allocate(Mx_zhigh_oldt(1,1))
            allocate(My_zhigh_oldt(1,1))
            allocate(Jx_zhigh_oldt(1,1))
            allocate(Jy_zhigh_oldt(1,1))
        end if
        !W and U is how we add fields together in time
        allocate(Wx(num_far_field_angles,len_far_field_arrays))
        allocate(Wy(num_far_field_angles,len_far_field_arrays))
        allocate(Wz(num_far_field_angles,len_far_field_arrays))
        allocate(Ux(num_far_field_angles,len_far_field_arrays))
        allocate(Uy(num_far_field_angles,len_far_field_arrays))
        allocate(Uz(num_far_field_angles,len_far_field_arrays))
        allocate(W_theta(num_far_field_angles,len_far_field_arrays))
        allocate(W_phi(num_far_field_angles,len_far_field_arrays))
        allocate(U_theta(num_far_field_angles,len_far_field_arrays))
        allocate(U_phi(num_far_field_angles,len_far_field_arrays))
        !Then we combine U and W in theta,phi format to E outputs
        !one of these will be cross pol term
        allocate(E_theta(num_far_field_angles,len_far_field_arrays))
        allocate(E_phi(num_far_field_angles,len_far_field_arrays))
        allocate(E_theta_out(num_far_field_angles,time_steps))
        allocate(E_phi_out(num_far_field_angles,time_steps))
        allocate(data_out_time(num_far_field_angles))
    else
        allocate(My_xlow(1,1))
        allocate(Mz_xlow(1,1))
        allocate(Jy_xlow(1,1))
        allocate(Jz_xlow(1,1))
        allocate(My_xlow_oldt(1,1))
        allocate(Mz_xlow_oldt(1,1))
        allocate(Jy_xlow_oldt(1,1))
        allocate(Jz_xlow_oldt(1,1))
        allocate(My_xhigh(1,1))
        allocate(Mz_xhigh(1,1))
        allocate(Jy_xhigh(1,1))
        allocate(Jz_xhigh(1,1))
        allocate(My_xhigh_oldt(1,1))
        allocate(Mz_xhigh_oldt(1,1))
        allocate(Jy_xhigh_oldt(1,1))
        allocate(Jz_xhigh_oldt(1,1))
        allocate(Mx_ylow(1,1))
        allocate(Mz_ylow(1,1))
        allocate(Jx_ylow(1,1))
        allocate(Jz_ylow(1,1))
        allocate(Mx_ylow_oldt(1,1))
        allocate(Mz_ylow_oldt(1,1))
        allocate(Jx_ylow_oldt(1,1))
        allocate(Jz_ylow_oldt(1,1))
        allocate(Mx_yhigh(1,1))
        allocate(Mz_yhigh(1,1))
        allocate(Jx_yhigh(1,1))
        allocate(Jz_yhigh(1,1))
        allocate(Mx_yhigh_oldt(1,1))
        allocate(Mz_yhigh_oldt(1,1))
        allocate(Jx_yhigh_oldt(1,1))
        allocate(Jz_yhigh_oldt(1,1))
        allocate(Mx_zlow(1,1))
        allocate(My_zlow(1,1))
        allocate(Jx_zlow(1,1))
        allocate(Jy_zlow(1,1))
        allocate(Mx_zlow_oldt(1,1))
        allocate(My_zlow_oldt(1,1))
        allocate(Jx_zlow_oldt(1,1))
        allocate(Jy_zlow_oldt(1,1))
        allocate(Mx_zhigh(1,1))
        allocate(My_zhigh(1,1))
        allocate(Jx_zhigh(1,1))
        allocate(Jy_zhigh(1,1))
        allocate(Mx_zhigh_oldt(1,1))
        allocate(My_zhigh_oldt(1,1))
        allocate(Jx_zhigh_oldt(1,1))
        allocate(Jy_zhigh_oldt(1,1))
        !W and U is how we add fields together in time
        allocate(Wx(1,1))
        allocate(Wy(1,1))
        allocate(Wz(1,1))
        allocate(Ux(1,1))
        allocate(Uy(1,1))
        allocate(Uz(1,1))
        allocate(W_theta(1,1))
        allocate(W_phi(1,1))
        allocate(U_theta(1,1))
        allocate(U_phi(1,1))
        !Then we combine U and W in theta,phi format to E outputs
        !one of these will be cross pol term
        allocate(E_theta(1,1))
        allocate(E_phi(1,1))
        allocate(E_theta_out(1,1))
        allocate(E_phi_out(1,1))
        allocate(data_out_time(1))
    end if

    !these are for far field phase centering if inc plane waves are used
    !currently we only use E and not H, but that could change and it's a small memory overhead
    if (num_far_field_angles>0 .and. plane_wave_amp>0) then
        allocate(Ex_ff_pc(time_steps))
        allocate(Ey_ff_pc(time_steps))
        allocate(Ez_ff_pc(time_steps))
        allocate(Hx_ff_pc(time_steps))
        allocate(Hy_ff_pc(time_steps))
        allocate(Hz_ff_pc(time_steps))
    else
        allocate(Ex_ff_pc(1))
        allocate(Ey_ff_pc(1))
        allocate(Ez_ff_pc(1))
        allocate(Hx_ff_pc(1))
        allocate(Hy_ff_pc(1))
        allocate(Hz_ff_pc(1))
    end if

    !now initialize all of these allocated arrays in the same approx. order
    Ez(:,:,:) = 0.0+0.0*imag
    Hz(:,:,:) = 0.0+0.0*imag
    Ex(:,:,:) = 0.0+0.0*imag
    Hx(:,:,:) = 0.0+0.0*imag
    Ey(:,:,:) = 0.0+0.0*imag
    Hy(:,:,:) = 0.0+0.0*imag

    J_source_x(:,:,:,:) = 0.0+0.0*imag
    J_source_y(:,:,:,:) = 0.0+0.0*imag
    J_source_z(:,:,:,:) = 0.0+0.0*imag
    Ex_oldt(:,:,:) = 0.0+0.0*imag
    Ey_oldt(:,:,:) = 0.0+0.0*imag
    Ez_oldt(:,:,:) = 0.0+0.0*imag

    Voltage(:,:)=0.0+0.0*imag
    Current(:,:)=0.0+0.0*imag
    Voltage_out(:,:)=0.0+0.0*imag
    Current_out(:,:)=0.0+0.0*imag
    V_inc(:)=0.0+0.0*imag
    C_inc(:)=0.0+0.0*imag

    p_loc(:)=0.0
    p_loc_d(:)=0.0

    Ex_video(:,:,:)=0.0+0.0*imag
    Hx_video(:,:,:)=0.0+0.0*imag
    Ey_video(:,:,:)=0.0+0.0*imag
    Hy_video(:,:,:)=0.0+0.0*imag
    Ez_video(:,:,:)=0.0+0.0*imag
    Hz_video(:,:,:)=0.0+0.0*imag

    incident(:)=0.0+0.0*imag
#ifndef use_kmax_version
    E_reflected(:)=0.0
    E_transmitted(:)=0.0
#endif
#ifdef use_kmax_version
    E_reflected_TE(:)=0.0+0.0*imag
    E_transmitted_TE(:)=0.0+0.0*imag
    E_reflected_TM(:)=0.0+0.0*imag
    E_transmitted_TM(:)=0.0+0.0*imag
#endif

    Ex_special(:,:,:)=0.0+0.0*imag
    Ey_special(:,:,:)=0.0+0.0*imag
    Ez_special(:,:,:)=0.0+0.0*imag

    My_xlow(:,:)=0.0+0.0*imag
    Mz_xlow(:,:)=0.0+0.0*imag
    Jy_xlow(:,:)=0.0+0.0*imag
    Jz_xlow(:,:)=0.0+0.0*imag
    My_xlow_oldt(:,:)=0.0+0.0*imag
    Mz_xlow_oldt(:,:)=0.0+0.0*imag
    Jy_xlow_oldt(:,:)=0.0+0.0*imag
    Jz_xlow_oldt(:,:)=0.0+0.0*imag
    My_xhigh(:,:)=0.0+0.0*imag
    Mz_xhigh(:,:)=0.0+0.0*imag
    Jy_xhigh(:,:)=0.0+0.0*imag
    Jz_xhigh(:,:)=0.0+0.0*imag    
    My_xhigh_oldt(:,:)=0.0+0.0*imag
    Mz_xhigh_oldt(:,:)=0.0+0.0*imag
    Jy_xhigh_oldt(:,:)=0.0+0.0*imag
    Jz_xhigh_oldt(:,:)=0.0+0.0*imag
    Mx_ylow(:,:)=0.0+0.0*imag
    Mz_ylow(:,:)=0.0+0.0*imag
    Jx_ylow(:,:)=0.0+0.0*imag
    Jz_ylow(:,:)=0.0+0.0*imag    
    Mx_ylow_oldt(:,:)=0.0+0.0*imag
    Mz_ylow_oldt(:,:)=0.0+0.0*imag
    Jx_ylow_oldt(:,:)=0.0+0.0*imag
    Jz_ylow_oldt(:,:)=0.0+0.0*imag
    Mx_yhigh(:,:)=0.0+0.0*imag
    Mz_yhigh(:,:)=0.0+0.0*imag
    Jx_yhigh(:,:)=0.0+0.0*imag
    Jz_yhigh(:,:)=0.0+0.0*imag    
    Mx_yhigh_oldt(:,:)=0.0+0.0*imag
    Mz_yhigh_oldt(:,:)=0.0+0.0*imag
    Jx_yhigh_oldt(:,:)=0.0+0.0*imag
    Jz_yhigh_oldt(:,:)=0.0+0.0*imag
    Mx_zlow(:,:)=0.0+0.0*imag
    My_zlow(:,:)=0.0+0.0*imag
    Jx_zlow(:,:)=0.0+0.0*imag
    Jy_zlow(:,:)=0.0+0.0*imag  
    Mx_zlow_oldt(:,:)=0.0+0.0*imag
    My_zlow_oldt(:,:)=0.0+0.0*imag
    Jx_zlow_oldt(:,:)=0.0+0.0*imag
    Jy_zlow_oldt(:,:)=0.0+0.0*imag
    Mx_zhigh(:,:)=0.0+0.0*imag
    My_zhigh(:,:)=0.0+0.0*imag
    Jx_zhigh(:,:)=0.0+0.0*imag
    Jy_zhigh(:,:)=0.0+0.0*imag    
    Mx_zhigh_oldt(:,:)=0.0+0.0*imag
    My_zhigh_oldt(:,:)=0.0+0.0*imag
    Jx_zhigh_oldt(:,:)=0.0+0.0*imag
    Jy_zhigh_oldt(:,:)=0.0+0.0*imag
    Wx(:,:)=0.0+0.0*imag
    Wy(:,:)=0.0+0.0*imag
    Wz(:,:)=0.0+0.0*imag
    Ux(:,:)=0.0+0.0*imag
    Uy(:,:)=0.0+0.0*imag
    Uz(:,:)=0.0+0.0*imag
    W_theta(:,:)=0.0+0.0*imag
    W_Phi(:,:)=0.0+0.0*imag
    U_theta(:,:)=0.0+0.0*imag
    U_Phi(:,:)=0.0+0.0*imag
    E_theta(:,:)=0.0+0.0*imag
    E_phi(:,:)=0.0+0.0*imag
    E_theta_out(:,:)=0.0+0.0*imag
    E_phi_out(:,:)=0.0+0.0*imag
    data_out_time(:)=1
    Ex_ff_pc(:)=0.0+0.0*imag
    Ey_ff_pc(:)=0.0+0.0*imag
    Ez_ff_pc(:)=0.0+0.0*imag
    Hx_ff_pc(:)=0.0+0.0*imag
    Hy_ff_pc(:)=0.0+0.0*imag
    Hz_ff_pc(:)=0.0+0.0*imag
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! END ALL SETUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_end)
    write(*, '("Setup is finished. Time to complete was ", F0.3, " seconds.")') &
          real(clock_time_end - clock_time_start) / real(clock_rate)
    write(*, '("Starting main FDTD algorithm...")')
    call system_clock(clock_time_start)    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! BEGIN MAIN FDTD ALGORITHM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !flags for acc or openmp are only used and recongized if we use the right compilation flags
    !$acc data copy(Hx, Hy, Hz, Ex, Ey, Ez) &
    !$acc      copy(Ex_video, Hx_video, Ey_video, Hy_video, Ez_video, Hz_video) &
#ifndef use_kmax_version
    !$acc      copy(E_reflected, E_transmitted, incident) &
#endif
#ifdef use_kmax_version
    !$acc      copy(E_reflected_TE, E_transmitted_TE,E_reflected_TM, E_transmitted_TM, incident) &
#endif
    !$acc      copy(Ux,Uy,Uz,Wx,Wy,Wz) &
    !$acc      copy(Voltage_out, Voltage) &
    !$acc      copy(Current_out, Current) &
#ifdef use_spice_version
    !$acc      copy(Spice_Current,Spice_Current_out,Spice_Voltage,Spice_Voltage_out,ports_spice,names_of_spice_ports_currents) &
#endif
    !$acc      copy(Ex_ff_pc,Ey_ff_pc,Ez_ff_pc,Hx_ff_pc,Hy_ff_pc,Hz_ff_pc) &
    !$acc      copyin(My_xlow, My_xlow_oldt, Mz_xlow, Mz_xlow_oldt, Jy_xlow, Jy_xlow_oldt, Jz_xlow, Jz_xlow_oldt) &
    !$acc      copyin(My_xhigh, My_xhigh_oldt, Mz_xhigh, Mz_xhigh_oldt, Jy_xhigh, Jy_xhigh_oldt, Jz_xhigh, Jz_xhigh_oldt) &
    !$acc      copyin(Mx_ylow, Mx_ylow_oldt, Mz_ylow, Mz_ylow_oldt, Jx_ylow, Jx_ylow_oldt, Jz_ylow, Jz_ylow_oldt) &
    !$acc      copyin(Mx_yhigh, Mx_yhigh_oldt, Mz_yhigh, Mz_yhigh_oldt, Jx_yhigh, Jx_yhigh_oldt, Jz_yhigh, Jz_yhigh_oldt) &
    !$acc      copyin(Mx_zlow, Mx_zlow_oldt, My_zlow, My_zlow_oldt, Jx_zlow, Jx_zlow_oldt, Jy_zlow, Jy_zlow_oldt) &
    !$acc      copyin(Mx_zhigh, Mx_zhigh_oldt, My_zhigh, My_zhigh_oldt, Jx_zhigh, Jx_zhigh_oldt, Jy_zhigh, Jy_zhigh_oldt) &
    !$acc      copyin(far_field_angles) &
    !$acc      copyin(p_loc,p_loc_d,ports) &
    !$acc      copyin(psi_Hxy_1, psi_Hxy_2, psi_Hxz_1, psi_Hxz_2) &
    !$acc      copyin(psi_Hyx_1, psi_Hyx_2, psi_Hyz_1, psi_Hyz_2) &
    !$acc      copyin(psi_Hzx_1, psi_Hzx_2, psi_Hzy_1, psi_Hzy_2) &
    !$acc      copyin(psi_Exy_1, psi_Exy_2, psi_Exz_1, psi_Exz_2) &
    !$acc      copyin(psi_Eyx_1, psi_Eyx_2, psi_Eyz_1, psi_Eyz_2) &
    !$acc      copyin(psi_Ezx_1, psi_Ezx_2, psi_Ezy_1, psi_Ezy_2) &
    !$acc      copyin(den_hy, den_hz, den_hx, den_ez, den_ex, den_ey) &
    !$acc      copyin(gax, gay, gaz, gbx, gby, gbz) &
    !$acc      copyin(J_plasma_ax,J_plasma_bx,J_plasma_ay,J_plasma_by,J_plasma_az,J_plasma_bz) &
    !$acc      copyin(J_source_x,J_source_y,J_source_z) &
    !$acc      copyin(Ex_oldt,Ey_oldt,Ez_oldt) &
    !$acc      copyin(bh_y_1, ch_y_1, bh_y_2, ch_y_2) &
    !$acc      copyin(bh_z_1, ch_z_1, bh_z_2, ch_z_2) &
    !$acc      copyin(bh_x_1, ch_x_1, bh_x_2, ch_x_2) & 
    !$acc      copyin(be_y_1, ce_y_1, be_y_2, ce_y_2) &
    !$acc      copyin(be_z_1, ce_z_1, be_z_2, ce_z_2) &
    !$acc      copyin(be_x_1, ce_x_1, be_x_2, ce_x_2) &
    !$acc      copyin(Ex_special, Ey_special, Ez_special) &
    !$acc      copyin(sheet_sig_x_x, sheet_ep_x_x, sheet_sig_y_y, sheet_ep_y_y, sheet_sig_z_z, sheet_ep_z_z) &
    !$acc      copyin(x_sheet_list, y_sheet_list, z_sheet_list) &
#ifdef use_spice_version
    !$acc      copyin(gridded_feed_all) &
    !$acc      copyin(gridded_feed_all_spice)
#endif
#ifndef use_spice_version
    !$acc      copyin(gridded_feed_all)
#endif
    !$omp parallel default(shared) private(i, j, k, ii, jj, kk, ang1, ang2, jj_real, time_f_var, step_vc_out, &
#ifndef use_kmax_version
    !$omp plasma_counter, temp_plasma_x, temp_plasma_y, temp_plasma_z, i_mirror, j_mirror, k_mirror)
#endif
#ifdef use_kmax_version
    !$omp plasma_counter, temp_plasma_x, temp_plasma_y, temp_plasma_z, i_mirror, j_mirror, k_mirror, &
    !$omp temp_plasma_x_c,temp_plasma_y_c,temp_plasma_z_c)
#endif

    do counter=1,  (time_steps)

        !!UPDATE HX!!
        !$acc parallel loop collapse(3) present(Hx, Ez, Ey, den_hy, den_hz)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do j = 1,y_size-1
                do i = 1,x_size-1
                    Hx(i,j,k) = DA * Hx(i,j,k) + &
                    DB*((Ez(i,j,k)-Ez(i,j+1,k))*den_hy(j)+(Ey(i,j,k+1) - Ey(i,j,k))*den_hz(k))
                end do
            end do 
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!UPDATE HY!!
        !$acc parallel loop collapse(3) present(Hy, Ez, Ex, den_hx, den_hz)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do j = 1,y_size-1
                do i = 1,x_size-1
                    Hy(i,j,k) = DA * Hy(i,j,k) + DB * &
                    ( (Ez(i+1,j,k) - Ez(i,j,k))*den_hx(i) + &
                    (Ex(i,j,k) - Ex(i,j,k+1))*den_hz(k) )
                end do 
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!UPDATE HZ!!
        !$acc parallel loop collapse(3) present(Hz, Ey, Ex, den_hx, den_hy)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do j = 1,y_size-1
                do i = 1,x_size-1
                    Hz(i,j,k) = DA * Hz(i,j,k) + DB &
                    * ((Ey(i,j,k) - Ey(i+1,j,k))*den_hx(i) + &
                    (Ex(i,j+1,k) - Ex(i,j,k))*den_hy(j))
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!SHEET IMPEDANCES!! - special sub-cell update        
        !!X-SHEETS!!
        !$acc parallel loop collapse(3) present(Hy, Hz, Ex, Ex_special, sheet_sig_x_x, sheet_ep_x_x, den_hz, den_hy, x_sheet_list)
        !$omp do collapse(2) schedule(static)
        do ii = 1, count_unique_sheets_x
            do k = 1, z_size-1
                do j = 1, y_size-1
                    i = x_sheet_list(ii)
                    
                    !x sheets if they exist
                    if ((sheet_sig_x_x(i,j,k) > 0) .or. (sheet_ep_x_x(i,j,k) > 1)) then 
                        Hy(i,j,k) = Hy(i,j,k) + DB*(sheet_thickness/del_x)* &
                        ((Ex(i,j,k+1)-Ex(i,j,k)+Ex_special(i,j,k)-Ex_special(i,j,k+1))* &
                        den_hz(k))
                        Hz(i,j,k) = Hz(i,j,k) + DB*(sheet_thickness/del_x)* &
                        ((Ex(i,j,k)-Ex(i,j+1,k)+Ex_special(i,j+1,k)-Ex_special(i,j,k))* &
                        den_hy(j))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!Y-SHEETS!!
        !$acc parallel loop collapse(3) present(Hx, Hz, Ey, Ey_special, sheet_sig_y_y, sheet_ep_y_y, den_hz, den_hx, y_sheet_list)
        !$omp do collapse(2) schedule(static)
        do jj = 1, count_unique_sheets_y
            do k = 1, z_size-1
                do i = 1, x_size-1
                    j = y_sheet_list(jj)
                    
                    !y sheets if they exist
                    if ((sheet_sig_y_y(i,j,k) > 0) .or. (sheet_ep_y_y(i,j,k) > 1)) then
                        Hx(i,j,k) = Hx(i,j,k) + DB*(sheet_thickness/del_y)* &
                        (Ey(i,j,k)-Ey(i,j,k+1)+Ey_special(i,j,k+1)-Ey_special(i,j,k))* &
                        den_hz(k)
                        Hz(i,j,k) = Hz(i,j,k) + DB*(sheet_thickness/del_y)* &
                        ((Ey(i+1,j,k)-Ey(i,j,k)-Ey_special(i+1,j,k)+Ey_special(i,j,k))* &
                        den_hx(i))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!Z-SHEETS!!
        !$acc parallel loop collapse(3) present(Hx, Hy, Ez, Ez_special, sheet_sig_z_z, sheet_ep_z_z, den_hy, den_hx, z_sheet_list)
        !$omp do collapse(2) schedule(static)
        do kk = 1, count_unique_sheets_z
            do j = 1, y_size-1
                do i = 1, x_size-1
                    k = z_sheet_list(kk)
                    
                    !z sheets if they exist
                    if ((sheet_sig_z_z(i,j,k) > 0) .or. (sheet_ep_z_z(i,j,k) > 1)) then
                        Hx(i,j,k) = Hx(i,j,k) + DB*(sheet_thickness/del_z)* &
                        (Ez(i,j+1,k)-Ez(i,j,k)-Ez_special(i,j+1,k)+Ez_special(i,j,k))* &
                        den_hy(j)
                        Hy(i,j,k) = Hy(i,j,k) + DB*(sheet_thickness/del_z)* &
                        ((Ez(i,j,k)-Ez(i+1,j,k)+Ez_special(i+1,j,k)-Ez_special(i,j,k))* &
                        den_hx(i))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do


        !!INJECT H PLANE WAVE SOURCES FOR TF/SF FORMULATION (if not kmax) INTO THE GRID!!
#ifdef use_kmax_version
        !if TE mode then add E fields to H fields - no longer TF/SF anymore
        !1/2 values still used so that TE and TM get launched from the same physical plane (I think) - though this is uncessary
        !differs from main fdtd program in terms of starting height locations for both E/H
        if ((mode_type==0) .and. (plane_wave_amp>0)) then

            !!X FACES!!
            if (pbc_y+pbc_z==2) then
                !$acc parallel loop collapse(2) present(Hz,Hy)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh
                    do j=ylow,yhigh
                        Hz(k_pl_start_H,j,k)=Hz(k_pl_start_H,j,k)+del_t/(mu_0*del_x)*& 
                        Inc(WEy,0.0,k_count_y,(j+0.5)*del_y,k_count_z,k*del_z,f_adj,&
                        t_spread,spread,counter-1.0,del_x,del_t,c,pi,imag)*& 
                        k_num_z_exception
                        Hy(k_pl_start_H,j,k)=Hy(k_pl_start_H,j,k)+del_t/(mu_0*del_x)*&
                        Inc(WEz,0.0,k_count_y,j*del_y,k_count_z,(k+0.5)*del_z,f_adj,&
                        t_spread,spread,counter-1.0,del_x,del_t,c,pi,imag)*& 
                        k_num_y_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !!Y FACES!!
            if (pbc_x+pbc_z==2) then
                !$acc parallel loop collapse(2) present(Hz,Hx)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh
                    do i=xlow,xhigh
                        Hz(i,k_pl_start_H,k)=Hz(i,k_pl_start_H,k)+del_t/(mu_0*del_y)*&
                        Inc(WEx,0.0,k_count_x,(i+0.5)*del_x,k_count_z,k*del_z,f_adj,&
                        t_spread,spread,counter-1.0,del_y,del_t,c,pi,imag)*& 
                        k_num_z_exception
                        Hx(i,k_pl_start_H,k)=Hx(i,k_pl_start_H,k)+del_t/(mu_0*del_y)*&
                        Inc(WEz,0.0,k_count_x,i*del_x,k_count_z,(k+0.5)*del_z,f_adj,&
                        t_spread,spread,counter-1.0,del_y,del_t,c,pi,imag)*& 
                        k_num_x_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !!Z FACES!!
            if (pbc_y+pbc_x==2) then
                !$acc parallel loop collapse(2) present(Hx,Hy)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh
                    do i=xlow,xhigh
                        Hx(i,j,k_pl_start_H)=Hx(i,j,k_pl_start_H)+del_t/(mu_0*del_z)*&
                        Inc(WEy,0.0,k_count_y,(j+0.5)*del_y,k_count_x,i*del_x,f_adj,&
                        t_spread,spread,counter-1.0,del_z,del_t,c,pi,imag)*& 
                        k_num_x_exception
                        Hy(i,j,k_pl_start_H)=Hy(i,j,k_pl_start_H)+del_t/(mu_0*del_z)*&
                        Inc(WEx,0.0,k_count_y,j*del_y,k_count_x,(i+0.5)*del_x,f_adj,&
                        t_spread,spread,counter-1.0,del_z,del_t,c,pi,imag)*& 
                        k_num_y_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

        end if
#endif

#ifndef use_kmax_version
        if (plane_wave_amp>0) then
            !!Y FACES!!
            if (pbc_y==0) then
                !$acc parallel loop collapse(2) present(Hz)
                !$omp do collapse(2) schedule(static)    
                do k=zlow,zhigh 
                    do i=xlow,xhigh-1+pbc_x
                        Hz(i,ylow-1,k)=Hz(i,ylow-1,k)-del_t/(mu_0*del_y)*Inc(WEx,i+0.5,ylow+0.0,k+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*ylow_wall
                        Hz(i,yhigh,k)=Hz(i,yhigh,k)+del_t/(mu_0*del_y)*Inc(WEx,i+0.5,yhigh+0.0,k+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx)
                !$omp do collapse(2) schedule(static)
                do k=zlow, zhigh-1+pbc_z
                    do i=xlow,xhigh
                        Hx(i,ylow-1,k)=Hx(i,ylow-1,k)+del_t/(mu_0*del_y)*Inc(WEz,i+0.0,ylow+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*ylow_wall
                        Hx(i,yhigh,k)=Hx(i,yhigh,k)-del_t/(mu_0*del_y)*Inc(WEz,i+0.0,yhigh+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!Z FACES!!
            if (pbc_z==0) then
                !$acc parallel loop collapse(2) present(Hy)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh
                    do i=xlow,xhigh-1+pbc_x
                        Hy(i,j,zlow-1)=Hy(i,j,zlow-1)+del_t/(mu_0*del_z)*Inc(WEx,i+0.5,j+0.0,zlow+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zlow_wall
                        Hy(i,j,zhigh)=Hy(i,j,zhigh)-del_t/(mu_0*del_z)*Inc(WEx,i+0.5,j+0.0,zhigh+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh-1+pbc_y
                    do i=xlow,xhigh
                        Hx(i,j,zlow-1)=Hx(i,j,zlow-1)-del_t/(mu_0*del_z)*Inc(WEy,i+0.0,j+0.5,zlow+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zlow_wall
                        Hx(i,j,zhigh)=Hx(i,j,zhigh)+del_t/(mu_0*del_z)*Inc(WEy,i+0.0,j+0.5,zhigh+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!X FACES!!
            if (pbc_x==0) then
                !$acc parallel loop collapse(2) present(Hz)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh
                    do j=ylow,yhigh-1+pbc_y
                        Hz(xlow-1,j,k)=Hz(xlow-1,j,k)+del_t/(mu_0*del_x)*Inc(WEy,xlow+0.0,j+0.5,k+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xlow_wall
                        Hz(xhigh,j,k)=Hz(xhigh,j,k)-del_t/(mu_0*del_x)*Inc(WEy,xhigh+0.0,j+0.5,k+0.0,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hy)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do j=ylow,yhigh
                        Hy(xlow-1,j,k)=Hy(xlow-1,j,k)-del_t/(mu_0*del_x)*Inc(WEz,xlow+0.0,j+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xlow_wall
                        Hy(xhigh,j,k)=Hy(xhigh,j,k)+del_t/(mu_0*del_x)*Inc(WEz,xhigh+0.0,j+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
        end if

        if ((is_mirror==1) .and. (plane_wave_amp>0)) then
            !!Y FACES!!
            if (pbc_y==0) then
                !$acc parallel loop collapse(2) present(Hz)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh
                    do i=xlow,xhigh-1+pbc_x
                        Hz(i,ylow-1,k)=Hz(i,ylow-1,k)-del_t/(mu_0*del_y)*Inc(WEx_mirror,i+0.5,ylow+0.0,k+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*ylow_wall
                        Hz(i,yhigh,k)=Hz(i,yhigh,k)+del_t/(mu_0*del_y)*Inc(WEx_mirror,i+0.5,yhigh+0.0,k+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx)
                !$omp do collapse(2) schedule(static)
                do k=zlow, zhigh-1+pbc_z
                    do i=xlow,xhigh
                        Hx(i,ylow-1,k)=Hx(i,ylow-1,k)+del_t/(mu_0*del_y)*Inc(WEz_mirror,i+0.0,ylow+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*ylow_wall
                        Hx(i,yhigh,k)=Hx(i,yhigh,k)-del_t/(mu_0*del_y)*Inc(WEz_mirror,i+0.0,yhigh+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!Z FACES!!
            if (pbc_z==0) then
                !$acc parallel loop collapse(2) present(Hy)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh
                    do i=xlow,xhigh-1+pbc_x
                        Hy(i,j,zlow-1)=Hy(i,j,zlow-1)+del_t/(mu_0*del_z)*Inc(WEx_mirror,i+0.5,j+0.0,zlow+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zlow_wall
                        Hy(i,j,zhigh)=Hy(i,j,zhigh)-del_t/(mu_0*del_z)*Inc(WEx_mirror,i+0.5,j+0.0,zhigh+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh-1+pbc_y
                    do i=xlow,xhigh
                        Hx(i,j,zlow-1)=Hx(i,j,zlow-1)-del_t/(mu_0*del_z)*Inc(WEy_mirror,i+0.0,j+0.5,zlow+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zlow_wall
                        Hx(i,j,zhigh)=Hx(i,j,zhigh)+del_t/(mu_0*del_z)*Inc(WEy_mirror,i+0.0,j+0.5,zhigh+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!X FACES!!
            if (pbc_x==0) then
                !$acc parallel loop collapse(2) present(Hz)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh
                    do j=ylow,yhigh-1+pbc_y
                        Hz(xlow-1,j,k)=Hz(xlow-1,j,k)+del_t/(mu_0*del_x)*Inc(WEy_mirror,xlow+0.0,j+0.5,k+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xlow_wall
                        Hz(xhigh,j,k)=Hz(xhigh,j,k)-del_t/(mu_0*del_x)*Inc(WEy_mirror,xhigh+0.0,j+0.5,k+0.0,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hy)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do j=ylow,yhigh
                        Hy(xlow-1,j,k)=Hy(xlow-1,j,k)-del_t/(mu_0*del_x)*Inc(WEz_mirror,xlow+0.0,j+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xlow_wall
                        Hy(xhigh,j,k)=Hy(xhigh,j,k)+del_t/(mu_0*del_x)*Inc(WEz_mirror,xhigh+0.0,j+0.0,k+0.5,&
                        t_spread,spread,counter-1.0,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
        end if
#endif

        !!INTERNAL PORTS IF APPLICABLE!! 
        !4 sections - I-H integrals/determination, I-V updates, Sources, gridded H attachement

        !$omp single
        !$acc serial present(Current, Hy, Hz, Hx, ports,p_loc,p_loc_d,gridded_feed_all)
        do i=1,num_ports

            if (ports(i,5,1)==0) then !then basic lumped ports are used
            
                !these are the starting and stopping index locations
                p_loc(1)=ports(i,3,1)
                p_loc(2)=ports(i,3,2)
                p_loc(3)=ports(i,3,3)
                p_loc_d(1)=ports(i,4,1)+p_loc(1)-1
                p_loc_d(2)=ports(i,4,2)+p_loc(2)-1
                p_loc_d(3)=ports(i,4,3)+p_loc(3)-1
                
                !reset the current source from fdtd before summing and determining here
                Current(i,1)=0.0+0.0*imag

                !next will pick version to use based on x,y,z
                select case (int(ports(i,1,1)))

                case (0)
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            Current(i,1)=Current(i,1)+&
                            del_y*(Hy(ii,jj,p_loc(3)-1)-Hy(ii,jj,p_loc_d(3)))
                        end do
                        do kk=p_loc(3),p_loc_d(3)
                            Current(i,1)=Current(i,1)+&
                            del_z*(Hz(ii,p_loc(2),kk)-Hz(ii,p_loc_d(2)-1,kk))
                        end do
                    end do
                    Current(i,1)=Current(i,1)*real(ports(i,4,2)*ports(i,4,3))/real(ports(i,4,1))

                case (1)
                    do jj=p_loc(2),p_loc_d(2)
                        do ii=p_loc(1),p_loc_d(1)
                            Current(i,1)=Current(i,1)+&
                            del_x*(Hx(ii,jj,p_loc(3))-Hx(ii,jj,p_loc_d(3)-1))
                        end do
                        do kk=p_loc(3),p_loc_d(3)
                            Current(i,1)=Current(i,1)+&
                            del_z*(Hz(p_loc(1)-1,jj,kk)-Hz(p_loc_d(1),jj,kk))
                        end do
                    end do
                    Current(i,1)=Current(i,1)*real(ports(i,4,1)*ports(i,4,3))/real(ports(i,4,2))

                case (2)
                    do kk=p_loc(3),p_loc_d(3)
                        do ii=p_loc(1),p_loc_d(1)
                            Current(i,1)=Current(i,1)+&
                            del_x*(Hx(ii,p_loc(2)-1,kk)-Hx(ii,p_loc_d(2),kk))
                        end do
                        do jj=p_loc(2),p_loc_d(2)
                            Current(i,1)=Current(i,1)+&
                            del_y*(Hy(p_loc(1),jj,kk)-Hy(p_loc_d(1)-1,jj,kk))
                        end do
                    end do
                    Current(i,1)=Current(i,1)*real(ports(i,4,1)*ports(i,4,2))/real(ports(i,4,3))

                end select
            end if

            if (ports(i,5,1)==1) then !then gridded lumped ports are used

                !center cell location
                p_loc(1)=ports(i,3,1)
                p_loc(2)=ports(i,3,2)
                p_loc(3)=ports(i,3,3)
                !grid file size in cells
                p_loc_d(1)=ports(i,4,1)
                p_loc_d(2)=ports(i,4,2)
                
                !reset the current source from fdtd before summing and determining here
                Current(i,1)=0.0+0.0*imag

                !similar, factor before looping each port and time step
                jj_real=0.0
                
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports(i,1,2)))

                !math is H=I*coeff but we will use a mode decomposition to ensure we only pull the TEM portion and ignore higher order and evanescent modes
                !some required switching the sign of the current for cross product reasonings
                !x_axis, y_axis with negative z_axis current direction so
                !x: y,z with negative x directed current - positive x directed port aligns - so case 0 needs a negative correction
                !y: x,z with positive y directed current - negative y directed port aligns - so case 3 needs a negative correction
                !z: x,y with negative z directed current - positive z directed port aligns - so case 4 needs a negative correction

                case(0)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hy(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,3,jj,kk) &
                            +Hz(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do
                    current(i,1)=current(i,1)*-1.0
                    
                case(1)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hy(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,3,jj,kk) &
                            +Hz(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do

                case(2)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hx(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,3,jj,kk) &
                            +Hz(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do

                case(3)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hx(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)*gridded_feed_all(i,3,jj,kk) &
                            +Hz(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do
                    current(i,1)=current(i,1)*-1.0

                case(4)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)*gridded_feed_all(i,3,jj,kk) &
                            +Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do
                    current(i,1)=current(i,1)*-1.0

                case(5)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            current(i,1)=current(i,1)+Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))*gridded_feed_all(i,3,jj,kk) &
                            +Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))*gridded_feed_all(i,4,jj,kk)
                        end do
                    end do

                end select

                !normalization for mode decomposition
                do jj=1, p_loc_d(1)+1
                    do kk=1, p_loc_d(2)+1
                        jj_real=jj_real+(gridded_feed_all(i,3,jj+1,kk))**2.0+&
                        (gridded_feed_all(i,4,jj,kk))**2.0
                    end do
                end do
                current(i,1)=current(i,1)/jj_real

            end if
        end do
        !$acc end serial
        !$omp end single

        !$acc parallel loop collapse(2) present(Current, Voltage, ports) private(step_vc_out)
        !$omp do collapse(2) schedule(static) 
        do i=1,num_ports
            do j=2, port_array_size-1

                !pick version to use based on x,y,z
                select case (int(ports(i,1,1)))
                case (0) 
                    step_vc_out=del_x
                case (1)
                    step_vc_out=del_y
                case (2)
                    step_vc_out=del_z
                end select

                !update the rest of the 1D transmission line
                Current(i,j)=Current(i,j)-&
                (1.0/ports(i,2,2))*(del_t/step_vc_out)*(Voltage(i,j+1)-Voltage(i,j))

            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !$omp single
        !$acc serial present(Current,ports)
        do i=1,num_ports
        
                !pick version to use based on x,y,z
                select case (int(ports(i,1,1)))
                case (0) 
                    step_vc_out=del_x
                case (1)
                    step_vc_out=del_y
                case (2)
                    step_vc_out=del_z
                end select

#ifndef use_kmax_version
                if (i==excitation_port_number) then
                    Current(i,2)=Current(i,2)-(1.0/ports(i,2,2))*(del_t/step_vc_out)*&
                    E_or_V_1D(antenna_amp,counter+0.0,0.0*step_vc_out,spread,t_spread,pulse_type,del_t,c)
                end if
#endif
#ifdef use_kmax_version                
                if (i==excitation_port_number) then
                    Current(i,2)=Current(i,2)-(1.0/ports(i,2,2))*(del_t/step_vc_out)*&
                    E_or_V_1D(antenna_amp,counter+0.0,0.0*step_vc_out,t_spread,del_t,c,spread,pi,f_adj,imag)
                end if
#endif
        end do
        !$acc end serial
        !$omp end single

        !Gridded H attachment to grid if needed
        !$omp single
        !$acc serial present(ports,p_loc,p_loc_d,Hx,Hy,Hz, gridded_feed_all)
        do i=1, num_ports
            if (ports(i,5,1)==1) then

                !center cell location
                p_loc(1)=ports(i,3,1)
                p_loc(2)=ports(i,3,2)
                p_loc(3)=ports(i,3,3)
                !grid file size in cells
                p_loc_d(1)=ports(i,4,1)
                p_loc_d(2)=ports(i,4,2)
                
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports(i,1,2)))

                !purposefully zero out the H field behind the port cross section
                !it doesn't contribute anywhere above port since E if enforced, but it keeps backscattering to zero

                case(0)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hy(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hz(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(1)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hy(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(2)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(3)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(4)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(5)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all(i,4,jj,kk)/=0) then
                                Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do                    
                    
                end select

            end if
        end do
        !$acc end serial
        !$omp end single

#ifdef use_spice_version
        !Gridded H attachment to grid for spice if needed
        !$omp single
        !$acc serial present(ports_spice,p_loc,p_loc_d,Hx,Hy,Hz,gridded_feed_all_spice)
        do i=1, num_spice_ports
            if (ports_spice(i,5,1)==1) then

                !center cell location
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                !grid file size in cells
                p_loc_d(1)=ports_spice(i,3,1)
                p_loc_d(2)=ports_spice(i,3,2)
                
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports_spice(i,1,2)))

                !purposefully zero out the H field behind the port cross section
                !it doesn't contribute anywhere above port since E if enforced, but it keeps backscattering to zero

                case(0)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hy(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hz(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(1)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hy(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(2)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(3)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hz(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do
                    
                case(4)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=0.0+0.0*imag
                            end if
                        end do
                    end do

                case(5)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,3,jj,kk)/=0) then
                                Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)=0.0+0.0*imag
                            end if
                            if (gridded_feed_all_spice(i,4,jj,kk)/=0) then
                                Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)=0.0+0.0*imag
                            end if
                        end do
                    end do                    
                    
                end select

            end if
        end do
        !$acc end serial
        !$omp end single
#endif

        !!CPML HX -Y and +Y DIRECTION!!
        !$acc parallel loop collapse(2) present(Hx, Ez, psi_Hxy_1, psi_Hxy_2, bh_y_1, ch_y_1, bh_y_2, ch_y_2)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do i = 1,x_size-1
                !!CPML HX -Y DIRECTION!!
                do j = 1,nyPML_1-1
                    psi_Hxy_1(i,j,k) = bh_y_1(j)*psi_Hxy_1(i,j,k) &
                    + ch_y_1(j) *(Ez(i,j,k) - Ez(i,j+1,k))/del_y
                    Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_1(i,j,k)
                end do
                
                !!CPML HX +Y DIRECTION!!
                jj = nyPML_2-1
                do j = y_size+1-nyPML_2,y_size-1
                    psi_Hxy_2(i,jj,k) = bh_y_2(jj)*psi_Hxy_2(i,jj,k) &
                        + ch_y_2(jj) *(Ez(i,j,k) - Ez(i,j+1,k))/del_y
                    Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_2(i,jj,k)
                    jj = jj-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML HX -Z and +Z DIRECTION!!
        !$acc parallel loop collapse(2) present(Hx, Ey, psi_Hxz_1, psi_Hxz_2, bh_z_1, ch_z_1, bh_z_2, ch_z_2)
        !$omp do collapse(2) schedule(static)
        do j = 1,y_size-1
            do i = 1,x_size-1 
                !!CPML HX -Z DIRECTION!!
                do k = 1,nzPML_1-1
                    psi_Hxz_1(i,j,k) = bh_z_1(k)*psi_Hxz_1(i,j,k) &
                    + ch_z_1(k) *(Ey(i,j,k+1) - Ey(i,j,k))/del_z
                    Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_1(i,j,k)
                end do
                
                !!CPML HX +Z DIRECTION!!
                kk = nzPML_2-1
                do k = z_size+1-nzPML_2,z_size-1
                    psi_Hxz_2(i,j,kk) = bh_z_2(kk)*psi_Hxz_2(i,j,kk) &
                    + ch_z_2(kk) *(Ey(i,j,k+1) - Ey(i,j,k))/del_z
                    Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_2(i,j,kk)
                    kk = kk-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML HY -X and +X DIRECTION!!
        !$acc parallel loop collapse(2) present(Hy, Ez, psi_Hyx_1, psi_Hyx_2, bh_x_1, ch_x_1, bh_x_2, ch_x_2)
        !$omp do collapse(2) schedule(static)        
        do k = 1,z_size-1
            do j = 1,y_size-1
                !!CPML HY -X DIRECTION!!
                do i = 1,nxPML_1-1
                    psi_Hyx_1(i,j,k) = bh_x_1(i)*psi_Hyx_1(i,j,k) &
                    + ch_x_1(i)*(Ez(i+1,j,k) - Ez(i,j,k))/del_x
                    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_1(i,j,k)
                end do
                
                !!CPML HY +X DIRECTION!!
                ii = nxPML_2-1
                do i = x_size+1-nxPML_2,x_size-1
                    psi_Hyx_2(ii,j,k) = bh_x_2(ii)*psi_Hyx_2(ii,j,k) &
                    + ch_x_2(ii)*(Ez(i+1,j,k) - Ez(i,j,k))/del_x
                    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_2(ii,j,k)
                    ii = ii-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML HY -Z and +Z DIRECTION!!
        !$acc parallel loop collapse(2) present(Hy, Ex, psi_Hyz_1, psi_Hyz_2, bh_z_1, ch_z_1, bh_z_2, ch_z_2)
        !$omp do collapse(2) schedule(static)        
        do j = 1,y_size-1
            do i = 1,x_size-1
                !!CPML HY -Z DIRECTION!!
                do k = 1,nzPML_1-1
                    psi_Hyz_1(i,j,k) = bh_z_1(k)*psi_Hyz_1(i,j,k) &
                    + ch_z_1(k)*(Ex(i,j,k) - Ex(i,j,k+1))/del_z
                    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_1(i,j,k)
                end do
                
                !!CPML HY +Z DIRECTION!!
                kk = nzPML_2-1
                do k = z_size+1-nzPML_2,z_size-1
                    psi_Hyz_2(i,j,kk) = bh_z_2(kk)*psi_Hyz_2(i,j,kk) &
                    + ch_z_2(kk)*(Ex(i,j,k) - Ex(i,j,k+1))/del_z
                    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_2(i,j,kk)
                    kk = kk-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML HZ -X and +X DIRECTION!!
        !$acc parallel loop collapse(2) present(Hz, Ey, psi_Hzx_1, psi_Hzx_2, bh_x_1, ch_x_1, bh_x_2, ch_x_2)
        !$omp do collapse(2) schedule(static)        
        do k = 1,z_size-1
            do j = 1,y_size-1
                !!CPML HZ -X DIRECTION!!
                do i = 1,nxPML_1-1
                    psi_Hzx_1(i,j,k) = bh_x_1(i)*psi_Hzx_1(i,j,k) &
                    + ch_x_1(i) *(Ey(i,j,k) - Ey(i+1,j,k))/del_x
                    Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_1(i,j,k)
                end do
                
                !!CPML HZ +X DIRECTION!!
                ii = nxPML_2-1
                do i = x_size+1-nxPML_2,x_size-1
                    psi_Hzx_2(ii,j,k) = bh_x_2(ii)*psi_Hzx_2(ii,j,k) &
                    + ch_x_2(ii) *(Ey(i,j,k) - Ey(i+1,j,k))/del_x
                    Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_2(ii,j,k)
                    ii = ii-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML HZ -Y and +Y DIRECTION!!
        !$acc parallel loop collapse(2) present(Hz, Ex, psi_Hzy_1, psi_Hzy_2, bh_y_1, ch_y_1, bh_y_2, ch_y_2)
        !$omp do collapse(2) schedule(static)        
        do k = 1,z_size-1
            do i = 1,x_size-1
                !!CPML HZ -Y DIRECTION!!
                do j = 1,nyPML_1-1
                    psi_Hzy_1(i,j,k) = bh_y_1(j)*psi_Hzy_1(i,j,k) &
                    + ch_y_1(j)*(Ex(i,j+1,k) - Ex(i,j,k))/del_y
                    Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_1(i,j,k)
                end do
                
                !!CPML HZ +Y DIRECTION!!
                jj = nyPML_2-1
                do j = y_size+1-nyPML_2,y_size-1
                    psi_Hzy_2(i,jj,k) = bh_y_2(jj)*psi_Hzy_2(i,jj,k) &
                    + ch_y_2(jj)*(Ex(i,j+1,k) - Ex(i,j,k))/del_y
                    Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_2(i,jj,k)
                    jj = jj-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!PBC H-FIELDS!!
        if (pbc_x==1) then
            !$acc parallel loop collapse(2) present(Hy, Hz)
            !$omp do collapse(2) schedule(static)
            do k = 1, z_size-1
                do j = 1, y_size-1
#ifndef use_kmax_version
                    Hz(x_size-1,j,k) = Hz(1,j,k)
                    Hy(x_size-1,j,k) = Hy(1,j,k)
#endif
#ifdef use_kmax_version
                    Hz(x_size-1,j,k) = Hz(1,j,k)*CEXP(-1.0*imag*k_count_x*(x_size-2)*del_x)
                    Hy(x_size-1,j,k) = Hy(1,j,k)*CEXP(-1.0*imag*k_count_x*(x_size-2)*del_x)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Y boundaries
        if (pbc_y==1) then
            !$acc parallel loop collapse(2) present(Hz, Hx)
            !$omp do collapse(2) schedule(static)
            do k = 1, z_size-1
                do i = 1, x_size-1
#ifndef use_kmax_version
                    Hz(i,y_size-1,k) = Hz(i,1,k)
                    Hx(i,y_size-1,k) = Hx(i,1,k)
#endif
#ifdef use_kmax_version
                    Hz(i,y_size-1,k) = Hz(i,1,k)*CEXP(-1.0*imag*k_count_y*(y_size-2)*del_y)
                    Hx(i,y_size-1,k) = Hx(i,1,k)*CEXP(-1.0*imag*k_count_y*(y_size-2)*del_y)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Z boundaries
        if (pbc_z==1) then
            !$acc parallel loop collapse(2) present(Hy, Hx)
            !$omp do collapse(2) schedule(static)
            do j = 1, y_size-1
                do i = 1, x_size-1
#ifndef use_kmax_version
                    Hy(i,j,z_size-1) = Hy(i,j,1)
                    Hx(i,j,z_size-1) = Hx(i,j,1)
#endif
#ifdef use_kmax_version
                    Hy(i,j,z_size-1) = Hy(i,j,1)*CEXP(-1.0*imag*k_count_z*(z_size-2)*del_z)
                    Hx(i,j,z_size-1) = Hx(i,j,1)*CEXP(-1.0*imag*k_count_z*(z_size-2)*del_z)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if

#ifdef use_spice_version
        !!!SPICE ROUTINE IF USED!!!
        !normally nothing should go after pbc update, but this doesn't interfere with H udpates in any way - it just pulls info from H fields
        !and it compiles slightly faster this way for some reason. Maybe because it sends data out of the simulation.
        !First will need to get current, like before, from H fields
        !Then we send the current to spice and return a voltage in a time synchronized manner

        !$acc serial present(Spice_Current,Hx,Hy,Hz,p_loc,p_loc_d,ports_spice,gridded_feed_all_spice)
        !$omp single
        do i=1, num_spice_ports

            if (ports_spice(i,5,1)==0) then !then basic lumped ports are used

                !these are the starting and stopping index locations
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                p_loc_d(1)=ports_spice(i,3,1)+p_loc(1)-1
                p_loc_d(2)=ports_spice(i,3,2)+p_loc(2)-1
                p_loc_d(3)=ports_spice(i,3,3)+p_loc(3)-1
                !reset the current source from fdtd before summing and determining here
                Spice_Current(i)=0.0+0.0*imag
                !next will pick version to use based on x,y,z
                select case (int(ports_spice(i,1,1)))

                case (0)
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_y*(Hy(ii,jj,p_loc(3)-1)-Hy(ii,jj,p_loc_d(3)))
                        end do
                        do kk=p_loc(3),p_loc_d(3)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_z*(Hz(ii,p_loc(2),kk)-Hz(ii,p_loc_d(2)-1,kk))
                        end do
                    end do
                    Spice_Current(i)=Spice_Current(i)*real(ports_spice(i,3,2)*ports_spice(i,3,3))/real(ports_spice(i,3,1))     
                case (1)
                    do jj=p_loc(2),p_loc_d(2)
                        do ii=p_loc(1),p_loc_d(1)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_x*(Hx(ii,jj,p_loc(3))-Hx(ii,jj,p_loc_d(3)-1))
                        end do
                        do kk=p_loc(3),p_loc_d(3)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_z*(Hz(p_loc(1)-1,jj,kk)-Hz(p_loc_d(1),jj,kk))
                        end do
                    end do
                    Spice_Current(i)=Spice_Current(i)*real(ports_spice(i,3,1)*ports_spice(i,3,3))/real(ports_spice(i,3,2))
                case (2)
                    do kk=p_loc(3),p_loc_d(3)
                        do ii=p_loc(1),p_loc_d(1)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_x*(Hx(ii,p_loc(2)-1,kk)-Hx(ii,p_loc_d(2),kk))
                        end do
                        do jj=p_loc(2),p_loc_d(2)
                            Spice_Current(i)=Spice_Current(i)+&
                            del_y*(Hy(p_loc(1),jj,kk)-Hy(p_loc_d(1)-1,jj,kk))
                        end do
                    end do
                    Spice_Current(i)=Spice_Current(i)*real(ports_spice(i,3,1)*ports_spice(i,3,2))/real(ports_spice(i,3,3))

                end select

            end if

            if (ports_spice(i,5,1)==1) then !then gridded lumped ports are used

                !center cell location
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                !grid file size in cells
                p_loc_d(1)=ports_spice(i,3,1)
                p_loc_d(2)=ports_spice(i,3,2)
                
                !reset the current source from fdtd before summing and determining here
                Spice_Current(i)=0.0+0.0*imag

                !similar, reset factor before looping each port and time step
                jj_real=0.0
                
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports_spice(i,1,2)))

                !math is H=I*coeff but we will use a mode decomposition to ensure we only pull the TEM portion and ignore higher order and evanescent modes
                !some required switching the sign of the current for cross product reasonings
                !x_axis, y_axis with negative z_axis current direction so
                !x: y,z with negative x directed current - positive x directed port aligns - so case 0 needs a negative correction
                !y: x,z with positive y directed current - negative y directed port aligns - so case 3 needs a negative correction
                !z: x,y with negative z directed current - positive z directed port aligns - so case 4 needs a negative correction

                case(0)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hy(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hz(p_loc(1)-1,jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do
                    spice_current(i)=spice_current(i)*-1.0
                    
                case(1)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hy(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hz(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do

                case(2)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hx(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hz(p_loc(1)+jj-1,p_loc(2)-1,kk+p_loc(3)-1)*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do

                case(3)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hx(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hz(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do
                    spice_current(i)=spice_current(i)*-1.0

                case(4)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)-1)*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do
                    spice_current(i)=spice_current(i)*-1.0

                case(5)
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            !mode decomposition
                            spice_current(i)=spice_current(i)+Hx(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))*gridded_feed_all_spice(i,3,jj,kk) &
                            +Hy(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))*gridded_feed_all_spice(i,4,jj,kk)
                        end do
                    end do

                end select

                !normalization for mode decomposition
                do jj=1, p_loc_d(1)+1
                    do kk=1, p_loc_d(2)+1
                        jj_real=jj_real+(gridded_feed_all_spice(i,3,jj+1,kk))**2.0+&
                        (gridded_feed_all_spice(i,4,jj,kk))**2.0
                    end do
                end do
                spice_current(i)=spice_current(i)/jj_real

            end if
        end do
        !$acc end serial
        !$omp end single
        
        !$acc update self(Spice_Current) 
        !$omp single
        do i=1, num_spice_ports
            !write(*,*) 'Time:', circuit%time, ' Altering at this time.'
            !write(*,*) Spice_Current(i)

            ! Construct the PWL command string using an internal write
            command_str = 'alter @'//trim(names_of_spice_ports_currents(i))//'[PWL] = ['

    #ifdef use_kmax_version
            !changed for kmax to deal with real and imaginary parts of the wave
            !doesn't really matter which one is odd and even - they just need to come in two pairs
            !additionally, I must be consistent when recombining them in voltage section
            if (mod(i, 2) == 0) then
                !even i values
                write(pair_str, '(E15.8, 1x, E15.8)') del_t*counter, aimag(Spice_Current(i))
            else
                !odd i values
                write(pair_str, '(E15.8, 1x, E15.8)') del_t*counter, real(Spice_Current(i))
            end if
    #endif
    #ifndef use_kmax_version
            write(pair_str, '(E15.8, 1x, E15.8)') del_t*counter, Spice_Current(i)
    #endif

            command_str = trim(command_str) // ' ' // trim(pair_str)
            command_str = trim(command_str) // ']'
            !write(*,*) trim(command_str) // char(0)
            
            ! --- Execute the command and step the simulation ---
            call command(trim(command_str) // char(0))
        end do
        !$omp end single

        !$omp single
        if (num_spice_ports>0) then
            call circuit%step()

            !write(*, '(A6, F20.15, A9, E15.8)') 'Time:', circuit%getTime(), ' V(out):', circuit%nodes%values(1)%voltage
            ! Check for time synchronization errors
            if (abs(circuit%getTime() - circuit%time) > 0.01) then
                error_cnt = error_cnt + 1
                write(*,*) 'Error at time:', circuit%getTime()
                write(*,*) 'Expected time:', circuit%time
            end if

            ! Update time for the next iteration
            circuit%time = circuit%time + circuit%dt

            !write(*,*) circuit%time, counter*del_t
        end if
        !$omp end single

        !$omp single
        do i=1, num_spice_ports
            Spice_Voltage(i)=circuit%nodes%values(i)%voltage-ports_spice(i,4,1)
            !write(*,*) Spice_Voltage(i)
        end do
        !$omp end single
        !$acc update device(Spice_Voltage)
        !!!END SPICE ROUTINE IF USED!!!
#endif

        !!UPDATE EX!!
        !$acc parallel loop collapse(3) present(Ex, Hz, Hy, gax, gbx, den_ey, den_ez)
        !$omp do collapse(2) schedule(static)
        do k = 2,z_size-1
            do j = 2,y_size-1
                do i = 1,x_size-1
                    Ex(i,j,k) = gax(i,j,k) * Ex(i,j,k) + gbx(i,j,k) * &
                    ((Hz(i,j,k) - Hz(i,j-1,k))*den_ey(j)  + &
                    (Hy(i,j,k-1) - Hy(i,j,k))*den_ez(k))
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!UPDATE EY!!
        !$acc parallel loop collapse(3) present(Ey, Hz, Hx, gay, gby, den_ex, den_ez)
        !$omp do collapse(2) schedule(static)
        do k = 2,z_size-1
            do j = 1,y_size-1
                do i = 2,x_size-1 
                    Ey(i,j,k) = gay(i,j,k) * Ey(i,j,k) + gby(i,j,k) * &
                    ((Hz(i-1,j,k) - Hz(i,j,k))*den_ex(i) + &
                    (Hx(i,j,k) - Hx(i,j,k-1))*den_ez(k))
                end do 
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!UPDATE EZ!!
        !$acc parallel loop collapse(3) present(Ez, Hy, Hx, gaz, gbz, den_ex, den_ey)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do j = 2,y_size-1
                do i = 2,x_size-1
                    Ez(i,j,k) = Gaz(i,j,k) * Ez(i,j,k) + Gbz(i,j,k)  &
                    * ((Hy(i,j,k) - Hy(i-1,j,k))*den_ex(i) + &
                    (Hx(i,j-1,k) - Hx(i,j,k))*den_ey(j))
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!SHEET IMPEDANCES!! - special sub-cell update
        !!X-SHEETS!!
        !$acc parallel loop collapse(3) present(Ex_special, Hz, Hy, sheet_sig_x_x, sheet_ep_x_x, den_ey, den_ez, x_sheet_list)
        !$omp do collapse(2) schedule(static)        
        do ii = 1, count_unique_sheets_x
            do k = 2, z_size-1
                do j = 2, y_size-1
                    i = x_sheet_list(ii)
                    
                    !x sheets if they exist
                    if ((sheet_sig_x_x(i,j,k) > 0) .or. (sheet_ep_x_x(i,j,k) > 1)) then
                        Ex_special(i,j,k) = ((1.0-sheet_sig_x_x(i,j,k)*del_t/(2.0*ep_0*sheet_ep_x_x(i,j,k))) / &
                        (1.0+sheet_sig_x_x(i,j,k)*del_t/(2.0*ep_0*sheet_ep_x_x(i,j,k))))*Ex_special(i,j,k) + &
                        ((del_t/(ep_0*sheet_ep_x_x(i,j,k)))/ &
                        (1.0+sheet_sig_x_x(i,j,k)*del_t/(2.0*ep_0*sheet_ep_x_x(i,j,k)))) * &
                        ((Hz(i,j,k) - Hz(i,j-1,k))*den_ey(j) + &
                        (Hy(i,j,k-1) - Hy(i,j,k))*den_ez(k))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!Y-SHEETS!!
        !$acc parallel loop collapse(3) present(Ey_special, Hz, Hx, sheet_sig_y_y, sheet_ep_y_y, den_ex, den_ez, y_sheet_list)
        !$omp do collapse(2) schedule(static)
        do jj = 1, count_unique_sheets_y
            do k = 2, z_size-1
                do i = 2, x_size-1
                    j = y_sheet_list(jj)
                    
                    !y sheets if they exist
                    if ((sheet_sig_y_y(i,j,k) > 0) .or. (sheet_ep_y_y(i,j,k) > 1)) then
                        Ey_special(i,j,k) = ((1.0-sheet_sig_y_y(i,j,k)*del_t/(2.0*ep_0*sheet_ep_y_y(i,j,k))) / &
                        (1.0+sheet_sig_y_y(i,j,k)*del_t/(2.0*ep_0*sheet_ep_y_y(i,j,k))))*Ey_special(i,j,k) + &
                        ((del_t/(ep_0*sheet_ep_y_y(i,j,k)))/ &
                        (1.0+sheet_sig_y_y(i,j,k)*del_t/(2.0*ep_0*sheet_ep_y_y(i,j,k)))) * &
                        ((Hz(i-1,j,k) - Hz(i,j,k))*den_ex(i) + &
                        (Hx(i,j,k) - Hx(i,j,k-1))*den_ez(k))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!Z-SHEETS!!
        !$acc parallel loop collapse(3) present(Ez_special, Hy, Hx, sheet_sig_z_z, sheet_ep_z_z, den_ex, den_ey, z_sheet_list)
        !$omp do collapse(2) schedule(static)        
        do kk = 1, count_unique_sheets_z
            do j = 2, y_size-1
                do i = 2, x_size-1
                    k = z_sheet_list(kk)
                    
                    !z sheets if they exist
                    if ((sheet_sig_z_z(i,j,k) > 0) .or. (sheet_ep_z_z(i,j,k) > 1)) then
                        Ez_special(i,j,k) = ((1.0-sheet_sig_z_z(i,j,k)*del_t/(2.0*ep_0*sheet_ep_z_z(i,j,k))) / &
                        (1.0+sheet_sig_z_z(i,j,k)*del_t/(2.0*ep_0*sheet_ep_z_z(i,j,k))))*Ez_special(i,j,k) + &
                        ((del_t/(ep_0*sheet_ep_z_z(i,j,k)))/ &
                        (1.0+sheet_sig_z_z(i,j,k)*del_t/(2.0*ep_0*sheet_ep_z_z(i,j,k)))) * &
                        ((Hy(i,j,k) - Hy(i-1,j,k))*den_ex(i) + &
                        (Hx(i,j-1,k) - Hx(i,j,k))*den_ey(j))
                    end if
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!DISPERSIVE MATERIAL UPDATE IF NEEDED!!
        !first filter if poles are present (any dispersive media), then filter on type.
        if (num_poles>0) then

            if (is_plasma==1) then

                !$acc parallel loop collapse(3) present(J_plasma_ax,J_plasma_bx,J_source_x,Ex,Ex_oldt,gbx) &
                !$acc private(temp_plasma_x, temp_plasma_x_c, plasma_counter)
                !$omp do collapse(2) schedule(static)
                do k=plasma_min_xfields_zpos, plasma_max_xfields_zpos
                    do j=plasma_min_xfields_ypos, plasma_max_xfields_ypos
                        do i=plasma_min_xfields_xpos, plasma_max_xfields_xpos
                            temp_plasma_x_c=0.0+0.0*imag
                            do plasma_counter=1, num_poles
                                temp_plasma_x_c=temp_plasma_x_c+(1.0+J_plasma_ax(plasma_counter,i,j,k))*J_source_x(plasma_counter,i,j,k)
                            end do
                            Ex(i,j,k)=Ex(i,j,k)-gbx(i,j,k)*temp_plasma_x_c/2.0
                            do plasma_counter=1, num_poles
                                J_source_x(plasma_counter,i,j,k)=J_plasma_ax(plasma_counter,i,j,k)*J_source_x(plasma_counter,i,j,k)+&
                                J_plasma_bx(plasma_counter,i,j,k)*(Ex(i,j,k)+Ex_oldt(i,j,k))
                            end do
                        end do
                    end do
                end do
                !$acc end parallel loop
                !$omp end do

                !$acc parallel loop collapse(3) present(J_plasma_ay,J_plasma_by,J_source_y,Ey,Ey_oldt,gby) &
                !$acc private(temp_plasma_y, temp_plasma_y_c, plasma_counter)
                !$omp do collapse(2) schedule(static)
                do k=plasma_min_yfields_zpos, plasma_max_yfields_zpos
                    do j=plasma_min_yfields_ypos, plasma_max_yfields_ypos
                        do i=plasma_min_yfields_xpos, plasma_max_yfields_xpos
                            temp_plasma_y_c=0.0+0.0*imag
                            do plasma_counter=1, num_poles
                                temp_plasma_y_c=temp_plasma_y_c+(1.0+J_plasma_ay(plasma_counter,i,j,k))*J_source_y(plasma_counter,i,j,k)
                            end do
                            Ey(i,j,k)=Ey(i,j,k)-gby(i,j,k)*temp_plasma_y_c/2.0
                            do plasma_counter=1, num_poles
                                J_source_y(plasma_counter,i,j,k)=J_plasma_ay(plasma_counter,i,j,k)*J_source_y(plasma_counter,i,j,k)+&
                                J_plasma_by(plasma_counter,i,j,k)*(Ey(i,j,k)+Ey_oldt(i,j,k))
                            end do
                        end do
                    end do
                end do
                !$omp end do
                !$acc end parallel loop

                !$acc parallel loop collapse(3) present(J_plasma_az,J_plasma_bz,J_source_z,Ez,Ez_oldt,gbz) &
                !$acc private(temp_plasma_z, temp_plasma_z_c, plasma_counter)
                !$omp do collapse(2) schedule(static)
                do k=plasma_min_zfields_zpos, plasma_max_zfields_zpos
                    do j=plasma_min_zfields_ypos, plasma_max_zfields_ypos
                        do i=plasma_min_zfields_xpos, plasma_max_zfields_xpos
                            temp_plasma_z_c=0.0+0.0*imag
                            do plasma_counter=1, num_poles
                                temp_plasma_z_c=temp_plasma_z_c+(1.0+J_plasma_az(plasma_counter,i,j,k))*J_source_z(plasma_counter,i,j,k)
                            end do
                            Ez(i,j,k)=Ez(i,j,k)-gbz(i,j,k)*temp_plasma_z_c/2.0
                            do plasma_counter=1, num_poles
                                J_source_z(plasma_counter,i,j,k)=J_plasma_az(plasma_counter,i,j,k)*J_source_z(plasma_counter,i,j,k)+&
                                J_plasma_bz(plasma_counter,i,j,k)*(Ez(i,j,k)+Ez_oldt(i,j,k))
                            end do
                        end do
                    end do
                end do
                !$omp end do
                !$acc end parallel loop

            end if

        end if
        !!END DISPERSIVE MEDIA UPDATES!!

        !!INJECT E PLANE WAVE SOURCES FOR TF/SF FORMULATION (if not kmax) INTO THE GRID!!
#ifdef use_kmax_version
        !if TM mode then add H fields to E fields - no longer TF/SF anymore
        !1/2 values still used so that TE and TM get launched from the same physical plane - though this is uncessary
        !differs from main fdtd program in terms of starting locations for both E/H
        if ((mode_type==1) .and. (plane_wave_amp>0)) then

            !!X FACES!!
            if (pbc_y+pbc_z==2) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=1,z_size-1
                    do j=2,y_size-1 
                        Ez(k_pl_start_E,j,k)=Ez(k_pl_start_E,j,k)+del_t/(ep_0*del_x)*&
                        Inc(WHy,-0.5,k_count_y,j*del_y,k_count_z,(k+0.5)*del_z,f_adj,&
                        t_spread,spread,counter-0.5,del_x,del_t,c,pi,imag)*& 
                        k_num_z_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do k=2,z_size-1
                    do j=1,y_size-1 
                        Ey(k_pl_start_E,j,k)=Ey(k_pl_start_E,j,k)+del_t/(ep_0*del_x)*&
                        Inc(WHz,-0.5,k_count_y,(j+0.5)*del_y,k_count_z,k*del_z,f_adj,&
                        t_spread,spread,counter-0.5,del_x,del_t,c,pi,imag)*& 
                        k_num_y_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !!Y FACES!!
            if (pbc_x+pbc_z==2) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=1,z_size-1
                    do i=2,x_size-1 
                        Ez(i,k_pl_start_E,k)=Ez(i,k_pl_start_E,k)+del_t/(ep_0*del_y)*&
                        Inc(WHx,-0.5,k_count_x,i*del_x,k_count_z,(k+0.5)*del_z,f_adj,&
                        t_spread,spread,counter-0.5,del_y,del_t,c,pi,imag)*& 
                        k_num_z_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do k=2,z_size-1
                    do i=1,x_size-1 
                        Ex(i,k_pl_start_E,k)=Ex(i,k_pl_start_E,k)+del_t/(ep_0*del_y)*&
                        Inc(WHz,-0.5,k_count_x,(i+0.5)*del_x,k_count_z,k*del_z,f_adj,&
                        t_spread,spread,counter-0.5,del_y,del_t,c,pi,imag)*& 
                        k_num_x_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !!Z FACES!!
            if (pbc_y+pbc_x==2) then
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do j=2,y_size-1
                    do i=1,x_size-1 
                        Ex(i,j,k_pl_start_E)=Ex(i,j,k_pl_start_E)+del_t/(ep_0*del_z)*&
                        Inc(WHy,-0.5,k_count_y,j*del_y,k_count_x,(i+0.5)*del_x,f_adj,&
                        t_spread,spread,counter-0.5,del_z,del_t,c,pi,imag)*& 
                        k_num_x_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do j=1,y_size-1
                    do i=2,x_size-1 
                        Ey(i,j,k_pl_start_E)=Ey(i,j,k_pl_start_E)+del_t/(ep_0*del_z)*&
                        Inc(WHx,-0.5,k_count_y,(j+0.5)*del_y,k_count_x,i*del_x,f_adj,&
                        t_spread,spread,counter-0.5,del_z,del_t,c,pi,imag)*& 
                        k_num_y_exception
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

        end if
#endif

#ifndef use_kmax_version
        if (plane_wave_amp>0) then
            !!Y FACES!!
            if (pbc_y==0) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do i=xlow+pbc_x,xhigh
                        Ez(i,ylow,k)=Ez(i,ylow,k)+del_t/(ep_0*del_y)*Inc(WHx,i+0.0,ylow-0.5,k+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*ylow_wall
                        Ez(i,yhigh,k)=Ez(i,yhigh,k)-del_t/(ep_0*del_y)*Inc(WHx,i+0.0,yhigh+0.5,k+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do k=zlow+pbc_z, zhigh
                    do i=xlow,xhigh-1+pbc_x
                        Ex(i,ylow,k)=Ex(i,ylow,k)-del_t/(ep_0*del_y)*Inc(WHz,i+0.5,ylow-0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*ylow_wall
                        Ex(i,yhigh,k)=Ex(i,yhigh,k)+del_t/(ep_0*del_y)*Inc(WHz,i+0.5,yhigh+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!Z FACES!!
            if (pbc_z==0) then
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh-1+pbc_y
                    do i=xlow+pbc_x,xhigh
                        Ey(i,j,zlow)=Ey(i,j,zlow)-del_t/(ep_0*del_z)*Inc(WHx,i+0.0,j+0.5,zlow-0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zlow_wall
                        Ey(i,j,zhigh)=Ey(i,j,zhigh)+del_t/(ep_0*del_z)*Inc(WHx,i+0.0,j+0.5,zhigh+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do j=ylow+pbc_y,yhigh
                    do i=xlow,xhigh-1+pbc_x
                        Ex(i,j,zlow)=Ex(i,j,zlow)+del_t/(ep_0*del_z)*Inc(WHy,i+0.5,j+0.0,zlow-0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zlow_wall
                        Ex(i,j,zhigh)=Ex(i,j,zhigh)-del_t/(ep_0*del_z)*Inc(WHy,i+0.5,j+0.0,zhigh+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!X FACES!!
            if (pbc_x==0) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do j=ylow+pbc_y,yhigh
                        Ez(xlow,j,k)=Ez(xlow,j,k)-del_t/(ep_0*del_x)*Inc(WHy,xlow-0.5,j+0.0,k+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xlow_wall
                        Ez(xhigh,j,k)=Ez(xhigh,j,k)+del_t/(ep_0*del_x)*Inc(WHy,xhigh+0.5,j+0.0,k+0.5,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do k=zlow+pbc_z,zhigh
                    do j=ylow,yhigh-1+pbc_y
                        Ey(xlow,j,k)=Ey(xlow,j,k)+del_t/(ep_0*del_x)*Inc(WHz,xlow-0.5,j+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xlow_wall
                        Ey(xhigh,j,k)=Ey(xhigh,j,k)-del_t/(ep_0*del_x)*Inc(WHz,xhigh+0.5,j+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta,phi,del_x,del_y,del_z,x_delay,y_delay,z_delay,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
        end if

        if ((is_mirror==1) .and. (plane_wave_amp>0)) then
            !!Y FACES!!
            if (pbc_y==0) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do i=xlow+pbc_x,xhigh
                        Ez(i,ylow,k)=Ez(i,ylow,k)+del_t/(ep_0*del_y)*Inc(WHx_mirror,i+0.0,ylow-0.5,k+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*ylow_wall
                        Ez(i,yhigh,k)=Ez(i,yhigh,k)-del_t/(ep_0*del_y)*Inc(WHx_mirror,i+0.0,yhigh+0.5,k+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do k=zlow+pbc_z, zhigh
                    do i=xlow,xhigh-1+pbc_x
                        Ex(i,ylow,k)=Ex(i,ylow,k)-del_t/(ep_0*del_y)*Inc(WHz_mirror,i+0.5,ylow-0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*ylow_wall
                        Ex(i,yhigh,k)=Ex(i,yhigh,k)+del_t/(ep_0*del_y)*Inc(WHz_mirror,i+0.5,yhigh+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*yhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!Z FACES!!
            if (pbc_z==0) then
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do j=ylow,yhigh-1+pbc_y
                    do i=xlow+pbc_x,xhigh
                        Ey(i,j,zlow)=Ey(i,j,zlow)-del_t/(ep_0*del_z)*Inc(WHx_mirror,i+0.0,j+0.5,zlow-0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zlow_wall
                        Ey(i,j,zhigh)=Ey(i,j,zhigh)+del_t/(ep_0*del_z)*Inc(WHx_mirror,i+0.0,j+0.5,zhigh+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ex)
                !$omp do collapse(2) schedule(static)
                do j=ylow+pbc_y,yhigh
                    do i=xlow,xhigh-1+pbc_x
                        Ex(i,j,zlow)=Ex(i,j,zlow)+del_t/(ep_0*del_z)*Inc(WHy_mirror,i+0.5,j+0.0,zlow-0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zlow_wall
                        Ex(i,j,zhigh)=Ex(i,j,zhigh)-del_t/(ep_0*del_z)*Inc(WHy_mirror,i+0.5,j+0.0,zhigh+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*zhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !!X FACES!!
            if (pbc_x==0) then
                !$acc parallel loop collapse(2) present (Ez)
                !$omp do collapse(2) schedule(static)
                do k=zlow,zhigh-1+pbc_z
                    do j=ylow+pbc_y,yhigh
                        Ez(xlow,j,k)=Ez(xlow,j,k)-del_t/(ep_0*del_x)*Inc(WHy_mirror,xlow-0.5,j+0.0,k+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xlow_wall
                        Ez(xhigh,j,k)=Ez(xhigh,j,k)+del_t/(ep_0*del_x)*Inc(WHy_mirror,xhigh+0.5,j+0.0,k+0.5,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present (Ey)
                !$omp do collapse(2) schedule(static)
                do k=zlow+pbc_z,zhigh
                    do j=ylow,yhigh-1+pbc_y
                        Ey(xlow,j,k)=Ey(xlow,j,k)+del_t/(ep_0*del_x)*Inc(WHz_mirror,xlow-0.5,j+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xlow_wall
                        Ey(xhigh,j,k)=Ey(xhigh,j,k)-del_t/(ep_0*del_x)*Inc(WHz_mirror,xhigh+0.5,j+0.5,k+0.0,&
                        t_spread,spread,counter-0.5,theta_mirror,phi_mirror,del_x,del_y,del_z,x_delay_mirror,y_delay_mirror,z_delay_mirror,pulse_type,del_t,c)*xhigh_wall
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
        end if
#endif

        !!INTERNAL PORTS IF APPLICABLE!!
        !3 sections - I-V updates, sources, attachment to E-H grid

        !$omp do collapse(2) schedule(static)
        !$acc parallel loop collapse(2) present (Current, Voltage, ports) private (step_vc_out)
        do i=1,num_ports
            do j=2,port_array_size-1

                !pick version to use based on x,y,z
                select case (int(ports(i,1,1)))
                case (0) 
                    step_vc_out=del_x
                case (1)
                    step_vc_out=del_y
                case (2)
                    step_vc_out=del_z
                end select

                Voltage(i,j)=Voltage(i,j)-&
                (1.0/ports(i,2,3))*(del_t/step_vc_out)*(Current(i,j)-Current(i,j-1))

            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !$omp single
        !$acc serial present(Voltage, ports)
        do i=1,num_ports

            !pick version to use based on x,y,z
            select case (int(ports(i,1,1)))
            case (0) 
                step_vc_out=del_x
            case (1)
                step_vc_out=del_y
            case (2)
                step_vc_out=del_z
            end select

#ifndef use_kmax_version
            if (i==excitation_port_number) then
                Voltage(i,2)=Voltage(i,2)+(1.0/sqrt(ports(i,2,3)*ports(i,2,2)))*(del_t/step_vc_out)*&
                E_or_V_1D(antenna_amp,counter+0.5,0.5*step_vc_out,spread,t_spread,pulse_type,del_t,c)
            end if
#endif
#ifdef use_kmax_version
            if (i==excitation_port_number) then
                Voltage(i,2)=Voltage(i,2)+(1.0/sqrt(ports(i,2,3)*ports(i,2,2)))*(del_t/step_vc_out)*&
                E_or_V_1D(antenna_amp,counter+0.5,0.5*step_vc_out,t_spread,del_t,c,spread,pi,f_adj,imag)
            end if
#endif
        end do
        !$acc end serial
        !$omp end single

        !$omp single
        !$acc serial present(Voltage, ports, Ey, Ez, Ex, p_loc,p_loc_d,gridded_feed_all)
        do i=1,num_ports

            if (ports(i,5,1)==0) then

                !these are the starting and stopping index locations
                p_loc(1)=ports(i,3,1)
                p_loc(2)=ports(i,3,2)
                p_loc(3)=ports(i,3,3)
                p_loc_d(1)=ports(i,4,1)+p_loc(1)-1
                p_loc_d(2)=ports(i,4,2)+p_loc(2)-1
                p_loc_d(3)=ports(i,4,3)+p_loc(3)-1

                select case (int(ports(i,1,1)))

                case (0)

                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Ex(ii,jj,kk)=1.0*(Voltage(i,2))/((ports(i,4,1))*(del_x))
                            end do
                        end do
                    end do
                
                case (1)
                    
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Ey(ii,jj,kk)=1.0*(Voltage(i,2))/((ports(i,4,2))*(del_y))
                            end do
                        end do
                    end do
                
                case (2)

                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
                                Ez(ii,jj,kk)=1.0*(Voltage(i,2))/((ports(i,4,3))*(del_z))
                            end do
                        end do
                    end do

                end select
            end if
            
            if (ports(i,5,1)==1) then

                !bottom left cell of the mask
                p_loc(1)=ports(i,3,1)
                p_loc(2)=ports(i,3,2)
                p_loc(3)=ports(i,3,3)
                !grid file size in cells
                p_loc_d(1)=ports(i,4,1)
                p_loc_d(2)=ports(i,4,2)
                
                !Uses formula for E=V*coeff to evaluate - no normal step_size needed here like for H
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports(i,1,2)))

                case(0,1)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,1,jj,kk)/=0) then
                                Ey(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Voltage(i,2)*gridded_feed_all(i,1,jj,kk)
                            end if
                            if (gridded_feed_all(i,2,jj,kk)/=0) then
                                Ez(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Voltage(i,2)*gridded_feed_all(i,2,jj,kk)
                            end if
                        end do
                    end do

                case(2,3)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,1,jj,kk)/=0) then
                                Ex(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Voltage(i,2)*gridded_feed_all(i,1,jj,kk)
                            end if
                            if (gridded_feed_all(i,2,jj,kk)/=0) then
                                Ez(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Voltage(i,2)*gridded_feed_all(i,2,jj,kk)
                            end if
                        end do
                    end do

                case(4,5)

                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all(i,1,jj,kk)/=0) then
                                Ex(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Voltage(i,2)*gridded_feed_all(i,1,jj,kk)
                            end if
                            if (gridded_feed_all(i,2,jj,kk)/=0) then
                                Ey(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Voltage(i,2)*gridded_feed_all(i,2,jj,kk)
                            end if
                        end do
                    end do

                end select
            end if

        end do
        !$acc end serial 
        !$omp end single

        !!!SPICE ATTACHMENT TO GRID - ALL OTHER SPICE DONE BETWEEN H AND E SECTIONS + H section if gridded!!!
#ifdef use_spice_version
        !$acc serial present(Spice_Voltage, ports_spice, Ey, Ez, Ex, p_loc,p_loc_d,gridded_feed_all_spice)
        !$omp single
        do i=1,num_spice_ports
            if (ports_spice(i,5,1)==0) then
                !these are the starting and stopping index locations
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                p_loc_d(1)=ports_spice(i,3,1)+p_loc(1)-1
                p_loc_d(2)=ports_spice(i,3,2)+p_loc(2)-1
                p_loc_d(3)=ports_spice(i,3,3)+p_loc(3)-1

                select case (int(ports_spice(i,1,1)))

                case (0)
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
    #ifdef use_kmax_version
                                if (mod(i, 2) == 0) then
                                    !i is even so imaginary to match current section
                                    !i is even comes after odd so must be added to current value
                                    Ex(ii,jj,kk)=Ex(ii,jj,kk)+imag*((Spice_Voltage(i)))/((ports_spice(i,3,1))*(del_x))
                                else
                                    !i is odd so real to match current section
                                    !is is odd will go first so it's a reset of equals
                                    Ex(ii,jj,kk)=1.0*((Spice_Voltage(i)))/((ports_spice(i,3,1))*(del_x))
                                end if
    #endif
    #ifndef use_kmax_version
                                Ex(ii,jj,kk)=1.0*((Spice_Voltage(i)))/((ports_spice(i,3,1))*(del_x))
    #endif
                            end do
                        end do
                    end do
                    
                case (1)
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
    #ifdef use_kmax_version
                                if (mod(i, 2) == 0) then
                                    !i is even so imaginary to match current section
                                    !i is even comes after odd so must be added to current value
                                    Ey(ii,jj,kk)=Ey(ii,jj,kk)+imag*(Spice_Voltage(i))/((ports_spice(i,3,2))*(del_y))
                                else
                                    !i is odd so real to match current section
                                    !is is odd will go first so it's a reset of equals
                                    Ey(ii,jj,kk)=1.0*((Spice_Voltage(i)))/((ports_spice(i,3,2))*(del_y))
                                end if
    #endif
    #ifndef use_kmax_version
                                Ey(ii,jj,kk)=1.0*((Spice_Voltage(i)))/((ports_spice(i,3,2))*(del_y))
    #endif
                            end do
                        end do
                    end do

                case (2)
                    do ii=p_loc(1),p_loc_d(1)
                        do jj=p_loc(2),p_loc_d(2)
                            do kk=p_loc(3),p_loc_d(3)
    #ifdef use_kmax_version
                                if (mod(i, 2) == 0) then
                                    !i is even so imaginary to match current section
                                    !i is even comes after odd so must be added to current value
                                    Ez(ii,jj,kk)=Ez(ii,jj,kk)+imag*(Spice_Voltage(i))/((ports_spice(i,3,3))*(del_z))
                                else
                                    !i is odd so real to match current section
                                    !is is odd will go first so it's a reset of equals
                                    Ez(ii,jj,kk)=1.0*(Spice_Voltage(i))/((ports_spice(i,3,3))*(del_z))
                                end if
    #endif
    #ifndef use_kmax_version
                                Ez(ii,jj,kk)=1.0*(Spice_Voltage(i))/((ports_spice(i,3,3))*(del_z))
    #endif
                            end do
                        end do
                    end do
                    
                end select
            end if

            if (ports_spice(i,5,1)==1) then
                !bottom left cell of the mask
                p_loc(1)=ports_spice(i,2,1)
                p_loc(2)=ports_spice(i,2,2)
                p_loc(3)=ports_spice(i,2,3)
                !grid file size in cells
                p_loc_d(1)=ports_spice(i,3,1)
                p_loc_d(2)=ports_spice(i,3,2)
                
                !Uses formula for E=V*coeff to evaluate - no normal step_size needed here like for H
                !next will pick version to use based on x,y,z (- or +)
                select case (int(ports_spice(i,1,2)))

                case(0,1)
    #ifdef use_kmax_version
                    if (mod(i, 2) == 0) then
                    !i is even so imaginary to match current section
                    !i is even comes after odd so must be added to current value
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ey(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Ey(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ez(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Ez(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    else
                    !i is odd so real to match current section
                    !is is odd will go first so it's a reset of equals
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ey(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ez(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    end if
    #endif
    #ifndef use_kmax_version
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                Ey(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                            end if
                            if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                Ez(p_loc(1),jj+p_loc(2)-1,kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                            end if
                        end do
                    end do
    #endif
                case(2,3)
    #ifdef use_kmax_version
                    if (mod(i, 2) == 0) then
                    !i is even so imaginary to match current section
                    !i is even comes after odd so must be added to current value
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ex(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Ex(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ez(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Ez(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    else
                    !i is odd so real to match current section
                    !is is odd will go first so it's a reset of equals
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ex(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ez(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    end if
    #endif
    #ifndef use_kmax_version
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                Ex(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                            end if
                            if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                Ez(p_loc(1)+jj-1,p_loc(2),kk+p_loc(3)-1)=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                            end if
                        end do
                    end do
    #endif
                case(4,5)
    #ifdef use_kmax_version
                    if (mod(i, 2) == 0) then
                    !i is even so imaginary to match current section
                    !i is even comes after odd so must be added to current value
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ex(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Ex(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ey(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Ey(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3)) + &
                                    imag*Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    else
                    !i is odd so real to match current section
                    !is is odd will go first so it's a reset of equals
                        do jj=1, p_loc_d(1)+1
                            do kk=1, p_loc_d(2)+1
                                if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                    Ex(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                                end if
                                if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                    Ey(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                                end if
                            end do
                        end do
                    end if
    #endif
    #ifndef use_kmax_version
                    do jj=1, p_loc_d(1)+1
                        do kk=1, p_loc_d(2)+1
                            if (gridded_feed_all_spice(i,1,jj,kk)/=0) then
                                Ex(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Spice_Voltage(i)*gridded_feed_all_spice(i,1,jj,kk)
                            end if
                            if (gridded_feed_all_spice(i,2,jj,kk)/=0) then
                                Ey(p_loc(1)+jj-1,p_loc(2)+kk-1,p_loc(3))=Spice_Voltage(i)*gridded_feed_all_spice(i,2,jj,kk)
                            end if
                        end do
                    end do
    #endif
                end select
            end if
        end do
        !$omp end single
        !$acc end serial
#endif

        !!CPML EX -Y and +Y DIRECTION!!
        !$acc parallel loop collapse(2) present(Ex, Hz, gbx, psi_Exy_1, psi_Exy_2, be_y_1, ce_y_1, be_y_2, ce_y_2)
        !$omp do collapse(2) schedule(static)
        do k = 2,z_size-1
            do i = 1,x_size-1
                !!CPML EX -Y DIRECTION!!
                do j = 2,nyPML_1
                    psi_Exy_1(i,j,k) = be_y_1(j)*psi_Exy_1(i,j,k) &
                    + ce_y_1(j) *(Hz(i,j,k) - Hz(i,j-1,k))/del_y
                    Ex(i,j,k) = Ex(i,j,k) + gbx(i,j,k)*psi_Exy_1(i,j,k)
                end do
                
                !!CPML EX +Y DIRECTION!!
                jj = nyPML_2
                do j = y_size+1-nyPML_2,y_size-1
                    psi_Exy_2(i,jj,k) = be_y_2(jj)*psi_Exy_2(i,jj,k) &
                    + ce_y_2(jj) *(Hz(i,j,k) - Hz(i,j-1,k))/del_y
                    Ex(i,j,k) = Ex(i,j,k) + gbx(i,j,k)*psi_Exy_2(i,jj,k)
                    jj = jj-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML EX -Z and +Z DIRECTION!!
        !$acc parallel loop collapse(2) present(Ex, Hy, gbx, psi_Exz_1, psi_Exz_2, be_z_1, ce_z_1, be_z_2, ce_z_2)
        !$omp do collapse(2) schedule(static)
        do j = 2,y_size-1
            do i = 1,x_size-1 
                !!CPML EX -Z DIRECTION!!
                do k = 2,nzPML_1
                    psi_Exz_1(i,j,k) = be_z_1(k)*psi_Exz_1(i,j,k) &
                    + ce_z_1(k) *(Hy(i,j,k-1) - Hy(i,j,k))/del_z
                    Ex(i,j,k) = Ex(i,j,k) + gbx(i,j,k)*psi_Exz_1(i,j,k)
                end do
                
                !!CPML EX +Z DIRECTION!!
                kk = nzPML_2
                do k = z_size+1-nzPML_2,z_size-1
                    psi_Exz_2(i,j,kk) = be_z_2(kk)*psi_Exz_2(i,j,kk) &
                    + ce_z_2(kk) *(Hy(i,j,k-1) - Hy(i,j,k))/del_z
                    Ex(i,j,k) = Ex(i,j,k) + gbx(i,j,k)*psi_Exz_2(i,j,kk)
                    kk = kk-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML EY -X and +X DIRECTION!!
        !$acc parallel loop collapse(2) present(Ey, Hz, gby, psi_Eyx_1, psi_Eyx_2, be_x_1, ce_x_1, be_x_2, ce_x_2)
        !$omp do collapse(2) schedule(static) 
        do k = 2,z_size-1
            do j = 1,y_size-1
                !!CPML EY -X DIRECTION!!
                do i = 2,nxPML_1
                    psi_Eyx_1(i,j,k) = be_x_1(i)*psi_Eyx_1(i,j,k) &
                    + ce_x_1(i)*(Hz(i-1,j,k) - Hz(i,j,k))/del_x
                    Ey(i,j,k) = Ey(i,j,k) + gby(i,j,k)*psi_Eyx_1(i,j,k)
                end do
                
                !!CPML EY +X DIRECTION!!
                ii = nxPML_2
                do i = x_size+1-nxPML_2,x_size-1
                    psi_Eyx_2(ii,j,k) = be_x_2(ii)*psi_Eyx_2(ii,j,k) &
                    + ce_x_2(ii)*(Hz(i-1,j,k) - Hz(i,j,k))/del_x
                    Ey(i,j,k) = Ey(i,j,k) + gby(i,j,k)*psi_Eyx_2(ii,j,k)
                    ii = ii-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML EY -Z and +Z DIRECTION!!
        !$acc parallel loop collapse(2) present(Ey, Hx, gby, psi_Eyz_1, psi_Eyz_2, be_z_1, ce_z_1, be_z_2, ce_z_2)
        !$omp do collapse(2) schedule(static)
        do j = 1,y_size-1
            do i = 2,x_size-1
                !!CPML EY -Z DIRECTION!!
                do k = 2,nzPML_1
                    psi_Eyz_1(i,j,k) = be_z_1(k)*psi_Eyz_1(i,j,k) &
                    + ce_z_1(k)*(Hx(i,j,k) - Hx(i,j,k-1))/del_z
                    Ey(i,j,k) = Ey(i,j,k) + gby(i,j,k)*psi_Eyz_1(i,j,k)
                end do
                
                !!CPML EY +Z DIRECTION!!
                kk = nzPML_2
                do k = z_size+1-nzPML_2,z_size-1
                    psi_Eyz_2(i,j,kk) = be_z_2(kk)*psi_Eyz_2(i,j,kk) &
                    + ce_z_2(kk)*(Hx(i,j,k) - Hx(i,j,k-1))/del_z
                    Ey(i,j,k) = Ey(i,j,k) + gby(i,j,k)*psi_Eyz_2(i,j,kk)
                    kk = kk-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML EZ -X and +X DIRECTION!!
        !$acc parallel loop collapse(2) present(Ez, Hy, gbz, psi_Ezx_1, psi_Ezx_2, be_x_1, ce_x_1, be_x_2, ce_x_2)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do j = 2,y_size-1
                !!CPML EZ -X DIRECTION!!
                do i = 2,nxPML_1
                    psi_Ezx_1(i,j,k) = be_x_1(i)*psi_Ezx_1(i,j,k) &
                    + ce_x_1(i) *(Hy(i,j,k) - Hy(i-1,j,k))/del_x
                    Ez(i,j,k) = Ez(i,j,k) + gbz(i,j,k)*psi_Ezx_1(i,j,k)
                end do
                
                !!CPML EZ +X DIRECTION!!
                ii = nxPML_2
                do i = x_size+1-nxPML_2,x_size-1
                    psi_Ezx_2(ii,j,k) = be_x_2(ii)*psi_Ezx_2(ii,j,k) &
                    + ce_x_2(ii) *(Hy(i,j,k) - Hy(i-1,j,k))/del_x
                    Ez(i,j,k) = Ez(i,j,k) + gbz(i,j,k)*psi_Ezx_2(ii,j,k)
                    ii = ii-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do
        
        !!CPML EZ -Y and +Y DIRECTION!!
        !$acc parallel loop collapse(2) present(Ez, Hx, gbz, psi_Ezy_1, psi_Ezy_2, be_y_1, ce_y_1, be_y_2, ce_y_2)
        !$omp do collapse(2) schedule(static)
        do k = 1,z_size-1
            do i = 2,x_size-1
                !!CPML EZ -Y DIRECTION!!
                do j = 2,nyPML_1
                    psi_Ezy_1(i,j,k) = be_y_1(j)*psi_Ezy_1(i,j,k) &
                    + ce_y_1(j)*(Hx(i,j-1,k) - Hx(i,j,k))/del_y
                    Ez(i,j,k) = Ez(i,j,k) + gbz(i,j,k)*psi_Ezy_1(i,j,k)
                end do
                
                !!CPML EZ +Y DIRECTION!!
                jj = nyPML_2
                do j = y_size+1-nyPML_2,y_size-1
                    psi_Ezy_2(i,jj,k) = be_y_2(jj)*psi_Ezy_2(i,jj,k) &
                    + ce_y_2(jj)*(Hx(i,j-1,k) - Hx(i,j,k))/del_y
                    Ez(i,j,k) = Ez(i,j,k) + gbz(i,j,k)*psi_Ezy_2(i,jj,k)
                    jj = jj-1
                end do
            end do
        end do
        !$acc end parallel loop
        !$omp end do

        !!PBC for E-FIELDS!!
        !X boundaries
        if (pbc_x==1) then
            !$acc parallel loop collapse(2) present(Ez, Ez_special, Ey, Ey_special)
            !$omp do collapse(2) schedule(static)
            do k = 1, z_size-1
                do j = 1, y_size-1
#ifndef use_kmax_version
                    Ez(1,j,k) = Ez(x_size-1,j,k)
                    Ez_special(1,j,k) = Ez_special(x_size-1,j,k)
                    Ey(1,j,k) = Ey(x_size-1,j,k)
                    Ey_special(1,j,k) = Ey_special(x_size-1,j,k)
#endif
#ifdef use_kmax_version
                    Ez(1,j,k) = Ez(x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
                    Ez_special(1,j,k) = Ez_special(x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
                    Ey(1,j,k) = Ey(x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
                    Ey_special(1,j,k) = Ey_special(x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Y boundaries
        if (pbc_y==1) then
            !$acc parallel loop collapse(2) present(Ez, Ez_special, Ex, Ex_special)
            !$omp do collapse(2) schedule(static)
            do k = 1, z_size-1
                do i = 1, x_size-1
#ifndef use_kmax_version
                    Ez(i,1,k) = Ez(i,y_size-1,k)
                    Ez_special(i,1,k) = Ez_special(i,y_size-1,k)
                    Ex(i,1,k) = Ex(i,y_size-1,k)
                    Ex_special(i,1,k) = Ex_special(i,y_size-1,k)
#endif
#ifdef use_kmax_version
                    Ez(i,1,k) = Ez(i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
                    Ez_special(i,1,k) = Ez_special(i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
                    Ex(i,1,k) = Ex(i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
                    Ex_special(i,1,k) = Ex_special(i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Z boundaries
        if (pbc_z==1) then
            !$acc parallel loop collapse(2) present(Ey, Ey_special, Ex, Ex_special)
            !$omp do collapse(2) schedule(static)
            do j = 1, y_size-1
                do i = 1, x_size-1
#ifndef use_kmax_version
                    Ey(i,j,1) = Ey(i,j,z_size-1)
                    Ey_special(i,j,1) = Ey_special(i,j,z_size-1)
                    Ex(i,j,1) = Ex(i,j,z_size-1)
                    Ex_special(i,j,1) = Ex_special(i,j,z_size-1)
#endif
#ifdef use_kmax_version
                    Ey(i,j,1) = Ey(i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
                    Ey_special(i,j,1) = Ey_special(i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
                    Ex(i,j,1) = Ex(i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
                    Ex_special(i,j,1) = Ex_special(i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
#endif
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        if ((num_poles>0) .and. (is_plasma==1)) then
            !X boundaries
            if (pbc_x==1) then
                !$acc parallel loop collapse(2) present(J_source_y,J_source_z)
                !$omp do collapse(2) schedule(static)
                do k = 1, z_size-1
                    do j = 1, y_size-1
                        do plasma_counter=1, num_poles
#ifndef use_kmax_version
                            J_source_z(plasma_counter,1,j,k)=J_source_z(plasma_counter,x_size-1,j,k)
                            J_source_y(plasma_counter,1,j,k)=J_source_y(plasma_counter,x_size-1,j,k)
#endif
#ifdef use_kmax_version
                            J_source_z(plasma_counter,1,j,k)=J_source_z(plasma_counter,x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
                            J_source_y(plasma_counter,1,j,k)=J_source_y(plasma_counter,x_size-1,j,k)*CEXP(imag*k_count_x*(x_size-2)*del_x)
#endif
                        end do
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !Y boundaries
            if (pbc_y==1) then
                !$acc parallel loop collapse(2) present(J_source_x,J_source_z)
                !$omp do collapse(2) schedule(static)
                do k = 1, z_size-1
                    do i = 1, x_size-1
                        do plasma_counter=1, num_poles
#ifndef use_kmax_version
                            J_source_z(plasma_counter,i,1,k)=J_source_z(plasma_counter,i,y_size-1,k)
                            J_source_x(plasma_counter,i,1,k)=J_source_x(plasma_counter,i,y_size-1,k)
#endif
#ifdef use_kmax_version
                            J_source_z(plasma_counter,i,1,k)=J_source_z(plasma_counter,i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
                            J_source_x(plasma_counter,i,1,k)=J_source_x(plasma_counter,i,y_size-1,k)*CEXP(imag*k_count_y*(y_size-2)*del_y)
#endif
                        end do
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
            !Z boundaries
            if (pbc_z==1) then
                !$acc parallel loop collapse(2) present(J_source_x,J_source_y)
                !$omp do collapse(2) schedule(static)
                do j = 1, y_size-1
                    do i = 1, x_size-1
                        do plasma_counter=1, num_poles
#ifndef use_kmax_version
                            J_source_y(plasma_counter,i,j,1)=J_source_y(plasma_counter,i,j,z_size-1)
                            J_source_x(plasma_counter,i,j,1)=J_source_x(plasma_counter,i,j,z_size-1)
#endif
#ifdef use_kmax_version
                            J_source_y(plasma_counter,i,j,1)=J_source_y(plasma_counter,i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
                            J_source_x(plasma_counter,i,j,1)=J_source_x(plasma_counter,i,j,z_size-1)*CEXP(imag*k_count_z*(z_size-2)*del_z)
#endif
                        end do
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if
        end if

        !!ADDED OLD TIME UPDATE FOR DISPERSIVE MEDIA!!
        !we need to save the E fields at this time step so we will have the n-1 time step saved always.
        !we could get fancy here and elminate a little waste since plasma uses these at only specific bounds, and each component is a little different.
        !but's it would be a small time savor and this is easier to read.
        if ((num_poles>0) .and. (is_plasma==1)) then
            !$acc parallel loop collapse(3) present(Ex,Ey,Ez,Ex_oldt,Ey_oldt,Ez_oldt)
            !$omp do collapse(2) schedule(static)
            do k=1, z_size-1
                do j=1, y_size-1
                    do i=1, x_size-1
                        Ex_oldt(i,j,k)=Ex(i,j,k)
                        Ey_oldt(i,j,k)=Ey(i,j,k)
                        Ez_oldt(i,j,k)=Ez(i,j,k)
                    end do
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!SEVERAL DATA SETS GET CONTRIBUTIONS AT EACH TIME STEP FOR EXPORTING AT END!! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !video/field arrays for exporting
        !current method is cell centered at the same time - so videos are same as material cells now
        !this leads to issues at first and last time step for H fields but this is okay - both should be zero to small.
        if (video_on==1) then
            if (slice==0) then
                i=slice_location
                !$acc parallel loop collapse(2) present(Ex_video, Ey_video, Ez_video, Ex, Ey, Ez)
                !$omp do collapse(2) schedule(static)
                do k = 1, vid_size2
                    do j = 1, vid_size1
                        !E fields are the base time step
                        Ex_video(j,k,counter) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0
                        Ey_video(j,k,counter) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0
                        Ez_video(j,k,counter) = (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                !$omp do collapse(2) schedule(static)
                do k = 1, vid_size2
                    do j = 1, vid_size1
                        !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                        Hx_video(j,k,counter) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                        Hy_video(j,k,counter) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                        Hz_video(j,k,counter) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                if (counter<time_steps) then
                    !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                    !$omp do collapse(2) schedule(static)
                    do k = 1, vid_size2
                        do j = 1, vid_size1
                            !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                            Hx_video(j,k,counter+1) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                            Hy_video(j,k,counter+1) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                            Hz_video(j,k,counter+1) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                        end do
                    end do
                    !$acc end parallel loop
                    !$omp end do
                end if
            else if (slice==1) then
                j=slice_location
                !$acc parallel loop collapse(2) present(Ex_video, Ey_video, Ez_video, Ex, Ey, Ez)
                !$omp do collapse(2) schedule(static)
                do k = 1, vid_size2
                    do i = 1, vid_size1
                        !E fields are the base time step
                        Ex_video(i,k,counter) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0
                        Ey_video(i,k,counter) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0
                        Ez_video(i,k,counter) = (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                !$omp do collapse(2) schedule(static)
                do k = 1, vid_size2
                    do i = 1, vid_size1
                        !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                        Hx_video(i,k,counter) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                        Hy_video(i,k,counter) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                        Hz_video(i,k,counter) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                if (counter<time_steps) then
                    !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                    !$omp do collapse(2) schedule(static)
                    do k = 1, vid_size2
                        do i = 1, vid_size1
                            !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                            Hx_video(i,k,counter+1) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                            Hy_video(i,k,counter+1) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                            Hz_video(i,k,counter+1) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                        end do
                    end do
                    !$acc end parallel loop
                    !$omp end do
                end if
            else if (slice==2) then
                k=slice_location
                !$acc parallel loop collapse(2) present(Ex_video, Ey_video, Ez_video, Ex, Ey, Ez)
                !$omp do collapse(2) schedule(static)
                do j = 1, vid_size2
                    do i = 1, vid_size1
                        !E fields are the base time step
                        Ex_video(i,j,counter) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0
                        Ey_video(i,j,counter) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0
                        Ez_video(i,j,counter) = (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                !$omp do collapse(2) schedule(static)
                do j = 1, vid_size2
                    do i = 1, vid_size1
                        !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                        Hx_video(i,j,counter) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                        Hy_video(i,j,counter) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                        Hz_video(i,j,counter) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                    end do
                end do
                !$acc end parallel loop
                !$omp end do
                if (counter<time_steps) then
                    !$acc parallel loop collapse(2) present(Hx_video, Hy_video, Hz_video, Hx, Hy, Hz)
                    !$omp do collapse(2) schedule(static)
                    do j = 1, vid_size2
                        do i = 1, vid_size1
                            !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                            Hx_video(i,j,counter+1) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                            Hy_video(i,j,counter+1) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                            Hz_video(i,j,counter+1) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
                        end do
                    end do
                    !$acc end parallel loop
                    !$omp end do
                end if
            end if
        end if

        !!!S parameters if a full pbc version if used!!!
#ifndef use_kmax_version
        !don't need to cell center here - faster and equivalent
        !X normal unit cell
        if (pbc_y+pbc_z==2) then
            !$acc parallel loop collapse(2) present(Ez, Ey, E_reflected, E_transmitted) &
            !$acc reduction(+:E_reflected(counter), E_transmitted(counter))
            !$omp do collapse(2) schedule(static) &
            !$omp reduction(+:E_reflected(counter), E_transmitted(counter))
            do k = ff_zlow, ff_zhigh 
                do j = ff_ylow, ff_yhigh
                    E_reflected(counter) = E_reflected(counter) +  Ez(xlow,j,k)*sin(pol) - &
                    Ey(xlow,j,k)*cos(pol)*cos(phi)
                    E_transmitted(counter) = E_transmitted(counter) + Ez(xhigh,j,k)*sin(pol) - &
                    Ey(xhigh,j,k)*cos(pol)*cos(phi)
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Y normal unit cell
        if (pbc_x+pbc_z==2) then
            !$acc parallel loop collapse(2) present(Ez, Ex, E_reflected, E_transmitted) &
            !$acc reduction(+:E_reflected(counter), E_transmitted(counter))
            !$omp do collapse(2) schedule(static) &
            !$omp reduction(+:E_reflected(counter), E_transmitted(counter))
            do k = ff_zlow, ff_zhigh
                do i = ff_xlow, ff_xhigh
                    E_reflected(counter) = E_reflected(counter) + Ez(i,ylow,k)*sin(pol) + &
                    Ex(i,ylow,k)*cos(pol)*sin(phi)
                    E_transmitted(counter) = E_transmitted(counter) + Ez(i,yhigh,k)*sin(pol) + &
                    Ex(i,yhigh,k)*cos(pol)*sin(phi)
                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if
        !Z normal unit cell
        if (pbc_x+pbc_y==2) then
            !$acc parallel loop collapse(2) present(Ey, Ex, E_reflected, E_transmitted) &
            !$acc reduction(+:E_reflected(counter), E_transmitted(counter))
            !$omp do collapse(2) schedule(static) &
            !$omp reduction(+:E_reflected(counter), E_transmitted(counter))
            do j = ff_ylow, ff_yhigh
                do i = ff_xlow, ff_xhigh
                    E_reflected(counter) = E_reflected(counter) + &
                    Ey(i,j,zlow)*(-1*cos(pol)*cos(phi)-sin(pol)*cos(theta)*sin(phi)) + &
                    Ex(i,j,zlow)*(cos(pol)*sin(phi)-sin(pol)*cos(theta)*cos(phi))
                    E_transmitted(counter) = E_transmitted(counter) + &
                    Ey(i,j,zhigh)*(-1*cos(pol)*cos(phi)-sin(pol)*cos(theta)*sin(phi)) + &
                    Ex(i,j,zhigh)*(cos(pol)*sin(phi)-sin(pol)*cos(theta)*cos(phi))
                end do 
            end do
            !$acc end parallel loop
            !$omp end do
        end if
#endif

#ifdef use_kmax_version
        !don't need to cell center here - faster and equivalent
        !to obtain incident field information correctly (and for dispersion reasons) a clear case is used for reference
        !because of this, the overall sign is not imported in front of TE,TM.
        !What is important is the sign between E1 and E2 (E1+-E2) in the combination relative to TE,TM
        !TE,TM and refl,trans combinations for 'correct' direction designation determined at export.

        if (pbc_y+pbc_z==2) then
            !$acc parallel loop collapse(2) present(Ez, Ey, E_reflected_TE,E_reflected_TM, E_transmitted_TE,E_transmitted_TM) &
            !$acc private(j,k,temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$acc reduction(+:E_reflected_TE(counter), E_transmitted_TE(counter),E_reflected_TM(counter), E_transmitted_TM(counter))
            !$omp do collapse(2) schedule(static) private(temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$omp reduction(+:E_reflected_TE(counter), E_reflected_TM(counter), E_transmitted_TE(counter), E_transmitted_TM(counter))
                
            do k = ff_zlow, ff_zhigh 
                do j = ff_ylow, ff_yhigh

                    temp_refl_TE=(Ez(xlow,j,k))*k_num_y_exception-(Ey(xlow,j,k))*k_num_z_exception
                    temp_refl_TM=(Ez(xlow,j,k))*k_num_z_exception+(Ey(xlow,j,k))*k_num_y_exception
                    temp_trans_TE=(Ez(xhigh,j,k))*k_num_y_exception-(Ey(xhigh,j,k))*k_num_z_exception
                    temp_trans_TM=(Ez(xhigh,j,k))*k_num_z_exception+(Ey(xhigh,j,k))*k_num_y_exception

                    E_reflected_TE(counter)=E_reflected_TE(counter) + temp_refl_TE
                    E_reflected_TM(counter)=E_reflected_TM(counter) + temp_refl_TM
                    E_transmitted_TE(counter)=E_transmitted_TE(counter) + temp_trans_TE
                    E_transmitted_TM(counter)=E_transmitted_TM(counter) + temp_trans_TM 

                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if

        if (pbc_x+pbc_z==2) then
            !$acc parallel loop collapse(2) present(Ez, Ex, E_reflected_TE,E_reflected_TM, E_transmitted_TE,E_transmitted_TM) &
            !$acc private(i,k,temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$acc reduction(+:E_reflected_TE(counter), E_transmitted_TE(counter),E_reflected_TM(counter), E_transmitted_TM(counter))
            !$omp do collapse(2) schedule(static) private(temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$omp reduction(+:E_reflected_TE(counter), E_reflected_TM(counter), E_transmitted_TE(counter), E_transmitted_TM(counter))
                
            do k = ff_zlow, ff_zhigh 
                do i = ff_xlow, ff_xhigh

                    temp_refl_TE=(Ez(i,ylow,k))*k_num_x_exception-(Ex(i,ylow,k))*k_num_z_exception
                    temp_refl_TM=(Ez(i,ylow,k))*k_num_z_exception+(Ex(i,ylow,k))*k_num_x_exception
                    temp_trans_TE=(Ez(i,yhigh,k))*k_num_x_exception-(Ex(i,yhigh,k))*k_num_z_exception
                    temp_trans_TM=(Ez(i,yhigh,k))*k_num_z_exception+(Ex(i,yhigh,k))*k_num_x_exception

                    E_reflected_TE(counter)=E_reflected_TE(counter) + temp_refl_TE
                    E_reflected_TM(counter)=E_reflected_TM(counter) + temp_refl_TM
                    E_transmitted_TE(counter)=E_transmitted_TE(counter) + temp_trans_TE
                    E_transmitted_TM(counter)=E_transmitted_TM(counter) + temp_trans_TM 

                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if

        if (pbc_y+pbc_x==2) then
            !$acc parallel loop collapse(2) present(Ex, Ey, E_reflected_TE,E_reflected_TM, E_transmitted_TE,E_transmitted_TM) &
            !$acc private(j,i,temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$acc reduction(+:E_reflected_TE(counter), E_transmitted_TE(counter),E_reflected_TM(counter), E_transmitted_TM(counter))
            !$omp do collapse(2) schedule(static) private(temp_refl_TE, temp_refl_TM, temp_trans_TE, temp_trans_TM) &
            !$omp reduction(+:E_reflected_TE(counter), E_reflected_TM(counter), E_transmitted_TE(counter), E_transmitted_TM(counter))
                
            do j = ff_ylow, ff_yhigh
                do i = ff_xlow, ff_xhigh 

                    temp_refl_TE=(Ex(i,j,zlow))*k_num_y_exception-(Ey(i,j,zlow))*k_num_x_exception
                    temp_refl_TM=(Ex(i,j,zlow))*k_num_x_exception+(Ey(i,j,zlow))*k_num_y_exception
                    temp_trans_TE=(Ex(i,j,zhigh))*k_num_y_exception-(Ey(i,j,zhigh))*k_num_x_exception
                    temp_trans_TM=(Ex(i,j,zhigh))*k_num_x_exception+(Ey(i,j,zhigh))*k_num_y_exception

                    E_reflected_TE(counter)=E_reflected_TE(counter) + temp_refl_TE
                    E_reflected_TM(counter)=E_reflected_TM(counter) + temp_refl_TM
                    E_transmitted_TE(counter)=E_transmitted_TE(counter) + temp_trans_TE
                    E_transmitted_TM(counter)=E_transmitted_TM(counter) + temp_trans_TM 

                end do
            end do
            !$acc end parallel loop
            !$omp end do
        end if

#endif
        
        !Now far field angles if requested - on the fly time domain version for far field calculations
        !You get broadband frequency info at each angle selected
        !this is not preferred for many angles and single frequencies - there is a better method for that that is forthcoming
        if (num_far_field_angles>0) then

            !first calculate all current sources as needed
            !they are the same for each angle of extraction at each time step
            
            !X faces
            if (pbc_x==0) then
                !$acc parallel loop collapse(2) present(My_xlow, My_xlow_oldt, Mz_xlow, Mz_xlow_oldt, &
                !$acc Jy_xlow, Jy_xlow_oldt, Jz_xlow, Jz_xlow_oldt, My_xhigh, My_xhigh_oldt, Mz_xhigh, &
                !$acc Mz_xhigh_oldt, Jy_xhigh, Jy_xhigh_oldt, Jz_xhigh, Jz_xhigh_oldt, Ez, Ey, Hz, Hy) &
                !$acc private(i)
                !$omp do collapse(2) schedule(static)
                do k = ff_zlow, ff_zhigh
                    do j = ff_ylow, ff_yhigh

                        i = ff_xlow
                        !M=-n hat cross E where n hat faces outward
                        My_xlow_oldt(j,k) = My_xlow(j,k)
                        My_xlow(j,k) = -1.0*(Ez(i,j,k)+Ez(i,j+1,k)+Ez(i+1,j,k)+Ez(i+1,j+1,k))/4.0*xlow_wall
                        Mz_xlow_oldt(j,k) = Mz_xlow(j,k) 
                        Mz_xlow(j,k) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0*xlow_wall
                        !J=n hat cross H where n hat faces outward
                        Jy_xlow_oldt(j,k) = Jy_xlow(j,k)
                        Jy_xlow(j,k) = (Hz(i,j,k)+Hz(i,j,k+1))/2.0*xlow_wall
                        Jz_xlow_oldt(j,k) = Jz_xlow(j,k)
                        Jz_xlow(j,k) = -1.0*(Hy(i,j,k)+Hy(i,j+1,k))/2.0*xlow_wall

                        i = ff_xhigh
                        !M=-n hat cross E where n hat faces outward
                        My_xhigh_oldt(j,k) = My_xhigh(j,k)
                        My_xhigh(j,k) = (Ez(i,j,k)+Ez(i,j+1,k)+Ez(i+1,j,k)+Ez(i+1,j+1,k))/4.0*xhigh_wall
                        Mz_xhigh_oldt(j,k) = Mz_xhigh(j,k)
                        Mz_xhigh(j,k) = -1.0*(Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0*xhigh_wall
                        !J=n hat cross H where n hat faces outward
                        Jy_xhigh_oldt(j,k) = Jy_xhigh(j,k)
                        Jy_xhigh(j,k) = -1.0*(Hz(i,j,k)+Hz(i,j,k+1))/2.0*xhigh_wall
                        Jz_xhigh_oldt(j,k) = Jz_xhigh(j,k)
                        Jz_xhigh(j,k) = (Hy(i,j,k)+Hy(i,j+1,k))/2.0*xhigh_wall

                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !Y faces
            if (pbc_y==0) then
                !$acc parallel loop collapse(2) present(Mx_ylow, Mx_ylow_oldt, Mz_ylow, Mz_ylow_oldt, &
                !$acc Jx_ylow, Jx_ylow_oldt, Jz_ylow, Jz_ylow_oldt, Mx_yhigh, Mx_yhigh_oldt, Mz_yhigh, &
                !$acc Mz_yhigh_oldt, Jx_yhigh, Jx_yhigh_oldt, Jz_yhigh, Jz_yhigh_oldt, Ez, Ex, Hz, Hx) &
                !$acc private(j)
                !$omp do collapse(2) schedule(static) 
                do k = ff_zlow, ff_zhigh
                    do i = ff_xlow, ff_xhigh

                        j = ff_ylow
                        !M=-n hat cross E where n hat faces outward
                        Mx_ylow_oldt(i,k) = Mx_ylow(i,k)
                        Mx_ylow(i,k) = (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0*ylow_wall
                        Mz_ylow_oldt(i,k) = Mz_ylow(i,k)
                        Mz_ylow(i,k) = -1.0*(Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0*ylow_wall
                        !J=n hat cross H where n hat faces outward
                        Jx_ylow_oldt(i,k) = Jx_ylow(i,k)
                        Jx_ylow(i,k) = -1.0*(Hz(i,j,k)+Hz(i,j,k+1))/2.0*ylow_wall
                        Jz_ylow_oldt(i,k) = Jz_ylow(i,k)
                        Jz_ylow(i,k) = (Hx(i,j,k)+Hx(i+1,j,k))/2.0*ylow_wall

                        j = ff_yhigh
                        !M=-n hat cross E where n hat faces outward
                        Mx_yhigh_oldt(i,k) = Mx_yhigh(i,k)
                        Mx_yhigh(i,k) = -1.0*(Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0*yhigh_wall
                        Mz_yhigh_oldt(i,k) = Mz_yhigh(i,k)
                        Mz_yhigh(i,k) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0*yhigh_wall
                        !J=n hat cross H where n hat faces outward
                        Jx_yhigh_oldt(i,k) = Jx_yhigh(i,k)
                        Jx_yhigh(i,k) = (Hz(i,j,k)+Hz(i,j,k+1))/2.0*yhigh_wall
                        Jz_yhigh_oldt(i,k) = Jz_yhigh(i,k)
                        Jz_yhigh(i,k) = -1.0*(Hx(i,j,k)+Hx(i+1,j,k))/2.0*yhigh_wall

                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !Z faces
            if (pbc_z==0) then
                !$acc parallel loop collapse(2) present(Mx_zlow, Mx_zlow_oldt, My_zlow, My_zlow_oldt, Jx_zlow, &
                !$acc Jx_zlow_oldt, Jy_zlow, Jy_zlow_oldt, Mx_zhigh, Mx_zhigh_oldt, My_zhigh, My_zhigh_oldt, &
                !$acc Jx_zhigh, Jx_zhigh_oldt, Jy_zhigh, Jy_zhigh_oldt, Ey, Ex, Hy, Hx) private(k)
                !$omp do collapse(2) schedule(static)
                do j = ff_ylow, ff_yhigh
                    do i = ff_xlow, ff_xhigh

                        k = ff_zlow
                        !M=-n hat cross E where n hat faces outward
                        Mx_zlow_oldt(i,j) = Mx_zlow(i,j)
                        Mx_zlow(i,j) = -1.0*(Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0*zlow_wall
                        My_zlow_oldt(i,j) = My_zlow(i,j)
                        My_zlow(i,j) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0*zlow_wall
                        !J=n hat cross H where n hat faces outward
                        Jx_zlow_oldt(i,j) = Jx_zlow(i,j)
                        Jx_zlow(i,j) = (Hy(i,j,k)+Hy(i,j+1,k))/2.0*zlow_wall
                        Jy_zlow_oldt(i,j) = Jy_zlow(i,j)
                        Jy_zlow(i,j) = -1.0*(Hx(i,j,k)+Hx(i+1,j,k))/2.0*zlow_wall

                        k = ff_zhigh
                        !M=-n hat cross E where n hat faces outward
                        Mx_zhigh_oldt(i,j) = Mx_zhigh(i,j)
                        Mx_zhigh(i,j) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0*zhigh_wall
                        My_zhigh_oldt(i,j) = My_zhigh(i,j)
                        My_zhigh(i,j) = -1.0*(Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0*zhigh_wall
                        !J=n hat cross H where n hat faces outward
                        Jx_zhigh_oldt(i,j) = Jx_zhigh(i,j)
                        Jx_zhigh(i,j) = -1.0*(Hy(i,j,k)+Hy(i,j+1,k))/2.0*zhigh_wall
                        Jy_zhigh_oldt(i,j) = Jy_zhigh(i,j)
                        Jy_zhigh(i,j) = (Hx(i,j,k)+Hx(i+1,j,k))/2.0*zhigh_wall

                    end do
                end do
                !$acc end parallel loop
                !$omp end do
            end if

            !now prep each angle of extraction

#ifdef use_kmax_version
            !$omp single
            !$acc update self(My_xlow, My_xlow_oldt, Mz_xlow, Mz_xlow_oldt, &
                    !$acc            Jy_xlow, Jy_xlow_oldt, Jz_xlow, Jz_xlow_oldt, &
                    !$acc            My_xhigh, My_xhigh_oldt, Mz_xhigh, Mz_xhigh_oldt, &
                    !$acc            Jy_xhigh, Jy_xhigh_oldt, Jz_xhigh, Jz_xhigh_oldt, &
                    !$acc            Mx_ylow, Mx_ylow_oldt, Mz_ylow, Mz_ylow_oldt, &
                    !$acc            Jx_ylow, Jx_ylow_oldt, Jz_ylow, Jz_ylow_oldt, &
                    !$acc            Mx_yhigh, Mx_yhigh_oldt, Mz_yhigh, Mz_yhigh_oldt, &
                    !$acc            Jx_yhigh, Jx_yhigh_oldt, Jz_yhigh, Jz_yhigh_oldt, &
                    !$acc            Mx_zlow, Mx_zlow_oldt, My_zlow, My_zlow_oldt, Jx_zlow, &
                    !$acc            Jx_zlow_oldt, Jy_zlow, Jy_zlow_oldt, &
                    !$acc            Mx_zhigh, Mx_zhigh_oldt, My_zhigh, My_zhigh_oldt, Jx_zhigh, &
                    !$acc            Jx_zhigh_oldt, Jy_zhigh, Jy_zhigh_oldt, &
                    !$acc            Ux, Uy, Uz, Wx, Wy, Wz, far_field_angles)
#endif

            !X faces contribution
            if (pbc_x==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(My_xlow, My_xlow_oldt, Mz_xlow, Mz_xlow_oldt, &
                !$acc                                    Jy_xlow, Jy_xlow_oldt, Jz_xlow, Jz_xlow_oldt, &
                !$acc                                    My_xhigh, My_xhigh_oldt, Mz_xhigh, Mz_xhigh_oldt, &
                !$acc                                    Jy_xhigh, Jy_xhigh_oldt, Jz_xhigh, Jz_xhigh_oldt, &
                !$acc                                    Uy, Uz, Wy, Wz, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, i)
                !$omp do collapse(3) schedule(static)
#endif
                do kk=1,num_far_field_angles
                    do k=ff_zlow,ff_zhigh
                        do j=ff_ylow,ff_yhigh

                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            i=ff_xlow
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xlow(j,k)-My_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                            
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xlow(j,k)-My_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xlow(j,k)-Mz_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xlow(j,k)-Mz_xlow_oldt(j,k))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xlow(j,k)-Jy_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xlow(j,k)-Jy_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xlow(j,k)-Jz_xlow_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xlow(j,k)-Jz_xlow_oldt(j,k))

                            i=ff_xhigh
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xhigh(j,k)-My_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xhigh(j,k)-My_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xhigh(j,k)-Mz_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif                        
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xhigh(j,k)-Mz_xhigh_oldt(j,k))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xhigh(j,k)-Jy_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xhigh(j,k)-Jy_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xhigh(j,k)-Jz_xhigh_oldt(j,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xhigh(j,k)-Jz_xhigh_oldt(j,k))

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

            !Y faces contribution
            if (pbc_y==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(Mx_ylow, Mx_ylow_oldt, Mz_ylow, Mz_ylow_oldt, &
                !$acc                                    Jx_ylow, Jx_ylow_oldt, Jz_ylow, Jz_ylow_oldt, &
                !$acc                                    Mx_yhigh, Mx_yhigh_oldt, Mz_yhigh, Mz_yhigh_oldt, &
                !$acc                                    Jx_yhigh, Jx_yhigh_oldt, Jz_yhigh, Jz_yhigh_oldt, &
                !$acc                                    Ux, Uz, Wx, Wz, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, j)
                !$omp do collapse(3) schedule(static)
#endif
                do kk=1,num_far_field_angles
                    do k=ff_zlow,ff_zhigh
                        do i=ff_xlow,ff_xhigh

                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            j=ff_ylow
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_ylow(i,k)-Mx_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_ylow(i,k)-Mx_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_ylow(i,k)-Mz_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_ylow(i,k)-Mz_ylow_oldt(i,k))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_ylow(i,k)-Jx_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_ylow(i,k)-Jx_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_ylow(i,k)-Jz_ylow_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_ylow(i,k)-Jz_ylow_oldt(i,k))

                            j=ff_yhigh
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_yhigh(i,k)-Mx_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_yhigh(i,k)-Mx_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_yhigh(i,k)-Mz_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_yhigh(i,k)-Mz_yhigh_oldt(i,k))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_yhigh(i,k)-Jx_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_yhigh(i,k)-Jx_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_yhigh(i,k)-Jz_yhigh_oldt(i,k))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_yhigh(i,k)-Jz_yhigh_oldt(i,k))

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

            !Z faces contribution
            if (pbc_z==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(Mx_zlow, Mx_zlow_oldt, My_zlow, My_zlow_oldt, &
                !$acc                                    Jx_zlow, Jx_zlow_oldt, Jy_zlow, Jy_zlow_oldt, &
                !$acc                                    Mx_zhigh, Mx_zhigh_oldt, My_zhigh, My_zhigh_oldt, &
                !$acc                                    Jx_zhigh, Jx_zhigh_oldt, Jy_zhigh, Jy_zhigh_oldt, &
                !$acc                                    Ux, Uy, Wx, Wy, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, k)
                !$omp do collapse(3) schedule(static)
#endif
                do kk=1,num_far_field_angles
                    do j=ff_ylow,ff_yhigh
                        do i=ff_xlow,ff_xhigh
                        
                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            k=ff_zlow
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zlow(i,j)-My_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zlow(i,j)-My_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zlow(i,j)-Mx_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zlow(i,j)-Mx_zlow_oldt(i,j))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zlow(i,j)-Jy_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zlow(i,j)-Jy_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zlow(i,j)-Jx_zlow_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zlow(i,j)-Jx_zlow_oldt(i,j))

                            k=ff_zhigh
                            time_f_var=(r_for_time_relay-((i-ic)*del_x*sin(ang1)*cos(ang2)+(j-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zhigh(i,j)-My_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zhigh(i,j)-My_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zhigh(i,j)-Mx_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zhigh(i,j)-Mx_zhigh_oldt(i,j))

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zhigh(i,j)-Jy_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zhigh(i,j)-Jy_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zhigh(i,j)-Jx_zhigh_oldt(i,j))
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zhigh(i,j)-Jx_zhigh_oldt(i,j))

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

#ifdef use_kmax_version
            !$omp end single
            !$acc update device(Ux,Uy,Uz,Wx,Wy,Wz)
#endif

        end if

        !is mirror is true and there are far field angles 
        !wall variables will take care of first 5 walls if mirrored, so these are the mirrored walls only
        !probably not needed but I have a third filter in case someone added an IGP for a 4 wall pbc condition
        if (pbc_x+pbc_y+pbc_z<2 .and. num_far_field_angles>0 .and. is_mirror==1) then

            !now prep each angle of extraction

#ifdef use_kmax_version
            !$omp single
#endif

            !X faces contribution (mirrored)
            if (pbc_x==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(My_xlow, My_xlow_oldt, Mz_xlow, Mz_xlow_oldt, &
                !$acc                                    Jy_xlow, Jy_xlow_oldt, Jz_xlow, Jz_xlow_oldt, &
                !$acc                                    My_xhigh, My_xhigh_oldt, Mz_xhigh, Mz_xhigh_oldt, &
                !$acc                                    Jy_xhigh, Jy_xhigh_oldt, Jz_xhigh, Jz_xhigh_oldt, &
                !$acc                                    Uy, Uz, Wy, Wz, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, i, i_mirror, j_mirror, k_mirror)
                !$omp do collapse(3) schedule(static)
#endif
                do kk=1,num_far_field_angles
                    do k=ff_zlow,ff_zhigh
                        do j=ff_ylow,ff_yhigh

                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            i=ff_xlow
                            i_mirror = (i) * (1.0 - use_x_mirror) + (x_mirror_offset - i) * use_x_mirror
                            j_mirror = (j) * (1.0 - use_y_mirror) + (y_mirror_offset - j) * use_y_mirror
                            k_mirror = (k) * (1.0 - use_z_mirror) + (z_mirror_offset - k) * use_z_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xlow(j,k)-My_xlow_oldt(j,k))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xlow(j,k)-My_xlow_oldt(j,k))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xlow(j,k)-Mz_xlow_oldt(j,k))*Mz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xlow(j,k)-Mz_xlow_oldt(j,k))*Mz_mirror

                            ii=int(counter+time_f_var-1.0) 
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xlow(j,k)-Jy_xlow_oldt(j,k))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xlow(j,k)-Jy_xlow_oldt(j,k))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xlow(j,k)-Jz_xlow_oldt(j,k))*Jz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xlow(j,k)-Jz_xlow_oldt(j,k))*Jz_mirror

                            i=ff_xhigh
                            i_mirror = (i) * (1.0 - use_x_mirror) + (x_mirror_offset - i) * use_x_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xhigh(j,k)-My_xhigh_oldt(j,k))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(My_xhigh(j,k)-My_xhigh_oldt(j,k))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xhigh(j,k)-Mz_xhigh_oldt(j,k))*Mz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Mz_xhigh(j,k)-Mz_xhigh_oldt(j,k))*Mz_mirror

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xhigh(j,k)-Jy_xhigh_oldt(j,k))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jy_xhigh(j,k)-Jy_xhigh_oldt(j,k))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xhigh(j,k)-Jz_xhigh_oldt(j,k))*Jz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_y*del_z/(4.0*pi*del_t*c))*(Jz_xhigh(j,k)-Jz_xhigh_oldt(j,k))*Jz_mirror

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

            !Y faces contribution (mirrored)
            if (pbc_y==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(Mx_ylow, Mx_ylow_oldt, Mz_ylow, Mz_ylow_oldt, &
                !$acc                                    Jx_ylow, Jx_ylow_oldt, Jz_ylow, Jz_ylow_oldt, &
                !$acc                                    Mx_yhigh, Mx_yhigh_oldt, Mz_yhigh, Mz_yhigh_oldt, &
                !$acc                                    Jx_yhigh, Jx_yhigh_oldt, Jz_yhigh, Jz_yhigh_oldt, &
                !$acc                                    Ux, Uz, Wx, Wz, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, j, i_mirror, j_mirror, k_mirror)
                !$omp do collapse(3) schedule(static)
#endif
                do kk=1,num_far_field_angles
                    do k=ff_zlow,ff_zhigh
                        do i=ff_xlow,ff_xhigh

                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            j=ff_ylow
                            i_mirror = (i) * (1.0 - use_x_mirror) + (x_mirror_offset - i) * use_x_mirror
                            j_mirror = (j) * (1.0 - use_y_mirror) + (y_mirror_offset - j) * use_y_mirror
                            k_mirror = (k) * (1.0 - use_z_mirror) + (z_mirror_offset - k) * use_z_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_ylow(i,k)-Mx_ylow_oldt(i,k))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_ylow(i,k)-Mx_ylow_oldt(i,k))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_ylow(i,k)-Mz_ylow_oldt(i,k))*Mz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_ylow(i,k)-Mz_ylow_oldt(i,k))*Mz_mirror

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_ylow(i,k)-Jx_ylow_oldt(i,k))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_ylow(i,k)-Jx_ylow_oldt(i,k))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_ylow(i,k)-Jz_ylow_oldt(i,k))*Jz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_ylow(i,k)-Jz_ylow_oldt(i,k))*Jz_mirror

                            j=ff_yhigh
                            j_mirror = (j) * (1.0 - use_y_mirror) + (y_mirror_offset - j) * use_y_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_yhigh(i,k)-Mx_yhigh_oldt(i,k))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mx_yhigh(i,k)-Mx_yhigh_oldt(i,k))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii)=Uz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_yhigh(i,k)-Mz_yhigh_oldt(i,k))*Mz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uz(kk,ii+1)=Uz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Mz_yhigh(i,k)-Mz_yhigh_oldt(i,k))*Mz_mirror

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_yhigh(i,k)-Jx_yhigh_oldt(i,k))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jx_yhigh(i,k)-Jx_yhigh_oldt(i,k))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii)=Wz(kk,ii)+(1.0-jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_yhigh(i,k)-Jz_yhigh_oldt(i,k))*Jz_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wz(kk,ii+1)=Wz(kk,ii+1)+(jj_real)*(del_x*del_z/(4.0*pi*del_t*c))*(Jz_yhigh(i,k)-Jz_yhigh_oldt(i,k))*Jz_mirror

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

            !Z faces contribution (mirrored)
                if (pbc_z==0) then
#ifndef use_kmax_version
                !$acc parallel loop collapse(3) present(Mx_zlow, Mx_zlow_oldt, My_zlow, My_zlow_oldt, &
                !$acc                                    Jx_zlow, Jx_zlow_oldt, Jy_zlow, Jy_zlow_oldt, &
                !$acc                                    Mx_zhigh, Mx_zhigh_oldt, My_zhigh, My_zhigh_oldt, &
                !$acc                                    Jx_zhigh, Jx_zhigh_oldt, Jy_zhigh, Jy_zhigh_oldt, &
                !$acc                                    Ux, Uy, Wx, Wy, far_field_angles) &
                !$acc private(ang1, ang2, time_f_var, ii, jj_real, k, i_mirror, j_mirror, k_mirror)
                !$omp do collapse(3) schedule (static)
#endif
                do kk=1,num_far_field_angles
                    do j=ff_ylow,ff_yhigh
                        do i=ff_xlow,ff_xhigh
                        
                            ang1=far_field_angles(kk,1)
                            ang2=far_field_angles(kk,2)

                            !now U and W w/ r=1 - cancels with other r anyway so value doesn't actually matter
                            !n+1 is counter in this notation - important for comparing to taflove

                            k=ff_zlow
                            i_mirror = (i) * (1.0 - use_x_mirror) + (x_mirror_offset - i) * use_x_mirror
                            j_mirror = (j) * (1.0 - use_y_mirror) + (y_mirror_offset - j) * use_y_mirror
                            k_mirror = (k) * (1.0 - use_z_mirror) + (z_mirror_offset - k) * use_z_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zlow(i,j)-My_zlow_oldt(i,j))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zlow(i,j)-My_zlow_oldt(i,j))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zlow(i,j)-Mx_zlow_oldt(i,j))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zlow(i,j)-Mx_zlow_oldt(i,j))*Mx_mirror

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zlow(i,j)-Jy_zlow_oldt(i,j))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zlow(i,j)-Jy_zlow_oldt(i,j))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zlow(i,j)-Jx_zlow_oldt(i,j))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zlow(i,j)-Jx_zlow_oldt(i,j))*Jx_mirror

                            k=ff_zhigh
                            k_mirror = (k) * (1.0 - use_z_mirror) + (z_mirror_offset - k) * use_z_mirror
                            time_f_var=(r_for_time_relay-((i_mirror-ic)*del_x*sin(ang1)*cos(ang2)&
                            +(j_mirror-jc)*del_y*sin(ang1)*sin(ang2)&
                            +(k_mirror-kc)*del_z*cos(ang1)))/(c*del_t)
                            ii=int(counter+0.5+time_f_var-1.0)
                            jj_real=(counter+0.5+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii)=Uy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zhigh(i,j)-My_zhigh_oldt(i,j))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Uy(kk,ii+1)=Uy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(My_zhigh(i,j)-My_zhigh_oldt(i,j))*My_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii)=Ux(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zhigh(i,j)-Mx_zhigh_oldt(i,j))*Mx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Ux(kk,ii+1)=Ux(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Mx_zhigh(i,j)-Mx_zhigh_oldt(i,j))*Mx_mirror

                            ii=int(counter+time_f_var-1.0)
                            jj_real=(counter+time_f_var-1.0)-ii
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii)=Wy(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zhigh(i,j)-Jy_zhigh_oldt(i,j))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wy(kk,ii+1)=Wy(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jy_zhigh(i,j)-Jy_zhigh_oldt(i,j))*Jy_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii)=Wx(kk,ii)+(1.0-jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zhigh(i,j)-Jx_zhigh_oldt(i,j))*Jx_mirror
#ifndef use_kmax_version
                            !$acc atomic update
                            !$omp atomic update
#endif 
                            Wx(kk,ii+1)=Wx(kk,ii+1)+(jj_real)*(del_x*del_y/(4.0*pi*del_t*c))*(Jx_zhigh(i,j)-Jx_zhigh_oldt(i,j))*Jx_mirror

                        end do
                    end do
                end do
#ifndef use_kmax_version
                !$acc end parallel loop
                !$omp end do
#endif
            end if

#ifdef use_kmax_version
            !$omp end single
            !$acc update device(Ux,Uy,Uz,Wx,Wy,Wz)
#endif
        end if


        !Voltage outputs if any relevant
        !$acc parallel loop present(Voltage_out, Voltage)
        !$omp do schedule(static)
        do i = 1, num_ports
            Voltage_out(i,counter) = Voltage(i,2)
        end do
        !$acc end parallel loop
        !$omp end do

        !added for spice
#ifdef use_spice_version
        !Voltage outputs if any relevant
        !$acc parallel loop present(Spice_Voltage_out, Spice_Voltage)
        !$omp do schedule(static)
        do i=1, num_spice_ports
            Spice_Voltage_out(i,counter)=Spice_Voltage(i)
        end do
        !$acc end parallel loop
        !$omp end do
#endif

        !lastly, I will save the fields at the cube corresponding to ic,jc,kc
        !this allow for easy phase centering of far fields when I use this from a clear case
        !so this goes un-needed and unused for a non-clear case - small overhead though
        !I don't currently use H fields from this, but it's a small overhead and we migth want them at some point
        !$acc serial present(Ex,Ey,Ez,Hx,Hy,Hz,Ex_ff_pc,Ey_ff_pc,Ez_ff_pc,Hx_ff_pc,Hy_ff_pc,Hz_ff_pc)
        !$omp single
        if (num_far_field_angles>0 .and. plane_wave_amp>0) then

            !locations
            i=ic
            j=jc
            k=kc

            !cell centered E fields at E field time step
            Ex_ff_pc(counter) = (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/4.0
            Ey_ff_pc(counter) = (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/4.0
            Ez_ff_pc(counter) = (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/4.0

            !H fields are at the base time step + 1/2 so use cell centered value * 0.5
            Hx_ff_pc(counter) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
            Hy_ff_pc(counter) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
            Hz_ff_pc(counter) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
            if (counter<time_steps) then
                !H fields are at the base time step + 1/2 so use cell centered value * 0.5
                Hx_ff_pc(counter+1) = (Hx(i,j,k)+Hx(i+1,j,k))/4.0
                Hy_ff_pc(counter+1) = (Hy(i,j,k)+Hy(i,j+1,k))/4.0
                Hz_ff_pc(counter+1) = (Hz(i,j,k)+Hz(i,j,k+1))/4.0
            end if

        end if
        !$omp end single
        !$acc end serial

        !$omp single
        !!!! writes out what time step we are out so we can track it
        if (mod(counter,100) == 0) then
            write(*, '(I0, " of ", I0, " time steps")') counter, time_steps
        end if
        !$omp end single
        
    end do

    !$acc end data
    !$omp end parallel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! END MAIN FDTD ALGORITHM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_end)
    write(*, '("Finished main FDTD algorithm. Time to complete was ", F0.3, " seconds.")') &
          real(clock_time_end - clock_time_start) / real(clock_rate)
    write(*, '("Average of ", F0.2, " Mcells per second.")') &
          (1.0E-6 * (x_size-1) * (y_size-1) * (z_size-1)) * time_steps / (real(clock_time_end - clock_time_start) / real(clock_rate))
    write(*, '("Post processing some data...")')
    call system_clock(clock_time_start)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! BEGIN POST PROCESSING AND EXPORT OF DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !S parameters final processing if pbc version and incident wave creation/determination
    !first nornal kmax, then kmax
#ifndef use_kmax_version
    !X Normal
    if (pbc_y+pbc_z==2) then
        do counter=1, time_steps
            E_reflected(counter) = E_reflected(counter) / ((y_size-2.0)*(z_size-2.0))
            E_transmitted(counter) = E_transmitted(counter) / ((y_size-2.0)*(z_size-2.0))
        end do
    end if
    !Y Normal
    if (pbc_x+pbc_z==2) then
        do counter=1, time_steps
            E_reflected(counter) = E_reflected(counter) / ((x_size-2.0)*(z_size-2.0))
            E_transmitted(counter) = E_transmitted(counter) / ((x_size-2.0)*(z_size-2.0))
        end do
    end if
    !Z Normal
    if (pbc_y+pbc_x==2) then
        do counter=1, time_steps
            E_reflected(counter) = E_reflected(counter) / ((y_size-2.0)*(x_size-2.0))
            E_transmitted(counter) = E_transmitted(counter) / ((y_size-2.0)*(x_size-2.0))
        end do
    end if

    !first plane waves excitations if used
    !this is currently only used in an effective aperature area calculation, I think.
    !Clear cases are used to generate equivalents for S-parameters and far fields
    if (plane_wave_amp>0) then
        do counter=1, time_steps
            incident(counter) = E_or_V_1D(plane_wave_amp,counter+0.0,0.0*del_x,&
                                                spread,t_spread,pulse_type,del_t,c)
        end do
    end if
    !then antenna excitation if used
    if (antenna_amp>0) then
        !make a clear grid so we need current and voltage equations with sources added at the right locations
        !first determine if x,y,or z so we know which step size to use
        if (ports(excitation_port_number,1,1)==0) then
            step_vc_out=del_x
        end if
        if (ports(excitation_port_number,1,1)==1) then
            step_vc_out=del_y
        end if
        if (ports(excitation_port_number,1,1)==2) then
            step_vc_out=del_z
        end if
        !start time loop - this sim time scales with time steps, not grid size directly
        do counter=1, time_steps
            !current grid
            do j=1, port_array_size*2-1
                C_inc(j)=C_inc(j)-&
                (1.0/ports(excitation_port_number,2,2))*(del_t/step_vc_out)*(V_inc(j+1)-V_inc(j))
            end do
            !current source in TF/SF
            C_inc(port_array_size)=C_inc(port_array_size)-&
            (1.0/ports(excitation_port_number,2,2))*(del_t/step_vc_out)*&
            E_or_V_1D(antenna_amp,counter+0.0,0.0*step_vc_out,spread,t_spread,pulse_type,del_t,c)
            !voltage grid
            do j=1, port_array_size*2-1
                V_inc(j)=V_inc(j)-&
                (1.0/ports(excitation_port_number,2,3))*(del_t/step_vc_out)*(C_inc(j)-C_inc(j-1))
            end do
            !voltage source in TF/SF
            V_inc(port_array_size)=V_inc(port_array_size)+&
            (1.0/sqrt(ports(excitation_port_number,2,3)*ports(excitation_port_number,2,2)))*(del_t/step_vc_out)*&
            E_or_V_1D(antenna_amp,counter+0.5,0.5*step_vc_out,spread,t_spread,pulse_type,del_t,c)
            !now get incident from the voltage clear case
            incident(counter)=V_inc(port_array_size)
        end do
    end if
#endif

#ifdef use_kmax_version
    if (pbc_y+pbc_z==2) then
        do counter=1, time_steps
            E_reflected_TE(counter)=E_reflected_TE(counter)/((y_size-2.0)*(z_size-2.0))
            E_transmitted_TE(counter)=E_transmitted_TE(counter)/((y_size-2.0)*(z_size-2.0))
            E_reflected_TM(counter)=E_reflected_TM(counter)/((y_size-2.0)*(z_size-2.0))
            E_transmitted_TM(counter)=E_transmitted_TM(counter)/((y_size-2.0)*(z_size-2.0))
        end do
    end if
    if (pbc_x+pbc_z==2) then
        do counter=1, time_steps
            E_reflected_TE(counter)=E_reflected_TE(counter)/((x_size-2.0)*(z_size-2.0))
            E_transmitted_TE(counter)=E_transmitted_TE(counter)/((x_size-2.0)*(z_size-2.0))
            E_reflected_TM(counter)=E_reflected_TM(counter)/((x_size-2.0)*(z_size-2.0))
            E_transmitted_TM(counter)=E_transmitted_TM(counter)/((x_size-2.0)*(z_size-2.0))
        end do
    end if
    if (pbc_x+pbc_y==2) then
        do counter=1, time_steps
            E_reflected_TE(counter)=E_reflected_TE(counter)/((x_size-2.0)*(y_size-2.0))
            E_transmitted_TE(counter)=E_transmitted_TE(counter)/((x_size-2.0)*(y_size-2.0))
            E_reflected_TM(counter)=E_reflected_TM(counter)/((x_size-2.0)*(y_size-2.0))
            E_transmitted_TM(counter)=E_transmitted_TM(counter)/((x_size-2.0)*(y_size-2.0))
        end do
    end if

    !plane wave excitation - saved and exported but not currently used except for limits of frequency range to use
    if (plane_wave_amp>0) then
        do counter=1, time_steps
            incident(counter) = E_or_V_1D(plane_wave_amp,counter+0.0,0.0,t_spread,del_t,c,spread,pi,f_adj,imag)
        end do
    end if
    !then antenna excitation if used then this will be the source used in post processor
    if (antenna_amp>0) then
        !make a clear grid so we need current and voltage equations with sources added at the right locations
        !first determine if x,y,or z so we know which step size to use
        if (ports(excitation_port_number,1,1)==0) then
            step_vc_out=del_x
        end if
        if (ports(excitation_port_number,1,1)==1) then
            step_vc_out=del_y
        end if
        if (ports(excitation_port_number,1,1)==2) then
            step_vc_out=del_z
        end if
        !start time loop - this sim time scales with time steps, not grid size directly
        do counter=1, time_steps
            !current grid
            do j=1, port_array_size*2-1
                C_inc(j)=C_inc(j)-&
                (1.0/ports(excitation_port_number,2,2))*(del_t/step_vc_out)*(V_inc(j+1)-V_inc(j))
            end do
            !current source in TF/SF
            C_inc(port_array_size)=C_inc(port_array_size)-&
            (1.0/ports(excitation_port_number,2,2))*(del_t/step_vc_out)*&
            E_or_V_1D(antenna_amp,counter+0.0,0.0*step_vc_out,t_spread,del_t,c,spread,pi,f_adj,imag)
            !voltage grid
            do j=1, port_array_size*2-1
                V_inc(j)=V_inc(j)-&
                (1.0/ports(excitation_port_number,2,3))*(del_t/step_vc_out)*(C_inc(j)-C_inc(j-1))
            end do
            !voltage source in TF/SF
            V_inc(port_array_size)=V_inc(port_array_size)+&
            (1.0/sqrt(ports(excitation_port_number,2,3)*ports(excitation_port_number,2,2)))*(del_t/step_vc_out)*&
            E_or_V_1D(antenna_amp,counter+0.5,0.5*step_vc_out,t_spread,del_t,c,spread,pi,f_adj,imag)
            !now get incident from the voltage clear case
            incident(counter)=V_inc(port_array_size)
        end do
    end if
#endif

    !Finalize any output arrays for far field information
    do i=1, num_far_field_angles
        ang1=far_field_angles(i,1)
        ang2=far_field_angles(i,2)
        do j=1,len_far_field_arrays
            W_theta(i,j)=Wx(i,j)*cos(ang1)*cos(ang2)+Wy(i,j)*cos(ang1)*sin(ang2)-Wz(i,j)*sin(ang1)
            W_phi(i,j)=-1.0*Wx(i,j)*sin(ang2)+Wy(i,j)*cos(ang2)
            U_theta(i,j)=Ux(i,j)*cos(ang1)*cos(ang2)+Uy(i,j)*cos(ang1)*sin(ang2)-Uz(i,j)*sin(ang1)
            U_phi(i,j)=-1.0*Ux(i,j)*sin(ang2)+Uy(i,j)*cos(ang2)
            E_theta(i,j)=-1.0*sqrt(mu_0/ep_0)*W_theta(i,j)-U_phi(i,j)
            E_phi(i,j)=-1.0*sqrt(mu_0/ep_0)*W_phi(i,j)+U_theta(i,j)
        end do
    end do

    !now loop through and create the output arrays with length time_steps.
    !this is purely for convenience so that output arrays are the same length.
    !this tries to indentify where the first contributions from the far field loops is added in.
    !this is the easiest, less computationally expensive way to do this. I tried several others.
    do i=1, num_far_field_angles

        ! find first non-zero index for E_phi
        data_out_phi = len_far_field_arrays
        do j=1, len_far_field_arrays
            if (E_phi(i,j) /= 0.0+0.0*imag) then
                data_out_phi = j
                exit
            end if
        end do

        ! find first non-zero index for E_theta
        data_out_theta = len_far_field_arrays
        do j=1, len_far_field_arrays
            if (E_theta(i,j) /= 0.0+0.0*imag) then
                data_out_theta = j
                exit
            end if
        end do

        ! choose earliest signal arrival
        data_out_time(i) = min(data_out_phi, data_out_theta)

        ! ensure we have enough samples remaining
        if (data_out_time(i) > (len_far_field_arrays - time_steps + 1)) then
            data_out_time(i) = len_far_field_arrays - time_steps + 1
        end if

        ! now safely copy output window
        do j=1, time_steps
            E_theta_out(i,j) = E_theta(i, j + data_out_time(i) - 1)
            E_phi_out(i,j)   = E_phi(i, j + data_out_time(i) - 1)
        end do

    end do

#ifndef use_kmax_version
    !Save any data we want from the subroutine into a txt like file
    open (unit=10,file=filename,action="write",status="replace")
        write(10,*)  "x steps = " , x_size-1-pbc_x
        write(10,*)  "y steps = " , y_size-1-pbc_y
        write(10,*)  "z steps = " , z_size-1-pbc_z
        write(10,*)  "delta x,y,z = " , step_size_x,step_size_y,step_size_z
        write(10,*)  "time step = " , del_t
        write(10,*)  "number of time steps = " , time_steps
        write(10,*)  "minimum time steps needed (but more recommended) = ", min_steps
        write(10,*)  "freq parameter = " , f_center
        write(10,*)  "pulse type = " , pulse_type
        write(10,*)  "slice = " , slice
        write(10,*)  "slice location = " , slice_location
        write(10,*)  "video size_a = " , vid_size1
        write(10,*)  "video size_b = " , vid_size2
        write(10,*)  "number of far field angles = ", num_far_field_angles
        write(10,*) "number of internal ports = ", num_ports
    #ifndef use_spice_version
        write(10,*) "number of spice ports = ", 0
    #endif
    #ifdef use_spice_version
        write(10,*) "number of spice ports = ", num_spice_ports
    #endif

        write(10,*)  "simulation type" , pbc_x, pbc_y, pbc_z

        if (trim(excitation_type)=='plane wave') then
            write(10,*)  "incident E-field wave:"
        end if
        if (trim(excitation_type)=='antenna') then
            write(10,*)  "incident voltage wave at port ", excitation_port_number
        end if
        do i=1,time_steps
            write(10,*)  incident(i)
        end do

        if ((pbc_x+pbc_y+pbc_z==2) .and. trim(excitation_type)=='plane wave') then
            !set conditions for directions - fields market as refl and trans are for low and high directions
            !We only need light constraints here - it's assumed it's aligned on axis if made it this far.
            if ((pbc_x==0 .and. phi==0) .or. (pbc_y==0 .and. (phi>89*pi/180 .and. phi<91*pi/180)) &
            .or. (pbc_z==0 .and. theta==0)) then
                write(10,*)  "reflected wave for PBC:"
                do i=1, time_steps
                    write(10,*)  E_reflected(i)
                end do
                write(10,*)  "transmitted wave for PBC:"
                do i=1, time_steps
                    write(10,*)  E_transmitted(i)
                end do
            else
                write(10,*)  "reflected wave for PBC:"
                do i=1, time_steps
                    write(10,*)  E_transmitted(i)
                end do
                write(10,*)  "transmitted wave for PBC:"
                do i=1, time_steps
                    write(10,*)  E_reflected(i)
                end do
            end if
        end if

        do i=1,num_far_field_angles
            write(10,*)  "far field rE theta polarization at theta,phi:"
            write(10,*)  180.0/pi*far_field_angles(i,1), 180.0/pi*far_field_angles(i,2)
            do j=1,time_steps
                write(10,*)  E_theta_out(i,j)
            end do
            write(10,*)  "far field rE phi polarization at theta,phi:"
            write(10,*)  180.0/pi*far_field_angles(i,1), 180.0/pi*far_field_angles(i,2)
            do j=1,time_steps
                write(10,*)  E_phi_out(i,j)
            end do
        end do

        do i=1,num_ports
            write(10,*) "recieved voltage at port ", i, "w/ Z"
            write(10,*) sqrt(ports(i,2,2)/ports(i,2,3))
            do j=1, time_steps
                write(10,*) Voltage_out(i,j)
            end do
        end do

    #ifdef use_spice_version
        do i=1, num_spice_ports
            write(10,*) "recieved voltage at spice port ", i, "w/ Z"
            write(10,*) 1.0 !spice impedance set to 1, will need normaliziation in outputs later on
            do j=1, time_steps
                write(10,*) Spice_Voltage_out(i,j)
            end do
        end do
    #endif

        if (num_far_field_angles>0) then
            write(10,*) "far field phase correction time information:"
            do i=1, num_far_field_angles
                write(10,*) (r_for_time_relay/c) - & 
                            (data_out_time(i) - 1)*del_t
            end do
        end if
        if (num_far_field_angles>0 .and. plane_wave_amp>0) then
            write(10,*) "incident for far field phase centering:"
            do i=1, time_steps
                    write(10,*) Ex_ff_pc(i)*(cos(pol)*sin(phi)-sin(pol)*cos(theta)*cos(phi)) + &
                    Ey_ff_pc(i)*(-1*cos(pol)*cos(phi)-sin(pol)*cos(theta)*sin(phi)) + &
                    Ez_ff_pc(i)*(sin(pol)*sin(theta))
            end do
        end if

    close (10)
#endif

#ifdef use_kmax_version
    open (unit=10,file=filename,action="write",status="replace")
        write(10,*)  "x steps = " , x_size-1-pbc_x
        write(10,*)  "y steps = " , y_size-1-pbc_y
        write(10,*)  "z steps = " , z_size-1-pbc_z
        write(10,*)  "delta x,y,z = " , step_size_x,step_size_y,step_size_z
        write(10,*)  "time step = " , del_t
        write(10,*)  "number of time steps = " , time_steps
        write(10,*)  "minimum time steps needed N/A"
        write(10,*)  "freq center = " , f_center
        write(10,*)  "pulse type = N/A"
        write(10,*)  "slice = " , slice
        write(10,*)  "slice location = " , slice_location
        write(10,*)  "video size_a = " , vid_size1
        write(10,*)  "video size_b = " , vid_size2
        write(10,*)  "number of far field angles = ", num_far_field_angles
        write(10,*) "number of internal ports = ", num_ports
    #ifndef use_spice_version
        write(10,*) "number of spice ports = ", 0
    #endif
    #ifdef use_spice_version
        write(10,*) "number of spice ports = ", num_spice_ports/2
    #endif

        write(10,*)  "simulation type" , pbc_x, pbc_y, pbc_z, k_direction

        write(10,*)  "simulation is k vector method w/ k1,k2:"
        if (pbc_y+pbc_z==2) then
            write(10,*)  k_count_y, k_count_z
        else if (pbc_x+pbc_z==2) then
            write(10,*)  k_count_x, k_count_z
        else if (pbc_x+pbc_y==2) then
            write(10,*)  k_count_x, k_count_y
        end if
        write(10,*)  "mode type is 0 for TE and 1 for TM:", mode_type
                
        if (trim(excitation_type)=='plane wave') then
            write(10,*)  "incident E-field wave: (not used in post processor bc r(k))"
        end if
        if (trim(excitation_type)=='antenna') then
            write(10,*)  "incident voltage wave at port ", excitation_port_number
        end if
        do i=1, time_steps
            write(10,*)  real(incident(i)), aimag(incident(i))
        end do

        if (k_direction==1) then
            write(10,*)  "TE Refl Wave:"
            do i = 1, time_steps
                write(10,*) real(E_reflected_TE(i)), aimag(E_reflected_TE(i))
            end do
            write(10,*)  "TM Refl Wave:"
            do i = 1, time_steps
                write(10,*) real(E_reflected_TM(i)), aimag(E_reflected_TM(i))
            end do
            write(10,*)  "TE Trans Wave:"
            do i = 1, time_steps
                write(10,*) real(E_transmitted_TE(i)), aimag(E_transmitted_TE(i))
            end do
            write(10,*)  "TM Trans Wave:"
            do i = 1, time_steps
                write(10,*) real(E_transmitted_TM(i)), aimag(E_transmitted_TM(i))
            end do
        end if
        if (k_direction==0) then
            write(10,*)  "TE Refl Wave:"
            do i = 1, time_steps
                write(10,*) real(E_transmitted_TE(i)), aimag(E_transmitted_TE(i))
            end do
            write(10,*)  "TM Refl Wave:"
            do i = 1, time_steps
                write(10,*) real(E_transmitted_TM(i)), aimag(E_transmitted_TM(i))
            end do
            write(10,*)  "TE Trans Wave:"
            do i = 1, time_steps
                write(10,*) real(E_reflected_TE(i)), aimag(E_reflected_TE(i))
            end do
            write(10,*)  "TM Trans Wave:"
            do i = 1, time_steps
                write(10,*) real(E_reflected_TM(i)), aimag(E_reflected_TM(i))
            end do
        end if

        do i = 1, num_far_field_angles
            write(10,*)  "far field rE theta polarization at theta,phi:"
            write(10,*)  180.0/pi*far_field_angles(i,1), 180.0/pi*far_field_angles(i,2)
            do j = 1, time_steps
                write(10,*) real(E_theta_out(i,j)), aimag(E_theta_out(i,j))
            end do
            write(10,*)  "far field rE phi polarization at theta,phi:"
            write(10,*)  180.0/pi*far_field_angles(i,1), 180.0/pi*far_field_angles(i,2)
            do j = 1, time_steps
                write(10,*) real(E_phi_out(i,j)), aimag(E_phi_out(i,j))
            end do
        end do

        do i = 1, num_ports
            write(10,*) "recieved voltage at port ", i, "w/ Z"
            write(10,*) sqrt(ports(i,2,2)/ports(i,2,3))
            do j = 1, time_steps
                write(10,*) real(Voltage_out(i,j)), aimag(Voltage_out(i,j))
            end do
        end do

    #ifdef use_spice_version
        do i=1, num_spice_ports/2
            !1,2 real,imag currents for same fdtd E-port that was duplicated by the user on purpose
            !These are all real quantities here though but correspond to real and imaginary parts of the current and thus voltage actual
            write(10,*) "recieved voltage at spice port ", i, "w/ Z"
            write(10,*) 1.0 !spice impedance set to 1, will need normaliziation in outputs later on
            do j=1, time_steps
                write(10,*) Spice_Voltage_out(2*i-1,j), Spice_Voltage_out(2*i,j)
            end do
        end do
    #endif

        if (num_far_field_angles>0) then
            write(10,*) "far field phase correction time information:"
            do i=1, num_far_field_angles
                write(10,*) (r_for_time_relay/c) - & 
                            (data_out_time(i) - 1)*del_t
            end do
        end if
        if (num_far_field_angles>0 .and. plane_wave_amp>0) then
            write(10,*) "incident for far field phase centering (Ex,Ey,Ez):"
            do i=1, time_steps
                write(10,*) real(Ex_ff_pc(i)),aimag(Ex_ff_pc(i))
            end do
            do i=1, time_steps
                write(10,*) real(Ey_ff_pc(i)),aimag(Ey_ff_pc(i))
            end do
            do i=1, time_steps
                write(10,*) real(Ez_ff_pc(i)),aimag(Ez_ff_pc(i))
            end do
        end if

    close(10)
#endif

    !Save the video arrays into binary fiels if we asked for videos
    if (video_on==1) then
        open(11, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Ex.bin", form="unformatted",action="write",status="replace")
            write(11) Ex_video
        close(11)
        open(12, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Hx.bin", form="unformatted",action="write",status="replace")
            write(12) Hx_video
        close(12)
        open(13, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Ey.bin", form="unformatted",action="write",status="replace")
            write(13) Ey_video
        close(13)
        open(14, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Hy.bin", form="unformatted",action="write",status="replace")
            write(14) Hy_video
        close(14)
        open(15, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Ez.bin", form="unformatted",action="write",status="replace")
            write(15) Ez_video
        close(15)
        open(16, file=filename(1:LEN_TRIM(filename)-4)//"_"//"video_Hz.bin", form="unformatted",action="write",status="replace")
            write(16) Hz_video
        close(16)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! END POST PROCESSING AND EXPORT OF DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call system_clock(clock_time_end)
    write(*, '("Finished post processing. Time to complete was ", F0.3, " seconds.")') &
          real(clock_time_end - clock_time_start) / real(clock_rate)

    !Define any functions below for convenience
    contains

#ifndef use_kmax_version
        function Inc(weight_f,x_f,y_f,z_f,t_spread_f,spread_f,counter_f,theta_f,phi_f,del_x_f,del_y_f,del_z_f,x_delay_f,y_delay_f,z_delay_f,pulse_type_f,del_t_f,c_f)
            !$acc routine seq
            implicit none
            real, intent(in) :: weight_f, x_f, y_f, z_f, t_spread_f, spread_f, theta_f, phi_f, del_x_f, del_y_f, del_z_f, del_t_f, c_f, counter_f
            integer, intent(in) :: x_delay_f, y_delay_f, z_delay_f, pulse_type_f
            real :: Inc
            real :: arg_f, term_x_f, term_y_f, term_z_f
            term_x_f = sin(theta_f) * cos(phi_f) * (x_f - x_delay_f) * del_x_f
            term_y_f = sin(theta_f) * sin(phi_f) * (y_f - y_delay_f) * del_y_f
            term_z_f = cos(theta_f) * (z_f - z_delay_f) * del_z_f
            arg_f = (t_spread_f - (counter_f * del_t_f) + (1.0/c_f) * (term_x_f + term_y_f + term_z_f)) / spread_f
            if (pulse_type_f == 1) then
                Inc = weight_f * exp(-0.5 * arg_f**2)
            end if
            if (pulse_type_f == 2) then
                Inc = -1.0 * weight_f * arg_f * exp(0.5) * exp(-0.5 * arg_f**2)
            end if
        end function Inc

        function E_or_V_1D(amp_f2,time_counter_f2,position_delay_f2,spread_f2,t_spread_f2,pulse_type_f2,del_t_f2,c_f2)
            !$acc routine seq
            implicit none
            real, intent(in) :: amp_f2, position_delay_f2,del_t_f2,c_f2,spread_f2,t_spread_f2,time_counter_f2
            integer, intent(in) :: pulse_type_f2
            real E_or_V_1D
            if (pulse_type_f2==1) then
                E_or_V_1D = amp_f2*EXP(-0.5*((t_spread_f2-((time_counter_f2-1.0)*del_t_f2)-1.0/c_f2*position_delay_f2)/spread_f2)**2)
            end if
            if (pulse_type_f2==2) then
                E_or_V_1D = amp_f2*-1.0*((t_spread_f2-((time_counter_f2-1.0)*del_t_f2)-1.0/c_f2*position_delay_f2)/spread_f2)*EXP(0.5)*&
                EXP(-0.5*((t_spread_f2-((time_counter_f2-1.0)*del_t_f2)-1.0/c_f2*position_delay_f2)/spread_f2)**2)
            end if
        end function E_or_V_1D
#endif

#ifdef use_kmax_version
        function Inc(weight_f,x_f,k_y_f,d_var_y_f,k_z_f,d_var_z_f,f_var_f,t_spread_f,spread_f,counter_f,del_x_f,del_t_f,c_f,pi_f,imag_f)

            !$acc routine seq
            implicit none

            real, intent(in) :: x_f
            real, intent(in) :: k_y_f
            real, intent(in) :: d_var_y_f
            real, intent(in) :: k_z_f
            real, intent(in) :: d_var_z_f
            real, intent(in) :: weight_f 
            real, intent(in) :: f_var_f

            real, intent(in) :: t_spread_f,spread_f,counter_f,del_x_f,del_t_f,c_f,pi_f
            complex, intent(in) :: imag_f

            complex Inc

            real :: arg1_f
            complex :: arg2_f

            arg1_f=weight_f*EXP(-0.5*((t_spread_f-((counter_f)*del_t_f)+1.0/c_f*(x_f*del_x_f))/spread_f)**2)
            arg2_f=CEXP(imag_f*2.0*pi_f*f_var_f*((counter_f)*del_t_f+1.0/c_f*(x_f*del_x_f))-1.0*imag_f*(k_y_f*d_var_y_f+k_z_f*d_var_z_f))

            Inc = arg1_f*arg2_f

        end function Inc

        function E_or_V_1D(antenna_amp_f,time_counter_f,position_delay_f,t_spread_f,del_t_f,c_f,spread_f,pi_f,f_adj_f,imag_f)
            !$acc routine seq
            implicit none
            complex E_or_V_1D
            real, intent(in) :: time_counter_f
            real, intent(in) :: position_delay_f
            
            real, intent(in) :: antenna_amp_f, t_spread_f,del_t_f,c_f,spread_f,pi_f,f_adj_f
            complex, intent(in) :: imag_f

            E_or_V_1D = antenna_amp_f*EXP(-0.5*((t_spread_f-((time_counter_f-1.0)*del_t_f)-1.0/c_f*position_delay_f)/spread_f)**2)*&
            CEXP(imag_f*2*pi_f*f_adj_f*((time_counter_f-1.0)*del_t_f-1.0/c_f*position_delay_f))

        end function E_or_V_1D
#endif

End Program fdtd
