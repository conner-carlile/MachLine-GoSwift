module panel_solver_mod

    use helpers_mod
    use json_mod
    use json_xtnsn_mod
    use panel_mod
    use vertex_mod
    use surface_mesh_mod
    use flow_mod
    use math_mod
    use linalg_mod
    use preconditioners_mod
    use sparsity_mod
    use sort_mod

    implicit none


    type panel_solver

        real :: corrected_M_inf 
        character(len=:),allocatable :: formulation, pressure_for_forces, matrix_solver, preconditioner
        logical :: incompressible_rule, isentropic_rule, second_order_rule, slender_rule, linear_rule
        logical :: morino, write_A_and_b
        logical :: compressible_correction, prandtl_glauert, karman_tsien, laitone
        type(dod),dimension(:,:),allocatable :: dod_info, wake_dod_info
        type(flow) :: freestream
        real :: norm_res, max_res, tol, rel
        real,dimension(3) :: C_F, inner_flow
        real,dimension(:,:),allocatable :: A
        real,dimension(:), allocatable :: b
        integer :: N, wake_start, N_cells, block_size, max_iterations, N_unknown, N_d_unknown, N_s_unknown

        contains

            ! Initialization
            procedure :: init => panel_solver_init
            procedure :: parse_solver_settings => panel_solver_parse_solver_settings
            procedure :: parse_processing_settings => panel_solver_parse_processing_settings
            procedure :: init_dirichlet => panel_solver_init_dirichlet
            procedure :: calc_domains_of_dependence => panel_solver_calc_domains_of_dependence

            ! Solve
            procedure :: solve => panel_solver_solve
            procedure :: calc_source_strengths => panel_solver_calc_source_strengths
            procedure :: update_system_row => panel_solver_update_system_row
            procedure :: calc_body_influences => panel_solver_calc_body_influences
            procedure :: calc_wake_influences => panel_solver_calc_wake_influences
            procedure :: solve_system => panel_solver_solve_system

            ! Post-processing
            procedure :: calc_cell_velocities => panel_solver_calc_cell_velocities
            procedure :: calc_point_velocities => panel_solver_calc_point_velocities
            procedure :: calc_surface_potentials => panel_solver_calc_surface_potentials
            procedure :: calc_pressures => panel_solver_calc_pressures
            procedure :: subsonic_pressure_correction => panel_solver_subsonic_pressure_correction
            procedure :: calc_crit_mach => panel_solver_calc_crit_mach
            procedure :: calc_forces => panel_solver_calc_forces
            procedure :: update_report => panel_solver_update_report
            procedure :: export_potential_slice => panel_solver_export_potential_slice
            procedure :: get_panel_induced_potentials => panel_solver_get_panel_induced_potentials

    end type panel_solver


contains


    subroutine panel_solver_init(this, solver_settings, processing_settings, body, freestream, control_point_file)

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings, processing_settings
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(inout) :: freestream
        character(len=:),allocatable,intent(in) :: control_point_file

        ! Store
        this%freestream = freestream

        ! Get solver settings
        call this%parse_solver_settings(solver_settings)

        ! Get post-processing settings
        call this%parse_processing_settings(processing_settings)

        ! Initialize based on formulation
        if (this%morino .or. this%formulation == 'source-free') then
            call this%init_dirichlet(solver_settings, body)
        end if
        
        ! Write out control point geometry
        if (control_point_file /= 'none') then
            call body%write_control_points(control_point_file, solved=.false.)
        end if

        ! Calculate domains of dependence
        call this%calc_domains_of_dependence(body)

    end subroutine panel_solver_init


    subroutine panel_solver_parse_solver_settings(this, solver_settings)
        ! Parses the solver settings from the input

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value),pointer, intent(in) :: solver_settings
        
        ! Parse settings
        call json_xtnsn_get(solver_settings, 'formulation', this%formulation, 'morino')        
        call json_xtnsn_get(solver_settings, 'matrix_solver', this%matrix_solver, 'LU')
        call json_xtnsn_get(solver_settings, 'block_size', this%block_size, 400)
        call json_xtnsn_get(solver_settings, 'tolerance', this%tol, 1e-10)
        call json_xtnsn_get(solver_settings, 'relaxation', this%rel, 0.9)
        call json_xtnsn_get(solver_settings, 'max_iterations', this%max_iterations, 1000)
        call json_xtnsn_get(solver_settings, 'preconditioner', this%preconditioner, 'NONE')
        call json_xtnsn_get(solver_settings, 'write_A_and_b', this%write_A_and_b, .false.)
        this%morino = this%formulation == 'morino'
        
    end subroutine panel_solver_parse_solver_settings


    subroutine panel_solver_parse_processing_settings(this, processing_settings)
        ! Parses the post-processing settings from the input

        implicit none
        
        class(panel_solver), intent(inout) :: this
        type(json_value),pointer, intent(in) :: processing_settings
        
        ! Get incompressible/isentropic pressure rules
        if (this%freestream%M_inf > 0.) then

            ! Isentropic is default for M > 0
            call json_xtnsn_get(processing_settings, 'pressure_rules.isentropic', this%isentropic_rule, .true.)

            ! Check for incompressible rule
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .false.)

            ! Notify user if we're throwing out the incompressible rule
            if (this%incompressible_rule) then
                write(*,*) "!!! The incompressible pressure rule cannot be used for M > 0. Switching to the isentropic rule."
                this%incompressible_rule = .false.
                this%isentropic_rule = .true.
            end if

        else if (this%freestream%M_inf == 0.) then

            ! Incompressible is deafult
            call json_xtnsn_get(processing_settings, 'pressure_rules.incompressible', this%incompressible_rule, .true.)

            ! Check for isentropic
            call json_xtnsn_get(processing_settings, 'pressure_rules.isentropic', this%isentropic_rule, .false.)

            ! Notify user if pressure rule applied is changed based on selected freestream Mach number
            if (this%isentropic_rule) then
                write(*,*) "!!! The isentropic pressure rule cannot be used for M = 0. Switching to the incompressible rule."
                this%isentropic_rule = .false.
                this%incompressible_rule = .true.
            end if

        else

            ! MachLine does not support a negative freestream Mach number
            write(*,*) "!!! A negative freestream Mach number is not allowed. Quitting..."
            stop

        end if

        ! Get other pressure rules (these are applicable for all Mach numbers)
        call json_xtnsn_get(processing_settings, 'pressure_rules.second-order', this%second_order_rule, .false.)
        call json_xtnsn_get(processing_settings, 'pressure_rules.slender-body', this%slender_rule, .false.)
        call json_xtnsn_get(processing_settings, 'pressure_rules.linear', this%linear_rule, .false.)
        
        ! Get information for pressure corrections
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.correction_mach_number', this%corrected_M_inf, 0.0)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.prandtl-glauert', this%prandtl_glauert, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.karman-tsien', this%karman_tsien, .false.)
        call json_xtnsn_get(processing_settings, 'subsonic_pressure_correction.laitone', this%laitone, .false.)
        
        ! Check the correction Mach number is subsonic
        if (this%corrected_M_inf < 0.0 .or. this%corrected_M_inf >= 1.0) then
            write(*,*) "!!! The pressure correction Mach number must be greater than zero and subsonic. Quitting..."
            stop
        end if

        ! Check freestream Mach number is set to 0 if pressure correction is selected
        if (((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) .and. (this%freestream%M_inf /= 0.0)) then
            write(*,*) "!!! In order to apply a subsonic pressure correction, the freestream Mach number must be set to '0'."
            write(*,*) "!!! Include the desired Mach number for the correction calculations as 'correction_mach_number' under"
            write(*,*) "!!! 'post-processing'-'subsonic_pressure_correction'. Quitting..."
            stop
        end if

        ! Check to see if user selected incompressible rule when applying the compressible pressure corrections
        if (((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) .and. .not. this%incompressible_rule) then
            this%incompressible_rule = .true.
        end if 

        ! Get which pressure rule will be used for force calculation
        if (this%incompressible_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'incompressible')
        else if (this%isentropic_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'isentropic')
        else if (this%second_order_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'second-order')
        else if (this%linear_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'linear')
        else if (this%slender_rule) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'slender-body')
        else if (this%prandtl_glauert) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'prandtl-glauert')
        else if (this%karman_tsien) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'karman-tsien')
        else if (this%laitone) then
            call json_xtnsn_get(processing_settings, 'pressure_for_forces', this%pressure_for_forces, 'laitone')
        end if
        
    end subroutine panel_solver_parse_processing_settings


    subroutine panel_solver_init_dirichlet(this, solver_settings, body)
        ! Initializes the solver to use one of the Dirichlet formulations

        implicit none

        class(panel_solver),intent(inout) :: this
        type(json_value),pointer,intent(in) :: solver_settings
        type(surface_mesh),intent(inout) :: body

        real :: offset

        ! Get offset
        call json_xtnsn_get(solver_settings, 'control_point_offset', offset, 1e-5)
        
        ! Place control points
        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Placing control points using offset of ", offset, "..."

        ! Place control points inside the body
        if (this%morino .or. this%formulation == 'source-free') then
            call body%place_interior_control_points(offset)
        end if

        ! Determine size of linear system
        if (body%asym_flow) then
            this%N = body%N_cp*2
        else
            this%N = body%N_cp
        end if

        ! Calculate target inner flow
        this%inner_flow = this%freestream%c_hat_g
        if (this%formulation == 'source-free') then
            this%inner_flow = this%inner_flow - matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        if (verbose) write(*,'(a, i6, a)') "Done. Placed", body%N_cp, " control points."
    
    end subroutine panel_solver_init_dirichlet


    subroutine panel_solver_calc_domains_of_dependence(this, body)
        ! Determines the domains of dependence for each control point based on the freestream condition

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, stat
        real,dimension(3) :: vert_loc, mirrored_vert_loc
        logical,dimension(:),allocatable :: wake_verts_in_dod, mirrored_wake_verts_in_dod, verts_in_dod, mirrored_verts_in_dod

        ! For asymmetric flow on a mirrored mesh, all domains of dependence must be calculated. There are no shortcuts.
        ! For symmetric flow on a mirrored mesh, domains of dependence will be the same between mirrored panels and mirrored
        ! control points. So, we just need to calculate the DoD for mirrored control points, and then we're good.

        if (this%freestream%supersonic .and. verbose) write(*,'(a)',advance='no') "     Calculating domains of dependence..."

        ! Allocate arrays for domain of dependence information for the body
        if (body%mirrored) then
            allocate(this%dod_info(2*body%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(2*body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")

            allocate(mirrored_verts_in_dod(2*body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        else
            allocate(this%dod_info(body%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        end if

        ! Allocate arrays for domain of dependence information for the wake
        if (body%mirrored .and. .not. body%asym_flow) then ! This is the only case where the wake is mirrored
            allocate(this%wake_dod_info(2*body%wake%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(wake_verts_in_dod(2*body%wake%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")

            allocate(mirrored_wake_verts_in_dod(2*body%wake%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        else
            allocate(this%wake_dod_info(body%wake%N_panels, this%N), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(wake_verts_in_dod(body%wake%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        end if

        ! If the freestream is subsonic, these don't need to be checked
        if (this%freestream%supersonic) then

            ! Loop through control points
            !$OMP parallel do private(i, vert_loc, mirrored_vert_loc, verts_in_dod, mirrored_verts_in_dod) &
            !$OMP & private(wake_verts_in_dod, mirrored_wake_verts_in_dod)
            do j=1,body%N_cp

                ! Initialize for this control point
                verts_in_dod = .false.
                wake_verts_in_dod = .false.
                if (body%asym_flow) then
                    mirrored_verts_in_dod = .false.
                end if
                if (body%mirrored .and. .not. body%asym_flow) then
                    mirrored_wake_verts_in_dod = .false.
                end if

                ! Loop through body vertices
                do i=1,body%N_verts

                    vert_loc = body%vertices(i)%loc

                    ! Original vertex and original control point
                    verts_in_dod(i) = this%freestream%point_in_dod(vert_loc, body%cp(:,j))

                    if (body%mirrored) then

                        mirrored_vert_loc = mirror_across_plane(vert_loc, body%mirror_plane)

                        ! Mirrored vertex and original control point
                        verts_in_dod(i+body%N_verts) = this%freestream%point_in_dod(mirrored_vert_loc, &
                                                                                           body%cp(:,j))

                        if (body%asym_flow) then

                            ! Original vertex and mirrored control point
                            mirrored_verts_in_dod(i) = this%freestream%point_in_dod(vert_loc, body%cp_mirrored(:,j))

                            ! Mirrored vertex and mirrored control point
                            mirrored_verts_in_dod(i+body%N_verts) = this%freestream%point_in_dod(mirrored_vert_loc, &
                                                                                                 body%cp_mirrored(:,j))

                        end if
                    end if
                end do

                ! Loop through wake vertices
                do i=1,body%wake%N_verts

                    vert_loc = body%wake%vertices(i)%loc

                    ! Original vertex and original control point
                    wake_verts_in_dod(i) = this%freestream%point_in_dod(vert_loc, body%cp(:,j))

                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Original vertex and mirrored control point
                            wake_verts_in_dod(i) = this%freestream%point_in_dod(vert_loc, body%cp_mirrored(:,j))

                        else

                            ! Mirrored vertex and original control point
                            mirrored_vert_loc = mirror_across_plane(vert_loc, body%mirror_plane)
                            mirrored_wake_verts_in_dod(i+body%wake%N_verts) = this%freestream%point_in_dod(mirrored_vert_loc, &
                                                                                                           body%cp(:,j))
                        end if
                    end if
                end do

                ! Loop through body panels
                do i=1,body%N_panels

                    ! Original panel and original control point
                    this%dod_info(i,j) = body%panels(i)%check_dod(body%cp(:,j), this%freestream, verts_in_dod)

                    if (body%mirrored) then

                        ! Check DoD for mirrored panel and original control point
                        this%dod_info(i+body%N_panels,j) = body%panels(i)%check_dod(body%cp(:,j), this%freestream, &
                                                                                    verts_in_dod, &
                                                                                    .true., body%mirror_plane)

                        if (body%asym_flow) then

                            ! Check DoD for original panel and mirrored control point
                            this%dod_info(i,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mirrored(:,j), this%freestream, &
                                                                                    mirrored_verts_in_dod)

                            ! Check DoD for mirrored panel and mirrored control point
                            this%dod_info(i+body%N_panels,j+body%N_cp) = body%panels(i)%check_dod(body%cp_mirrored(:,j), &
                                                                                                  this%freestream, &
                                                                                                  mirrored_verts_in_dod, &
                                                                                                  .true., body%mirror_plane)

                        end if
                    end if
                end do

                ! Loop through wake panels
                do i=1,body%wake%N_panels

                    ! Check DoD for panel and original control point
                    this%wake_dod_info(i,j) = body%wake%panels(i)%check_dod(body%cp(:,j), this%freestream, wake_verts_in_dod)

                    if (body%mirrored) then

                        if (body%asym_flow) then

                            ! Check DoD for panel and mirrored control point
                            this%wake_dod_info(i,j+body%N_cp) = body%wake%panels(i)%check_dod(body%cp_mirrored(:,j), &
                                                                                              this%freestream, &
                                                                                              wake_verts_in_dod)

                        else

                            ! Check DoD for mirrored panel and original control point
                            this%wake_dod_info(i+body%wake%N_panels,j) = body%wake%panels(i)%check_dod(body%cp(:,j), &
                                                                                                       this%freestream, &
                                                                                                       mirrored_wake_verts_in_dod, &
                                                                                                       .true., body%mirror_plane)

                        end if
                    end if
                end do

            end do

            if (verbose) write(*,*) "Done."
        end if
    
    end subroutine panel_solver_calc_domains_of_dependence


    subroutine panel_solver_solve(this, body)
        ! Calls the relevant subroutine to solve the case based on the selected formulation

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        ! Calculate source strengths
        call this%calc_source_strengths(body)

        ! Calculate body influences
        call this%calc_body_influences(body)

        ! Calculate wake influences
        if (body%wake%N_panels > 0) then
            call this%calc_wake_influences(body)
        end if

        ! Solve the linear system
        call this%solve_system(body)

        ! Calculate velocities
        call this%calc_cell_velocities(body)
        if (doublet_order == 2) call this%calc_point_velocities(body)
        
        ! Calculate potentials
        call this%calc_surface_potentials(body)

        ! Calculate pressures
        call this%calc_pressures(body)

        ! Calculate forces
        call this%calc_forces(body)

    end subroutine panel_solver_solve


    subroutine panel_solver_calc_source_strengths(this, body)
        ! Calculates the necessary source strengths

        implicit none

        class(panel_solver),intent(in) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: N_sigma, i, stat

        ! Determine number of source strengths
        if (source_order == 0) then

            if (body%asym_flow) then
                N_sigma = body%N_panels*2
            else
                N_sigma = body%N_panels
            end if

        else if (source_order == 1) then

            if (body%asym_flow) then
                N_sigma = body%N_verts*2
            else
                N_sigma = body%N_verts
            end if

        end if

        ! Allocate source strength array (yes this does need to be allocated for the source-free formulation)
        allocate(body%sigma(N_sigma), source=0., stat=stat)
        call check_allocation(stat, "source strength array")

        ! Set source strengths
        if (this%morino) then

            if (verbose) write(*,'(a)',advance='no') "     Calculating source strengths..."

            ! Use panel normals
            if (source_order == 0) then

                ! Loop through panels
                do i=1,body%N_panels

                    ! Existing panels
                    body%sigma(i) = -inner(body%panels(i)%n_g, this%freestream%c_hat_g)

                    ! Mirrored panels for asymmetric flow
                    if (body%asym_flow) then
                        body%sigma(i+body%N_panels) = -inner(body%panels(i)%n_g_mir, this%freestream%c_hat_g)
                    end if
                end do

            ! Use vertex normals
            else if (source_order == 1) then

                ! Loop through vertices
                do i=1,body%N_verts

                    ! Existing vertices
                    body%sigma(i) = -inner(body%vertices(i)%n_g, this%freestream%c_hat_g)

                    ! Mirrored panels for asymmetric flow
                    if (body%asym_flow) then
                        body%sigma(i+body%N_verts) = -inner(body%vertices(i)%n_g_mir, this%freestream%c_hat_g)
                    end if
                end do

            end if

            if (verbose) write(*,*) "Done."
        end if
    
    end subroutine panel_solver_calc_source_strengths


    subroutine panel_solver_update_system_row(this, body, A_row, phi_cp_s, i_panel, source_inf, doublet_inf, mirrored_panel)
        ! Updates the linear system with the source and doublet influences

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        real,dimension(:),allocatable,intent(inout) :: A_row
        real,intent(inout) :: phi_cp_s
        integer,intent(in) :: i_panel
        real,dimension(:),allocatable,intent(in) :: source_inf, doublet_inf
        logical,intent(in) :: mirrored_panel

        integer :: k

        ! Add source influence (if sources are present)
        if (this%morino) then

            ! Constant
            if (source_order == 0) then
                if (mirrored_panel) then
                    phi_cp_s = phi_cp_s + source_inf(1)*body%sigma(i_panel+body%N_panels)
                else
                    phi_cp_s = phi_cp_s + source_inf(1)*body%sigma(i_panel)
                end if

            ! Linear
            else
                do k=1,size(body%panels(i_panel)%i_vert_s)
                    if (mirrored_panel) then
                        phi_cp_s = phi_cp_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k)+body%N_cp)
                    else
                        phi_cp_s = phi_cp_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k))
                    end if
                end do
            end if
        end if

        ! Add doublet influence
        ! This method is the same for linear and quadratic doublets
        ! Loop through panel vertices
        do k=1,size(body%panels(i_panel)%i_vert_d)
            if (mirrored_panel) then
                A_row(body%panels(i_panel)%i_vert_d(k)+body%N_cp) = A_row(body%panels(i_panel)%i_vert_d(k)+body%N_cp) &
                                                                    + doublet_inf(k)
            else
                A_row(body%panels(i_panel)%i_vert_d(k)) = A_row(body%panels(i_panel)%i_vert_d(k)) + doublet_inf(k)
            end if
        end do
    
    end subroutine panel_solver_update_system_row


    subroutine panel_solver_calc_body_influences(this, body)
        ! Calculates the influence of the body on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k, stat
        real,dimension(:),allocatable :: source_inf, doublet_inf, A_i, A_i_mir
        real :: phi_cp_s, phi_cp_s_mir
        real,dimension(3) :: x

        ! Allocate space for inner potential calculations
        allocate(body%phi_cp_sigma(this%N), source=0., stat=stat)
        call check_allocation(stat, "induced potential vector")

        ! Allocate AIC matrix
        allocate(this%A(this%N, this%N), source=0., stat=stat)
        call check_allocation(stat, "AIC matrix")

        ! Allocate b vector
        allocate(this%b(this%N), source=0., stat=stat)
        call check_allocation(stat, "b vector")

        ! Allocate row of A
        allocate(A_i(this%N))
        if (body%asym_flow) then
            allocate(A_i_mir(this%N))
        end if

        if (verbose) write(*,'(a)',advance='no') "     Calculating body influences..."
        
        ! Parameter used for calculating inner potential for the source-free formulation
        if (this%formulation == 'source-free') then
            x = matmul(this%freestream%B_mat_g_inv, this%freestream%c_hat_g)
        end if

        ! Calculate source and doublet influences from body on each control point
        !$OMP parallel do private(j, source_inf, doublet_inf, k, A_i, A_i_mir, phi_cp_s, phi_cp_s_mir) &
        !$OMP schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            phi_cp_s = 0.
            if (body%asym_flow) then
                A_i_mir = 0.
                phi_cp_s_mir = 0.
            end if

            ! Loop through panels
            do j=1,body%N_panels

                ! Existing panel on existing control point
                if (this%dod_info(j,i)%in_dod) then
                    
                    ! Calculate influence
                    call body%panels(j)%calc_potentials(body%cp(:,i), this%freestream, this%dod_info(j,i), .false., &
                                                        source_inf, doublet_inf)

                    ! Add influence
                    call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .false.)

                end if

                ! Calculate mirrored influences
                if (body%mirrored) then

                    if (body%asym_flow) then

                        if (body%vertices(i)%mirrored_is_unique) then

                            ! Recalculate influence of mirrored panel on mirrored control point for compressible flow
                            ! If the flow is symmetric or incompressible, this will be the same as already calculated
                            if (this%dod_info(j+body%N_panels,i+body%N_cp)%in_dod) then
                                if (.not. this%freestream%incompressible) then
                                    call body%panels(j)%calc_potentials(body%cp_mirrored(:,i), this%freestream, &
                                                                        this%dod_info(j+body%N_panels,i+body%N_cp), .true., &
                                                                        source_inf, doublet_inf)
                                end if

                                ! Add influence of mirrored panel on mirrored control point
                                call this%update_system_row(body, A_i_mir, phi_cp_s_mir, j, source_inf, doublet_inf, .true.)
                            end if

                        end if

                        ! Calculate influence of existing panel on mirrored control point
                        if (this%dod_info(j,i+body%N_cp)%in_dod) then
                            if (body%vertices(i)%mirrored_is_unique .or. this%freestream%incompressible) then
                                call body%panels(j)%calc_potentials(body%cp_mirrored(:,i), this%freestream, &
                                                                    this%dod_info(j,i+body%N_cp), .false., &
                                                                    source_inf, doublet_inf)
                            end if

                            ! Add influence of existing panel on mirrored control point
                            if (body%vertices(i)%mirrored_is_unique) then
                                call this%update_system_row(body, A_i_mir, phi_cp_s_mir, j, source_inf, doublet_inf, .false.)
                            end if
                        end if

                        ! Recalculate mirrored->existing influences for compressible flow
                        if (this%dod_info(j+body%N_panels,i)%in_dod) then
                            if (.not. this%freestream%incompressible) then
                                call body%panels(j)%calc_potentials(body%cp(:,i), this%freestream, &
                                                                    this%dod_info(j+body%N_panels,i), .true., source_inf, &
                                                                    doublet_inf)
                            end if

                            ! Add influence of mirrored panel on existing control point
                            call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .true.)
                        end if

                    else

                        ! Calculate influence of existing panel on mirrored control point
                        ! This is the same as the influence of the mirrored panel on the existing control point,
                        ! even for compressible flow, since we know the flow is symmetric here
                        if (this%dod_info(j+body%N_panels,i)%in_dod) then
                            call body%panels(j)%calc_potentials(body%cp_mirrored(:,i), this%freestream, &
                                                                this%dod_info(j+body%N_panels,i), &
                                                                .false., source_inf, doublet_inf)

                            ! Add influence of mirrored panel on existing control point
                            call this%update_system_row(body, A_i, phi_cp_s, j, source_inf, doublet_inf, .false.)
                        end if

                    end if
                end if

            end do

            ! Update A matrix with rows
            !$OMP critical
            this%A(i,:) = A_i

            ! Update for mirrored points
            if (body%asym_flow) then

                ! Enforce doublet strength matching (i.e. for non-unique, mirrored control points, the
                ! doublet strengths must be the same). The RHS for these rows should still be zero.
                if (.not. body%vertices(i)%mirrored_is_unique) then
                    this%A(i+body%N_cp,i) = 1.
                    this%A(i+body%N_cp,i+body%N_cp) = -1.

                ! Otherwise, add influences
                else
                    this%A(i+body%N_cp,:) = A_i_mir
                end if
            end if

            ! Set potential for morino formulation
            if (this%morino) then
                body%phi_cp_sigma(i) = phi_cp_s

                ! Set for unique mirrored control  points
                if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then
                    body%phi_cp_sigma(i+body%N_cp) = phi_cp_s_mir
                end if

            ! Set target potential for source-free formulation
            else
                this%b(i) = -inner(x, body%cp(:,i))

                ! Set for unique mirrored control  points
                if (body%asym_flow .and. body%vertices(i)%mirrored_is_unique) then
                    this%b(i+body%N_cp) = -inner(x, body%cp_mirrored(:,i))
                end if
            end if
            !$OMP end critical

        end do

        ! Clean up memory
        deallocate(this%dod_info)

        if (verbose) write(*,*) "Done."
    
    end subroutine panel_solver_calc_body_influences


    subroutine panel_solver_calc_wake_influences(this, body)
        ! Calculates the influence of the wake on the control points

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, j, k
        real,dimension(:),allocatable ::  doublet_inf, source_inf, A_i, A_i_mir

        ! Calculate influence of wake
        if (verbose) write(*,'(a)',advance='no') "     Calculating wake influences..."

        ! Allocate A rows
        allocate(A_i(this%N))
        if (body%asym_flow) then
            allocate(A_i_mir(this%N))
        end if

        ! Loop through control points
        !$OMP parallel do private(j, source_inf, doublet_inf, k, A_i, A_i_mir) schedule(dynamic)
        do i=1,body%N_cp

            ! Initialize
            A_i = 0.
            if (body%asym_flow) then
                A_i_mir = 0.
            end if

            ! Get doublet influence from wake
            ! Note that for the wake, in the case of mirrored meshes with asymmetric flow, the mirrored wake panels have actually been created.
            ! In this case, there are technically no mirrored panels, and this loop will cycle through both existing and mirrored panels.
            ! For symmetric flow, mirrored panels still need to be added as before.
            do j=1,body%wake%N_panels

                ! Caclulate influence of existing panel on existing control point
                call body%wake%panels(j)%calc_potentials(body%cp(:,i), this%freestream, this%wake_dod_info(j,i), .false., &
                                                         source_inf, doublet_inf)

                ! Add influence
                do k=1,size(body%wake%panels(j)%i_vert_d)
                    A_i(body%wake%panels(j)%i_vert_d(k)) = A_i(body%wake%panels(j)%i_vert_d(k)) + doublet_inf(k)
                end do

                ! Get influence on mirrored control point
                if (body%mirrored) then

                    if (body%asym_flow) then

                        ! Calculate influence of existing panel on mirrored point
                        call body%wake%panels(j)%calc_potentials(body%cp_mirrored(:,i), this%freestream, &
                                                                 this%wake_dod_info(j,i+body%N_cp), .false., &
                                                                 source_inf, doublet_inf)

                        ! Add influence
                        if (body%vertices(i)%mirrored_is_unique) then
                            do k=1,size(body%wake%panels(j)%i_vert_d)
                                A_i_mir(body%wake%panels(j)%i_vert_d(k)) = A_i_mir(body%wake%panels(j)%i_vert_d(k)) + doublet_inf(k)
                            end do
                        end if

                    else

                        ! Calculate influence of existing panel on mirrored control point
                        ! This is the same as the influence of a mirrored panel on an existing control point
                        ! even for compressible flow, since we know the flow is symmetric here
                        call body%wake%panels(j)%calc_potentials(body%cp_mirrored(:,i), this%freestream, &
                                                                 this%wake_dod_info(j+body%wake%N_panels,i), .false., & ! No, this is not the DoD for this computation; yes, it is equivalent
                                                                 source_inf, doublet_inf)

                        ! Add influence
                        do k=1,size(body%wake%panels(j)%i_vert_d)
                            A_i(body%wake%panels(j)%i_vert_d(k)) = A_i(body%wake%panels(j)%i_vert_d(k)) + doublet_inf(k)
                        end do

                    end if

                end if
            end do

            ! Update rows of A
            !$OMP critical
            this%A(i,:) = this%A(i,:) + A_i
            if (body%asym_flow) then
                this%A(i+body%N_cp,:) = this%A(i+body%N_cp,:) + A_i_mir
            end if
            !$OMP end critical

        end do

        ! Clean up memory
        deallocate(this%wake_dod_info)

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_wake_influences


    subroutine panel_solver_solve_system(this, body)
        ! Solves the linear system for the singularity strengths

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        real,dimension(:,:),allocatable :: A_p
        real,dimension(:),allocatable :: b_p
        integer :: stat, i, j, unit
        real,dimension(:),allocatable :: A_ii_inv

        ! Set b vector for Morino formulation
        if (this%formulation == "morino") then
            this%b = -body%phi_cp_sigma
        end if

        if (run_checks) then

            if (verbose) write(*,'(a)',advance='no') "     Checking validity of linear system..."

            ! Check for NaNs
            if (any(isnan(this%A))) then
                write(*,*) "!!! Invalid value detected in A matrix. Quitting..."
                stop
            end if
            if (any(isnan(this%b))) then
                write(*,*) "!!! Invalid value detected in b vector. Quitting..."
                stop
            end if

            ! Check for uninfluenced/ing points
            do i=1,this%N
                if (all(this%A(i,:) == 0.)) then
                    write(*,*) "!!! Control point ", i, " is not influenced. Quitting..."
                    stop
                end if
                if (all(this%A(:,i) == 0.)) then
                    write(*,*) "!!! Vertex ", i, " exerts no influence. Quitting..."
                    stop
                end if
            end do

            if (verbose) write(*,*) "Done."

        end if

        ! Write A and b to file
        if (this%write_A_and_b) then

            if (verbose) write(*,'(a)',advance='no') "     Writing linear system to file..."

            ! A
            open(newunit=unit, file="A_mat.txt")
            do i=1,this%N
                write(unit,*) this%A(i,:)
            end do
            close(unit)

            ! b
            open(newunit=unit, file="b_vec.txt")
            do i=1,this%N
                write(unit,*) this%b(i)
            end do
            close(unit)

            if (verbose) write(*,*) "Done."

        end if

        if (verbose) write(*,'(a)',advance='no') "     Solving linear system..."

        ! Precondition
        select case(this%preconditioner)

        ! Diagonal preconditioning
        case ('DIAG')
            call diagonal_preconditioner(this%N, this%A, this%b, A_p, b_p)

        ! No preconditioning
        case default
            allocate(A_p, source=this%A, stat=stat)
            call check_allocation(stat, "solver copy of AIC matrix")
            allocate(b_p, source=this%b, stat=stat)
            call check_allocation(stat, "solver copy of b vector")

        end select

        ! Solve
        select case(this%matrix_solver)

        ! LU decomposition
        case ('LU')
            call lu_solve(this%N, A_p, b_p, body%mu)

        ! Purcell's method
        case ('PURC')
            call purcell_solve(this%N, A_p, b_p, body%mu)

        ! Block successive over-relaxation
        case ('BSOR')
            call block_sor_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, verbose, body%mu)

        ! Adaptive block SOR
        case ('ABSOR')
            this%rel = -1.
            call block_sor_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, verbose, body%mu)
        
        ! Block Jacobi
        case ('BJAC')
            call block_jacobi_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, verbose, body%mu)

        ! Optimally relaxed block Jacobi
        case ('ORBJ')
            this%rel = -1.
            call block_jacobi_solve(this%N, A_p, b_p, this%block_size, this%tol, this%rel, &
                                 this%max_iterations, verbose, body%mu)
        ! Improper specification
        case default
            write(*,*) "!!! ", this%matrix_solver, " is not a valid option. Defaulting to LU decomposition."
            call lu_solve(this%N, A_p, b_p, body%mu)

        end select
        if (verbose) write(*,*) "Done."

        ! Clean up memory
        deallocate(A_p)
        deallocate(b_p)

        ! Calculate potential at control points
        body%phi_cp_mu = matmul(this%A, body%mu)
        deallocate(this%A)
        body%phi_cp = body%phi_cp_mu+body%phi_cp_sigma

        ! Calculate residual parameters
        this%max_res = maxval(abs(body%phi_cp_mu-this%b))
        this%norm_res = sqrt(sum((body%phi_cp_mu-this%b)**2))
        deallocate(this%b)
        if (verbose) then
            write(*,*) "        Maximum residual:", this%max_res
            write(*,*) "        Norm of residual:", this%norm_res
        end if

    end subroutine panel_solver_solve_system


    subroutine panel_solver_calc_cell_velocities(this, body)
        ! Calculates the surface velocities on the mesh cells

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat, N_neighbors, j, i_panel, N_cycle
        real,dimension(3) :: v_jump, loc_mir, cent
        real,dimension(:,:),allocatable :: Vs

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface velocities..."

        ! Determine number of surface cells
        this%N_cells = body%N_panels
        N_cycle = 1
        if (body%asym_flow) then
            this%N_cells = this%N_cells*2
        end if
        if (doublet_order == 2) then
            this%N_cells = this%N_cells*4
            N_cycle = 4
        end if

        ! Allocate cell velocity storage
        allocate(body%V_cells(3,this%N_cells), stat=stat)
        call check_allocation(stat, "surface velocity vectors")

        ! Calculate the surface velocity on each existing panel
        !$OMP parallel do private(Vs, N_neighbors, v_jump, j, i_panel, loc_mir, cent) schedule(static)
        do i=1,body%N_panels

            ! Get velocity jump at centroid
            v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false.)

            ! Calculate velocity
            body%V_cells(:,(i-1)*N_cycle+1) = this%freestream%U*(this%inner_flow + v_jump)

            ! Calculate for subpanels
            if (doublet_order == 2) then
                do j=1,body%panels(i)%N

                    ! Get centroid
                    cent = body%panels(i)%get_subpanel_centroid(j)

                    ! Calculate velocity
                    v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .false., cent)
                    body%V_cells(:,(i-1)*N_cycle+1+j) = this%freestream%U*(this%inner_flow + v_jump)

                end do
            end if

            ! Calculate surface velocity on each mirrored panel
            if (body%asym_flow) then

                ! Get velocity jump
                v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .true.)

                ! Calculate velocity
                body%V_cells(:,(i-1)*N_cycle+1+body%N_panels) = this%freestream%U*(this%inner_flow + v_jump)

                ! Calculate for subpanels
                if (doublet_order == 2) then
                    do j=1,body%panels(i)%N

                        ! Get centroid
                        cent = mirror_across_plane(body%panels(i)%get_subpanel_centroid(j), body%mirror_plane)

                        ! Calculate velocity
                        v_jump = body%panels(i)%get_velocity_jump(body%mu, body%sigma, .true., cent)
                        body%V_cells(:,(i-1)*N_cycle+1+j+body%N_panels) = this%freestream%U*(this%inner_flow + v_jump)

                    end do
                end if

            end if
        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_cell_velocities


    subroutine panel_solver_calc_point_velocities(this, body)
        ! Calculates the surface velocities the mesh vertices

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat, N_neighbors, j, i_panel, N_cycle
        real,dimension(3) :: v_jump, loc_mir, cent
        real,dimension(:,:),allocatable :: Vs

        if (verbose) write(*,'(a)',advance='no') "     Calculating point velocities..."

        ! Allocate vertex velocity storage
        allocate(body%V_verts_avg(3,this%N))
        allocate(body%V_verts_std(3,this%N), source=0.)

        ! Calculate total potential on outside of mesh and quadratic-doublet velocity discrepancies
        !$OMP parallel do private(Vs, N_neighbors, v_jump, j, i_panel, loc_mir, cent) schedule(dynamic)
        do i=1,body%N_verts

            ! Allocate storage
            N_neighbors = body%vertices(i)%panels_not_across_wake_edge%len()
            allocate(Vs(3,N_neighbors))

            ! Get velocities
            do j=1,N_neighbors

                ! Get panel index
                call body%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                ! Get velocity jump at the vertex location
                v_jump = body%panels(i_panel)%get_velocity_jump(body%mu, body%sigma, .false., body%vertices(i)%loc)

                ! Store
                Vs(:,j) = this%freestream%U*(this%inner_flow + v_jump)
            end do

            ! Calculate stats
            body%V_verts_avg(:,i) = sum(Vs, dim=2)/N_neighbors
            if (N_neighbors > 1) then
                do j=1,N_neighbors
                    body%V_verts_std(:,i) = body%V_verts_std(:,i) + (Vs(:,j) - body%V_verts_avg(:,i))**2
                end do
                body%V_verts_std(:,i) = body%V_verts_std(:,i)/(N_neighbors - 1)
            end if

            ! Mirrored vertices
            if (body%asym_flow) then

                ! Get mirrored location
                loc_mir = mirror_across_plane(body%vertices(i)%loc, body%mirror_plane)

                ! Get velocities
                do j=1,N_neighbors

                    ! Get panel index
                    call body%vertices(i)%panels_not_across_wake_edge%get(j, i_panel)

                    ! Get velocity jump at the vertex location
                    v_jump = body%panels(i_panel)%get_velocity_jump(body%mu, body%sigma, .true., loc_mir)

                    ! Store
                    Vs(:,j) = this%freestream%U*(this%inner_flow + v_jump)
                end do

                ! Calculate stats
                body%V_verts_avg(:,i+body%N_verts) = sum(Vs, dim=2)/N_neighbors
                if (N_neighbors > 1) then
                    do j=1,N_neighbors
                        body%V_verts_std(:,i+body%N_verts) = body%V_verts_std(:,i+body%N_verts) &
                                                             + (Vs(:,j) - body%V_verts_avg(:,i))**2
                    end do
                    body%V_verts_std(:,i+body%N_verts) = body%V_verts_std(:,i+body%N_verts)/(N_neighbors - 1)
                end if

            end if

            deallocate(Vs)

        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_point_velocities


    subroutine panel_solver_calc_surface_potentials(this, body)
        ! Calculates the total outer potential at mesh vertices

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: loc_mir

        if (verbose) write(*,'(a)',advance='no') "     Calculating upper surface potentials..."

        ! Allocate storage
        allocate(body%Phi_u, source=body%mu, stat=stat)
        call check_allocation(stat, "total outer surface potential")

        ! Calculate total potential on outside of mesh and quadratic-doublet velocity discrepancies
        !$OMP parallel do private(loc_mir) schedule(dynamic)
        do i=1,body%N_verts

            ! Existing points
            body%Phi_u(i) = body%Phi_u(i) + inner(this%inner_flow, body%vertices(i)%loc)

            ! Mirrored points
            if (body%asym_flow) then
                loc_mir = mirror_across_plane(body%vertices(i)%loc, body%mirror_plane)
                body%Phi_u(i+body%N_verts) = body%phi_u(i+body%N_verts) + inner(this%inner_flow, loc_mir)
            end if
        end do

        if (verbose) write(*,*) "Done."

    end subroutine panel_solver_calc_surface_potentials


    subroutine panel_solver_calc_pressures(this, body)
        ! Calculates the surface pressures

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body

        integer :: i, stat
        real,dimension(3) :: V_pert
        real :: a, b, c, C_p_vac, lin, sln

        if (verbose) write(*,'(a)',advance='no') "     Calculating surface pressures..."

        ! Calculate vacuum pressure coefficient
        C_p_vac = -2./(this%freestream%gamma*this%freestream%M_inf**2)

        ! Allocate storage
        if (this%incompressible_rule) then
            allocate(body%C_p_inc(this%N_cells), stat=stat)
            call check_allocation(stat, "incompressible surface pressures")
        end if
        
        if (this%isentropic_rule) then
            allocate(body%C_p_ise(this%N_cells), stat=stat)
            call check_allocation(stat, "isentropic surface pressures")
        end if
        
        if (this%second_order_rule) then
            allocate(body%C_p_2nd(this%N_cells), stat=stat)
            call check_allocation(stat, "second-order surface pressures")
        end if
        
        if (this%slender_rule) then
            allocate(body%C_p_sln(this%N_cells), stat=stat)
            call check_allocation(stat, "slender-body surface pressures")
        end if
        
        if (this%linear_rule) then
            allocate(body%C_p_lin(this%N_cells), stat=stat)
            call check_allocation(stat, "linear surface pressures")
        end if

        ! Calculate reusable terms for the isentropic rule
        if (this%isentropic_rule) then
            a = 2./(this%freestream%gamma*this%freestream%M_inf**2)
            b = 0.5*(this%freestream%gamma-1.)*this%freestream%M_inf**2
            c = this%freestream%gamma/(this%freestream%gamma-1.)
        end if

        ! Calculate pressures
        !$OMP parallel do private(V_pert, lin, sln) schedule(static)
        do i=1,this%N_cells

            ! Incompressible rule
            if (this%incompressible_rule) then
                body%C_p_inc(i) = 1.-inner(body%V_cells(:,i), body%V_cells(:,i))*this%freestream%U_inv**2
            end if
        
            ! Isentropic rule
            if (this%isentropic_rule) then
                
                body%C_p_ise(i) = 1. - inner(body%V_cells(:,i), body%V_cells(:,i))*this%freestream%U_inv**2
                body%C_p_ise(i) = a*( (1. + b*body%C_p_ise(i))**c - 1.)

                ! Check for NaN and replace with vacuum pressure
                if (isnan(body%C_p_ise(i))) then
                    body%C_p_ise(i) = C_p_vac
                end if
            end if

            ! Get perturbation velocity in the compressible frame
            V_pert = matmul(this%freestream%A_g_to_c, body%V_cells(:,i)-this%freestream%v_inf)

            ! Linear term
            if (this%linear_rule .or. this%slender_rule .or. this%second_order_rule) then
                lin = -2.*V_pert(1)*this%freestream%U_inv
            end if

            ! Slender-body term
            if (this%slender_rule .or. this%second_order_rule) then
                sln = lin - (V_pert(2)**2 + V_pert(3)**2)*this%freestream%U_inv**2
            end if
        
            ! Second-order rule
            if (this%second_order_rule) then
                body%C_p_2nd(i) = sln - (1.-this%freestream%M_inf**2)*V_pert(1)**2*this%freestream%U_inv**2
            end if
        
            ! Slender-body rule
            if (this%slender_rule) then
                body%C_p_sln(i) = sln
            end if
        
            ! Linear rule
            if (this%linear_rule) then
                body%C_p_lin(i) = lin
            end if

        end do

        if (verbose) write(*,*) "Done."
        
        ! Report min and max pressure coefficients
        if (this%incompressible_rule .and. verbose) then
            write(*,*) "        Maximum incompressible pressure coefficient:", maxval(body%C_p_inc)
            write(*,*) "        Minimum incompressible pressure coefficient:", minval(body%C_p_inc)
        end if
        
        if (this%isentropic_rule .and. verbose) then
            write(*,*) "        Maximum isentropic pressure coefficient:", maxval(body%C_p_ise)
            write(*,*) "        Minimum isentropic pressure coefficient:", minval(body%C_p_ise)
        end if
        
        if (this%second_order_rule .and. verbose) then
            write(*,*) "        Maximum second-order pressure coefficient:", maxval(body%C_p_2nd)
            write(*,*) "        Minimum second-order pressure coefficient:", minval(body%C_p_2nd)
        end if
        
        if (this%slender_rule .and. verbose) then
            write(*,*) "        Maximum slender-body pressure coefficient:", maxval(body%C_p_sln)
            write(*,*) "        Minimum slender-body pressure coefficient:", minval(body%C_p_sln)
        end if
        
        if (this%linear_rule .and. verbose) then
            write(*,*) "        Maximum linear pressure coefficient:", maxval(body%C_p_lin)
            write(*,*) "        Minimum linear pressure coefficient:", minval(body%C_p_lin)
        end if
        
        ! Report vacuum pressure coefficient
        if (this%freestream%M_inf > 0.) then
            if (verbose) write(*,*) "        Vacuum pressure coefficient:", C_p_vac
        end if
        
        ! Apply subsonic pressure corrections
        if ((this%prandtl_glauert) .or. (this%karman_tsien) .or. (this%laitone)) then
            call this%subsonic_pressure_correction(body)
        end if

        ! Check if critical Mach number has been exceeded
        if (this%incompressible_rule) then
            if ((this%corrected_M_inf < 1) .and. (this%corrected_M_inf /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        else
            if ((this%freestream%M_inf < 1) .and. (this%freestream%M_inf /= 0.)) then
                call this%calc_crit_mach(body)
            end if
        end if
        
    end subroutine panel_solver_calc_pressures

    
    subroutine panel_solver_subsonic_pressure_correction(this, body)
        ! Apply selected method of correcting subsonic pressures
        
        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer :: i,stat
        real :: val_holder_L, val_holder_KT

        if (verbose) write(*,'(a)',advance='no') "     Calculating compressibility pressure corrections..."        
        
        ! Prandtl-Glauert rule
        if (this%prandtl_glauert) then
            ! Allocate storage
            allocate(body%C_p_pg(this%N_cells), stat=stat)
            call check_allocation(stat, "Prandtl-Glauert corrected surface pressures")
            
            ! Perform calculations for Prandtl-Glauert Rule (Modern Compressible Flow by John Anderson EQ 9.36)
            body%C_p_pg = body%C_p_inc / (sqrt(1 - (this%corrected_M_inf**2)))
        end if
        
        ! Laitone rule
        if (this%laitone) then
            ! Allocate storage
            allocate(body%C_p_lai(this%N_cells), stat=stat)
            call check_allocation(stat, "Laitone corrected surface pressures")
            
            ! Perform calculations (Modern Compressible Flow by John Anderson EQ 9.39)
            val_holder_L = this%corrected_M_inf**2 * (1 + (0.5 * (this%freestream%gamma - 1) * this%corrected_M_inf**2)) &
            / (2 * sqrt(1 - this%corrected_M_inf**2))
            
            body%C_p_lai = body%C_p_inc / &
            (sqrt(1 - this%corrected_M_inf**2) + (val_holder_L * body%C_p_inc))
        end if
        
        ! Karman-Tsien rule
        if (this%karman_tsien) then
            ! Allocate storage
            allocate(body%C_p_kt(this%N_cells), stat=stat)
            call check_allocation(stat, "Karman-Tsien corrected surface pressures")
            
            ! Perform calculations (Modern Compressible Flow by John Anderson EQ 9.40)
            val_holder_KT = this%corrected_M_inf**2 / (1 + sqrt(1 - this%corrected_M_inf**2))
            
            body%C_p_kt = body%C_p_inc / &
            (sqrt(1 - this%corrected_M_inf**2) + val_holder_KT * (0.5 * body%C_p_inc))
        end if

        if (verbose) write(*,*) "Done. "  
        
    end subroutine panel_solver_subsonic_pressure_correction
    

    subroutine panel_solver_calc_crit_mach(this, body)
        ! Warn the user if body is in transonic flow

        implicit none

        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        integer,dimension(1) :: min_loc
        real :: C_p_crit, C_p_min, numerator, denominator, M_inf_selected

        ! Locate minimum pressure location on body based on pressure for forces selection
        select case (this%pressure_for_forces)
            case ("incompressible")
                min_loc = MINLOC(body%C_p_inc)
                C_p_min = MINVAL(body%C_p_inc)
                M_inf_selected = this%corrected_M_inf  
            
            case ("prandtl-glauert")
                min_loc = MINLOC(body%C_p_pg)
                C_p_min = MINVAL(body%C_p_pg)
                M_inf_selected = this%corrected_M_inf

            case ("karman-tsien")
                min_loc = MINLOC(body%C_p_kt)
                C_p_min = MINVAL(body%C_p_kt)
                M_inf_selected = this%corrected_M_inf

            case ("laitone")
                min_loc = MINLOC(body%C_p_lai)
                C_p_min = MINVAL(body%C_p_lai)
                M_inf_selected = this%corrected_M_inf

            case ("isentropic")
                min_loc = MINLOC(body%C_p_ise)
                C_p_min = MINVAL(body%C_p_ise)
                M_inf_selected = this%freestream%M_inf

            case ("second-order")
                min_loc = MINLOC(body%C_p_2nd)
                C_p_min = MINVAL(body%C_p_2nd)
                M_inf_selected = this%freestream%M_inf

            case ("linear")
                min_loc = MINLOC(body%C_p_lin)
                C_p_min = MINVAL(body%C_p_lin)
                M_inf_selected = this%freestream%M_inf

            case ("slender-body")
                min_loc = MINLOC(body%C_p_sln)
                C_p_min = MINVAL(body%C_p_sln)
                M_inf_selected = this%freestream%M_inf
                
            case default
                write(*,*) "!!! ", this%pressure_for_forces," pressure for forces is not available. Quitting..."
                stop
        end select

        ! Calculate critical pressure coefficient for the selected Mach number (Modern Compressible Flow by John Anderson EQ 9.55)
        numerator = 1 + ((this%freestream%gamma - 1) / 2) * M_inf_selected**2 
        denominator = 1 + (this%freestream%gamma -1) / 2

        C_p_crit = (2 / (this%freestream%gamma * M_inf_selected**2)) &
                    * (((numerator/denominator) ** (this%freestream%gamma / (this%freestream%gamma - 1))) - 1)
    
        ! Report the critical mach pressure coefficient and the minimum pressure coefficients to the terminal
        if (verbose) then
        write(*,*) "        Critical Mach C_p = ", C_p_crit
        write(*,*) "        Minimum C_p       = ", C_p_min
        end if

        ! Throw warning message indicating body is in transonic flow
        if ((C_p_min <= C_p_crit) .and. (verbose)) then
            write(*,*) "!!! The critical Mach number has been exceeded over panel ", min_loc, ". Results may not be reliable."
            
        end if

    end subroutine panel_solver_calc_crit_mach
    

    subroutine panel_solver_calc_forces(this, body)
        ! Calculates the forces
        
        implicit none
        
        class(panel_solver),intent(inout) :: this
        type(surface_mesh),intent(inout) :: body
        
        integer :: i, stat

        if (verbose) write(*,'(a, a, a)',advance='no') "     Calculating forces using the ", & 
                                                       this%pressure_for_forces, " pressure rule..."

        ! Allocate force storage
        allocate(body%dC_f(3,this%N_cells), stat=stat)
        call check_allocation(stat, "forces")

        ! Calculate total forces
        !$OMP parallel do schedule(static)
        do i=1,body%N_panels

            select case (this%pressure_for_forces)

            case ('incompressible')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = -body%C_p_inc(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = -body%C_p_inc(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            case ('isentropic')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = -body%C_p_ise(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = -body%C_p_ise(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            case ('second-order')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = -body%C_p_2nd(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = -body%C_p_2nd(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            case ('prandtl-glauert')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = body%C_p_pg(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = body%C_p_pg(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            case ('karman-tsien')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = body%C_p_kt(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = body%C_p_kt(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            case ('laitone')

                ! Discrete force coefficient acting on panel
                body%dC_f(:,i) = body%C_p_lai(i)*body%panels(i)%A*body%panels(i)%n_g

                ! Mirror
                if (body%mirrored .and. body%asym_flow) then
                    body%dC_f(:,i+body%N_panels) = body%C_p_lai(i+body%N_panels)*body%panels(i)%A*body%panels(i)%n_g_mir
                end if

            end select

        end do

        ! Sum discrete forces
        this%C_F(:) = sum(body%dC_f, dim=2)/body%S_ref

        ! Add contributions from mirrored half
        if (body%mirrored .and. .not. body%asym_flow) then
            this%C_F = 2.*this%C_F
            this%C_F(body%mirror_plane) = 0. ! We know this
        end if

        if (verbose) then
            write(*,*) "Done."
            write(*,*) "        Cx:", this%C_F(1)
            write(*,*) "        Cy:", this%C_F(2)
            write(*,*) "        Cz:", this%C_F(3)
        end if
    
    end subroutine panel_solver_calc_forces


    subroutine panel_solver_update_report(this, p_json, body)
        ! Updates the report JSON with the information relevant to the solver

        implicit none

        class(panel_solver),intent(in) :: this
        type(json_value),pointer,intent(inout) :: p_json
        type(surface_mesh),intent(in) :: body

        type(json_value),pointer :: p_parent, p_child
        integer :: i_unit

        ! Write solver results
        call json_value_create(p_parent)
        call to_object(p_parent, 'solver_results')
        call json_value_add(p_json, p_parent)
        call json_value_create(p_child)
        call to_object(p_child, 'residual')
        call json_value_add(p_parent, p_child)
        call json_value_add(p_child, 'max', this%max_res)
        call json_value_add(p_child, 'norm', this%norm_res)
        nullify(p_parent)
        nullify(p_child)

        ! Write pressure results
        call json_value_create(p_parent)
        call to_object(p_parent, 'pressure_calculations')
        call json_value_add(p_json, p_parent)

        ! Incompressible rule
        if (this%incompressible_rule) then
            call json_value_create(p_child)
            call to_object(p_child, 'incompressible_rule')
            call json_value_add(p_parent, p_child)
            call json_value_add(p_child, 'max', maxval(body%C_p_inc))
            call json_value_add(p_child, 'min', minval(body%C_p_inc))
            nullify(p_child)
        end if

        ! Isentropic rule
        if (this%isentropic_rule) then
            call json_value_create(p_child)
            call to_object(p_child, 'isentropic_rule')
            call json_value_add(p_parent, p_child)
            call json_value_add(p_child, 'max', maxval(body%C_p_ise))
            call json_value_add(p_child, 'min', minval(body%C_p_ise))
            nullify(p_child)
        end if

        ! Second-order rule
        if (this%second_order_rule) then
            call json_value_create(p_child)
            call to_object(p_child, 'second_order_rule')
            call json_value_add(p_parent, p_child)
            call json_value_add(p_child, 'max', maxval(body%C_p_2nd))
            call json_value_add(p_child, 'min', minval(body%C_p_2nd))
            nullify(p_child)
        end if
        nullify(p_parent)

        ! Write forces
        call json_value_create(p_parent)
        call to_object(p_parent, 'total_forces')
        call json_value_add(p_json, p_parent)
        call json_value_add(p_parent, 'Cx', this%C_F(1))
        call json_value_add(p_parent, 'Cy', this%C_F(2))
        call json_value_add(p_parent, 'Cz', this%C_F(3))
        nullify(p_parent)

    end subroutine panel_solver_update_report


    subroutine panel_solver_get_panel_induced_potentials(this, i_panel, source_inf, doublet_inf, mirror_panel, body, phi_d, phi_s)
        ! Calculates the induced potential from the panel

        implicit none

        class(panel_solver),intent(in) :: this
        integer,intent(in) :: i_panel
        real,dimension(:),allocatable,intent(in) :: source_inf, doublet_inf
        logical,intent(in) :: mirror_panel
        type(surface_mesh),intent(in) :: body
        real,intent(out) :: phi_d, phi_s

        integer :: k

        phi_d = 0.
        phi_s = 0.

        ! Add source influence (if sources are present)
        if (this%morino) then

            ! Constant
            if (source_order == 0) then
                if (mirror_panel) then
                    phi_s = source_inf(1)*body%sigma(i_panel + body%N_panels)
                else
                    phi_s = source_inf(1)*body%sigma(i_panel)
                end if

            ! Linear
            else
                do k=1,size(body%panels(i_panel)%i_vert_s)
                    if (mirror_panel) then
                        phi_s = phi_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k) + body%N_cp)
                    else
                        phi_s = phi_s + source_inf(k)*body%sigma(body%panels(i_panel)%i_vert_s(k))
                    end if
                end do
            end if
        end if

        ! Add doublet influence
        ! This method is the same for linear and quadratic doublets
        ! Loop through panel vertices
        do k=1,size(body%panels(i_panel)%i_vert_d)
            if (mirror_panel) then
                phi_d = phi_d + doublet_inf(k)*body%mu(body%panels(i_panel)%i_vert_d(k) + body%N_cp)
            else
                phi_d = phi_d + doublet_inf(k)*body%mu(body%panels(i_panel)%i_vert_d(k))
            end if
        end do
        
    end subroutine panel_solver_get_panel_induced_potentials


    subroutine panel_solver_export_potential_slice(this, slice_file, output_settings, body, freestream)
        ! Writes out a csv file of potentials on a user-specified slice of the flow

        implicit none

        class(panel_solver),intent(inout) :: this
        character(len=:),allocatable :: slice_file
        type(json_value),pointer,intent(in) :: output_settings
        type(surface_mesh),intent(inout) :: body
        type(flow),intent(inout) :: freestream

        integer :: i, j, k, unit, i_fixed, i_var_1, i_var_2, N_var_1, N_var_2, N_points, stat
        character(len=:),allocatable :: parallel_plane
        real,dimension(2) :: var_lims_1, var_lims_2
        real,dimension(3) :: vert_loc, mirrored_vert_loc
        real :: fixed_dim, d1, d2, var_1_min, var_2_min, var_1_max, var_2_max, phi_inf, phi_d, phi_s, phi_d_tot, phi_s_tot
        real,dimension(:,:,:),allocatable :: points
        real,dimension(:),allocatable :: var_1_space, var_2_space, source_inf, doublet_inf
        type(dod),dimension(:,:,:),allocatable :: dod_info
        logical,dimension(:),allocatable :: verts_in_dod

        if (verbose) write(*,'(a)',advance='no') "    Calculating potential on user-specified slice..."

        ! Get settings
        call json_xtnsn_get(output_settings, 'potential_slice.parallel_to', parallel_plane, 'xz')        

        select case(parallel_plane)

        case ('xz')

            ! Specify which axes are fixed and which are variable
            i_fixed = 2
            i_var_1 = 1
            i_var_2 = 3

            ! Get location of slice
            call json_xtnsn_get(output_settings, 'potential_slice.y_loc', fixed_dim, 0.0)
            call json_xtnsn_get(output_settings, 'potential_slice.x_min', var_lims_1(1))
            call json_xtnsn_get(output_settings, 'potential_slice.z_min', var_lims_2(1))
            call json_xtnsn_get(output_settings, 'potential_slice.x_max', var_lims_1(2))
            call json_xtnsn_get(output_settings, 'potential_slice.z_max', var_lims_2(2))

            ! Get number of points in the slice
            call json_xtnsn_get(output_settings, 'potential_slice.Nx', N_var_1)
            call json_xtnsn_get(output_settings, 'potential_slice.Nx', N_var_2)

        end select

        ! Initialize ranges of variable locations
        allocate(var_1_space(N_var_1))
        d1 = (var_lims_1(2) - var_lims_1(1))/(N_var_1-1)
        do i=1,N_var_1
            var_1_space(i) = var_lims_1(1) + (i-1)*d1
        end do

        allocate(var_2_space(N_var_2))
        d2 = (var_lims_2(2) - var_lims_2(1))/(N_var_2-1)
        do i=1,N_var_2
            var_2_space(i) = var_lims_2(1) + (i-1)*d2
        end do

        ! Initialize points
        N_points = N_var_1*N_var_2
        allocate(points(3, N_var_2, N_var_1))
        points(i_fixed,:,:) = fixed_dim
        do i=1,N_var_1
            do j=1,N_var_2
                points(i_var_1,j,i) = var_1_space(i)
                points(i_var_2,j,i) = var_2_space(j)
            end do
        end do
        deallocate(var_1_space)
        deallocate(var_2_space)

        ! Allocate DoD storage
        if (body%mirrored) then
            allocate(dod_info(2*body%N_panels,N_var_2,N_var_1), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(2*body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        else
            allocate(dod_info(body%N_panels,N_var_2,N_var_1), stat=stat)
            call check_allocation(stat, "domain of dependence storage")

            allocate(verts_in_dod(body%N_verts), stat=stat)
            call check_allocation(stat, "vertex domain of dependence storage")
        end if

        ! If the freestream is supersonic, calculate domain of dependence info
        if (this%freestream%supersonic) then

            ! Loop through control points
            !$OMP parallel do private(i, vert_loc, mirrored_vert_loc, verts_in_dod)
            do i=1,N_var_1
                do j=1,N_var_2

                    ! Initialize for this point
                    verts_in_dod = .false.

                    ! Loop through body vertices
                    do k=1,body%N_verts

                        ! Get this vertex location
                        vert_loc = body%vertices(k)%loc

                        ! Check if this vertex is in the DoD
                        verts_in_dod(k) = this%freestream%point_in_dod(vert_loc, points(:,j,i))

                        if (body%mirrored) then

                            ! Get the mirrored vertex location
                            mirrored_vert_loc = mirror_across_plane(vert_loc, body%mirror_plane)

                            ! Check if this vertex is in the DoD
                            verts_in_dod(k+body%N_verts) = this%freestream%point_in_dod(mirrored_vert_loc, &
                                                                                        points(:,j,i))

                        end if
                    end do

                    ! Loop through body panels
                    do k=1,body%N_panels

                        ! Check if the panel is in the DoD
                        dod_info(k,j,i) = body%panels(k)%check_dod(points(:,j,i), this%freestream, verts_in_dod)

                        if (body%mirrored) then

                            ! Check DoD for mirrored panel
                            dod_info(k+body%N_panels,j,i) = body%panels(i)%check_dod(points(:,j,i), this%freestream, &
                                                                                     verts_in_dod, &
                                                                                     .true., body%mirror_plane)
                        end if
                    end do

                end do
            end do

        end if
        deallocate(verts_in_dod)

        ! Open slice file
        100 format(e20.12, ', ', e20.12, ', ', e20.12, ', ', e20.12, ', ', e20.12, ', ', e20.12, ', ', e20.12, ', ', e20.12)
        open(newunit=unit, file=slice_file)

        ! Write header
        write(unit,*) 'x,y,z,phi_inf,phi_d,phi_s,phi,Phi'

        ! Calculate potentials
        do i=1,N_var_1
            do j=1,N_var_2

                ! Calculate freestream potential
                phi_inf = this%freestream%U*inner(points(:,j,i), freestream%c_hat_g)

                ! Loop through panels
                phi_s_tot = 0.
                phi_d_tot = 0.
                do k=1,body%N_panels

                    ! Check DoD
                    if (dod_info(k,j,i)%in_dod) then
                    
                        ! Calculate influence
                        call body%panels(k)%calc_potentials(points(:,j,i), this%freestream, dod_info(k,j,i), .false., &
                                                            source_inf, doublet_inf)

                        ! Add influence
                        call this%get_panel_induced_potentials(k, source_inf, doublet_inf, .false., body, phi_d, phi_s)
                        phi_d_tot = phi_d_tot + phi_d
                        phi_s_tot = phi_s_tot + phi_s

                    end if

                    ! Calculate mirrored influences
                    if (body%mirrored) then

                        if (dod_info(k+body%N_panels,j,i)%in_dod) then

                            if (body%asym_flow) then

                                ! Get influence of mirrored panel
                                call body%panels(k)%calc_potentials(points(:,j,i), this%freestream, &
                                                                    dod_info(k+body%N_panels,j,i), .true., &
                                                                    source_inf, doublet_inf)

                                ! Add influence of mirrored panel using asymmetrical strengths
                                call this%get_panel_induced_potentials(k, source_inf, doublet_inf, .true., body, phi_d, phi_s)
                                phi_d_tot = phi_d_tot + phi_d
                                phi_s_tot = phi_s_tot + phi_s

                            else

                                ! Get influence of panel on mirrored point
                                call body%panels(k)%calc_potentials(mirror_across_plane(points(:,j,i), body%mirror_plane), &
                                                                    this%freestream, &
                                                                    dod_info(k+body%N_panels,j,i), .false., &
                                                                    source_inf, doublet_inf)

                                ! Add influence of mirrored panel using symmetrical strengths
                                call this%get_panel_induced_potentials(k, source_inf, doublet_inf, .false., body, phi_d, phi_s)
                                phi_d_tot = phi_d_tot + phi_d
                                phi_s_tot = phi_s_tot + phi_s

                            end if

                        end if

                    end if

                end do

                ! Write to file
                phi_d_tot = phi_d_tot*this%freestream%U
                phi_s_tot = phi_s_tot*this%freestream%U
                write(unit,100) points(1,j,i), points(2,j,i), points(3,j,i), phi_inf, phi_d_tot, phi_s_tot, phi_d_tot + phi_s_tot, &
                                phi_inf + phi_d_tot + phi_s_tot

            end do
        end do

        close(unit)

        write(*,*) "Done."


    end subroutine panel_solver_export_potential_slice



end module panel_solver_mod