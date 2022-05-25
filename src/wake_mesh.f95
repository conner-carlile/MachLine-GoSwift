! Class for modeling wake meshes
module wake_mesh_mod

    use json_mod
    use json_xtnsn_mod
    use linked_list_mod
    use helpers_mod
    use vertex_mod
    use panel_mod
    use math_mod
    use flow_mod
    use edge_mod

    implicit none


    type wake_mesh

        type(vertex),allocatable,dimension(:) :: vertices
        type(panel),allocatable,dimension(:) :: panels
        integer :: N_verts, N_panels

        contains

            procedure :: init => wake_mesh_init

    end type wake_mesh


contains


    subroutine wake_mesh_init(this, body_edges, body_verts, N_body_panels, freestream, asym_flow, mirror_plane, &
                              N_panels_streamwise, trefftz_dist)
        ! Creates the vertices and panels. Handles vertex association.

        implicit none

        class(wake_mesh),intent(inout) :: this
        type(edge),allocatable,dimension(:),intent(in) :: body_edges
        type(vertex),allocatable,dimension(:),intent(inout) :: body_verts
        integer,intent(in) :: N_body_panels
        type(flow),intent(in) :: freestream
        logical,intent(in) :: asym_flow
        integer,intent(in) :: mirror_plane, N_panels_streamwise
        real,intent(in) :: trefftz_dist

        real :: distance, vertex_separation, mirrored_distance, mirrored_vertex_separation
        real,dimension(3) :: loc, start, mirrored_start
        integer :: i, j, k, i_vert, i_mirrored_vert, i_panel, i_top_parent, i_bot_parent, i_start, i_stop, i1, i2, i3
        integer :: N_wake_edge_verts
        integer :: N_wake_edges
        integer,dimension(:),allocatable :: wake_edge_indices, wake_edge_verts
        logical,dimension(:),allocatable :: is_wake_edge_vertex

        if (verbose) write(*,'(a ES10.4 a)',advance='no') "     Initializing wake with a Trefftz distance of ", trefftz_dist, "..."

        ! Count up wake-shedding edges
        N_wake_edges = 0
        do i=1,size(body_edges)
            if (body_edges(i)%sheds_wake) N_wake_edges = N_wake_edges + 1
        end do

        ! Get indices of wake-shedding edges and vertices
        allocate(wake_edge_indices(N_wake_edges))
        allocate(is_wake_edge_vertex(size(body_verts)))
        j = 0
        do i=1,size(body_edges)
            if (body_edges(i)%sheds_wake) then

                ! Store edge index
                j = j + 1
                wake_edge_indices(j) = i

                ! Store that the vertices are wake-shedding
                is_wake_edge_vertex(body_edges(i)%verts(1)) = .true.
                is_wake_edge_vertex(body_edges(i)%verts(2)) = .true.

            end if
        end do

        ! Determine how many vertices are along the wake-shedding edges
        N_wake_edge_verts = count(is_wake_edge_vertex)

        ! Get wake-shedding edge vertex indices
        allocate(wake_edge_verts(N_wake_edge_verts))
        j = 0
        do i=1,size(body_verts)
            if (is_wake_edge_vertex(i)) then
                j = j + 1
                wake_edge_verts(j) = i
                body_verts(i)%index_in_wake_vertices = j
            end if
        end do

        ! Determine necessary number of vertices
        if (asym_flow) then
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)*2
        else
            this%N_verts = N_wake_edge_verts*(N_panels_streamwise+1)
        end if

        ! Allocate vertex storage
        allocate(this%vertices(this%N_verts))

        ! Determine vertex placement
        i_vert = 0
        do i=1,N_wake_edge_verts

            ! Check this is not a midpoint vertex
            i_top_parent = wake_edge_verts(i)
            if (body_verts(i_top_parent)%vert_type == 1) then

                ! Get indices
                i_bot_parent = body_verts(i_top_parent)%i_wake_partner

                ! Determine distance from origin to wake-shedding vertex in the direction of the freestream flow
                start = body_verts(i_top_parent)%loc
                distance = trefftz_dist - inner(start, freestream%c_hat_g)

                ! Double-check
                if (dist(start, body_verts(i_bot_parent)%loc) > 1.e-12) then
                    write(*,*) "!!! Wake edge vertices are not identical. Quitting..."
                    stop
                end if

                ! Determine vertex separation
                vertex_separation = distance/N_panels_streamwise

                ! Same for mirror
                if (asym_flow) then

                    ! Determine start location
                    mirrored_start = mirror_across_plane(start, mirror_plane)
                    mirrored_distance = trefftz_dist-inner(mirrored_start, freestream%c_hat_g)

                    ! Determine vertex separation
                    mirrored_vertex_separation = mirrored_distance/N_panels_streamwise

                end if

                ! Loop down the streamwise direction to place vertices
                do j=1,N_panels_streamwise+1

                    ! Determine location
                    i_vert = i_vert + 1
                    loc = start + vertex_separation*(j-1)*freestream%c_hat_g

                    ! Initialize vertex
                    call this%vertices(i_vert)%init(loc, i_vert, 1)

                    ! Set parent index
                    this%vertices(i_vert)%top_parent = i_top_parent
                    this%vertices(i_vert)%bot_parent = i_bot_parent

                    ! Initialize mirror
                    if (asym_flow) then

                        ! Determine location
                        i_mirrored_vert = i_vert + this%N_verts/2
                        mirrored_start = mirror_across_plane(start, mirror_plane)
                        loc = mirrored_start + mirrored_vertex_separation*(j-1)*freestream%c_hat_g

                        ! Initialize vertex
                        call this%vertices(i_mirrored_vert)%init(loc, i_mirrored_vert, 1)

                        ! Set parent index
                        this%vertices(i_mirrored_vert)%top_parent = i_top_parent + size(body_verts)
                        this%vertices(i_mirrored_vert)%bot_parent = i_bot_parent + size(body_verts)

                    end if

                end do
            end if
        end do

        ! Determine necessary number of panels
        if (asym_flow) then
            this%N_panels = N_wake_edges*N_panels_streamwise*4
        else
            this%N_panels = N_wake_edges*N_panels_streamwise*2
        end if
        allocate(this%panels(this%N_panels))

        ! Loop through wake-shedding edges to intialize panels
        do k=1,N_wake_edges

            ! Get index of edge
            i = wake_edge_indices(k)

            ! Determine which wake-shedding vertices this panel lies between
            i_start = body_verts(body_edges(i)%verts(1))%index_in_wake_vertices
            i_stop = body_verts(body_edges(i)%verts(2))%index_in_wake_vertices

            ! Create panels heading downstream
            do j=1,N_panels_streamwise

                ! Determine index of first triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j-1

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_start-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1

                ! Initialize
                call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)

                ! Specify this panel is in the wake
                this%panels(i_panel)%in_wake = .true.
                this%panels(i_panel)%top_parent = body_edges(i)%panels(1)
                this%panels(i_panel)%bot_parent = body_edges(i)%panels(2)

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_start-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2

                    ! Initialize (order of vertices is reversed to maintain panel orientation through mirror)
                    call this%panels(i_panel)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_panel)

                    ! Specify this panel is in the wake
                    this%panels(i_panel)%in_wake = .true.
                    this%panels(i_panel)%top_parent = body_edges(i)%panels(1)+N_body_panels
                    this%panels(i_panel)%bot_parent = body_edges(i)%panels(2)+N_body_panels

                end if

                ! Determine index of second triangular panel
                i_panel = (k-1)*N_panels_streamwise*2+2*j

                ! Determine vertex indices
                i1 = (i_start-1)*(N_panels_streamwise+1)+j
                i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1
                i3 = (i_stop-1)*(N_panels_streamwise+1)+j

                ! Initialize
                call this%panels(i_panel)%init(this%vertices(i1), this%vertices(i2), this%vertices(i3), i_panel)

                ! Specify this panel is in the wake
                this%panels(i_panel)%in_wake = .true.
                this%panels(i_panel)%top_parent = body_edges(i)%panels(1)
                this%panels(i_panel)%bot_parent = body_edges(i)%panels(2)

                ! Create mirror
                if (asym_flow) then

                    ! Determine index
                    i_panel = i_panel + this%N_panels/2

                    ! Determine vertex indices
                    i1 = (i_start-1)*(N_panels_streamwise+1)+j+this%N_verts/2
                    i2 = (i_stop-1)*(N_panels_streamwise+1)+j+1+this%N_verts/2
                    i3 = (i_stop-1)*(N_panels_streamwise+1)+j+this%N_verts/2

                    ! Initialize (again, order is reversed)
                    call this%panels(i_panel)%init(this%vertices(i3), this%vertices(i2), this%vertices(i1), i_panel)

                    ! Specify this panel is in the wake
                    this%panels(i_panel)%in_wake = .true.
                    this%panels(i_panel)%top_parent = body_edges(i)%panels(1)+N_body_panels
                    this%panels(i_panel)%bot_parent = body_edges(i)%panels(2)+N_body_panels

                end if

            end do
        end do

        ! Initialize freestream-dependent properties
        ! The mirror of wake panels will never need to be initialized
        do i=1,this%N_panels
            call this%panels(i)%init_with_flow(freestream, .false., mirror_plane)
        end do

        if (verbose) write(*,'(a, i7, a, i7, a)') "Done. Created ", this%N_verts, " wake vertices and ", &
                                                  this%N_panels, " wake panels."

    end subroutine wake_mesh_init


end module wake_mesh_mod