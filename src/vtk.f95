! Subroutines for I/O with VTK files
module vtk_mod

    use panel_mod
    use vertex_mod

    implicit none

    
contains

    subroutine load_surface_vtk(mesh_file, N_verts, N_panels, vertices, panels)

        implicit none

        character(len=:),allocatable,intent(in) :: mesh_file
        integer,intent(out) :: N_verts, N_panels
        character(len=200) :: dummy_read
        type(vertex),dimension(:),allocatable,intent(inout) :: vertices
        type(panel),dimension(:),allocatable,intent(inout) :: panels
        integer :: i, j, N, i1, i2, i3, i4

        ! Open file
        open(1, file=mesh_file)

            ! Determine number of vertices
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) ! Header
            read(1,*) dummy_read, N_verts, dummy_read

            ! Allocate vertex array
            allocate(vertices(N_verts))

            ! Store vertices
            do i=1,N_verts
                read(1,*) vertices(i)%loc(1), vertices(i)%loc(2), vertices(i)%loc(3)
                vertices(i)%index = i
            end do

            ! Determine number of panels
            read(1,*) dummy_read, N_panels, dummy_read

            ! Allocate panel array
            allocate(panels(N_panels))

            ! Initialize panels
            do i=1,N_panels

                ! Get data
                read(1,'(a)') dummy_read
                
                ! Determine size of panel
                if (dummy_read(1:2) == '3 ') then
                    read(dummy_read,*) N, i1, i2, i3
                else if (dummy_read(1:2) == '4 ') then
                    read(dummy_read,*) N, i1, i2, i3, i4
                else
                    write(*,*) "MFTran supports only triangular and quadrilateral panels."
                    stop
                end if

                ! Initialize triangular panel
                if (N == 3) then
                    call panels(i)%init(vertices(i1+1), vertices(i2+1), vertices(i3+1)) ! Need +1 because VTK uses 0-based indexing

                    ! Add panel index to vertices
                    call vertices(i1+1)%panels%append(i)
                    call vertices(i2+1)%panels%append(i)
                    call vertices(i3+1)%panels%append(i)

                ! Initialize quadrilateral panel
                else
                    call panels(i)%init(vertices(i1+1), vertices(i2+1), vertices(i3+1), vertices(i4+1))

                    ! Add panel index to vertices
                    call vertices(i1+1)%panels%append(i)
                    call vertices(i2+1)%panels%append(i)
                    call vertices(i3+1)%panels%append(i)
                    call vertices(i4+1)%panels%append(i)
                end if

            end do

        close(1)
    
    end subroutine load_surface_vtk


    subroutine write_surface_vtk(output_file, vertices, panels)

        implicit none

        character(len=:),allocatable,intent(in) :: output_file
        type(vertex),dimension(:),intent(in) :: vertices
        type(panel),dimension(:),intent(in) :: panels
        integer :: i, N_verts, N_panels, panel_info_size, j

        ! Open file
        open(1, file=output_file)

            ! Write header
            write(1,'(a)') "# vtk DataFile Version 3.0"
            write(1,'(a)') "MFTran results file. Generated by MFTran, USU AeroLab (c) 2021."
            write(1,'(a)') "ASCII"

            ! Write vertex information
            N_verts = size(vertices)
            write(1,'(a)') "DATASET POLYDATA"
            write(1,'(a i20 a)') "POINTS", N_verts, " float"

            ! Write out vertices
            100 format(f20.12, ' ', f20.12, ' ', f20.12) ! Vertices
            do i=1,N_verts
                write(1,100) vertices(i)%loc(1), vertices(i)%loc(2), vertices(i)%loc(3)
            end do

            ! Determine panel info size
            panel_info_size = 0
            N_panels = size(panels)
            do i=1,N_panels
                panel_info_size = panel_info_size + panels(i)%N + 1
            end do
            
            ! Write out panels
            write(1,'(a i20 i20)') "POLYGONS", N_panels, panel_info_size
            do i=1,N_panels

                ! Number of vertices
                write(1,'(i1) ',advance='no') panels(i)%N

                ! Indices of each vertex
                do j=1,panels(i)%N
                    write(1,'(i20) ',advance='no') panels(i)%vertices(j)%ptr%index-1
                end do
                write(1,*)

            end do

            ! Write out panel normals
            write(1,'(a i20)') "CELL_DATA", N_panels
            write(1,'(a)') "NORMALS panel_normals float"
            do i=1,N_panels
                write(1,100) panels(i)%normal(1), panels(i)%normal(2), panels(i)%normal(3)
            end do

            ! Panel source strengths
            write(1,'(a)') "SCALARS phi_n float 1"
            write(1,'(a)') "LOOKUP_TABLE default"
            do i=1,N_panels
                write(1,'(f20.12)') panels(i)%phi_n
            end do

            ! Vertex doublet strengths
            write(1, '(a i20)') "POINT_DATA", N_verts
            write(1,'(a)') "SCALARS phi float 1"
            write(1,'(a)') "LOOKUP_TABLE default"
            do i=1,N_verts
                write(1,'(f20.12)') vertices(i)%phi
            end do

        close(1)
        
    
    end subroutine write_surface_vtk

    
end module vtk_mod