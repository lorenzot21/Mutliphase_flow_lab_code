subroutine generate_output(nstep)

use commondata

integer :: nstep,numx,numy,numz
integer :: i,k,j
integer :: x_start,x_end,dnx
integer :: y_start,y_end,dny
integer :: z_start,z_end,dnz
character(len=40) :: namefile
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

! end of line character
lf=achar(10)

! fields included
nfields=uflag + vflag + wflag + phiflag

!input??
x_start=1
x_end=nx
dnx=1
y_start=1
y_end=ny
dny=1
z_start=1
z_end=nz
dnz=1

numx=0
numy=0
numz=0
do i=x_start,x_end,dnx
 numx=numx+1
enddo
do j=y_start,y_end,dny
 numy=numy+1
enddo
do k=z_start,z_end,dnz
 numz=numz+1
enddo



write(namefile,'(a,i8.8,a)') './output/OUTPAR_',nstep,'.vtk'
open(66,file=trim(namefile),status='new',form='unformatted',access='stream',convert='big_endian')
! start writing vtk file
! write header
buffer='# vtk DataFile Version 3.0'//lf
write(66) trim(buffer)
buffer='Phase field'//lf
write(66) trim(buffer)
buffer='BINARY'//lf
write(66) trim(buffer)
buffer='DATASET RECTILINEAR_GRID'//lf
 write(66) trim(buffer)

 !write grid
 write(str1(1:8),'(i8)') numx
 write(str2(1:8),'(i8)') numy
 write(str3(1:8),'(i8)') numz
 buffer='DIMENSIONS '//str1//str2//str3//lf
 write(66) trim(buffer)
 buffer='X_COORDINATES '//str1//'  float'//lf ;
 write(66) trim(buffer)
 do i=x_start,x_end,dnx
  write(66) real(x(i))
 enddo
 buffer='Y_COORDINATES '//str2//'  float'//lf ;
 write(66) trim(buffer)
 do j=y_start,y_end,dny
  write(66) real(y(j))
 enddo
 buffer='Z_COORDINATES '//str3//'  float'//lf ;
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  write(66) real(z(k))
 enddo

 ! write content (data format)
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer='POINT_DATA '//str4//lf
 write(66) trim(buffer)
 write(str1(1:8),'(i8)') nfields
 buffer='FIELD FieldData '//str1//lf
 write(66) trim(buffer)


! write u field
if (uflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'u 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(u(i,j,k))
            enddo
        enddo
    enddo
endif

! write v field
if (vflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'v 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(v(i,j,k))
            enddo
        enddo
    enddo
endif

! write w field
if (wflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'w 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(w(i,j,k))
           enddo
        enddo
    enddo
endif

! write phi field
if (phiflag .eq. 1) then
    write(str4(1:16),'(i16)') numx*numy*numz
    buffer = 'phi 1 '//str4//' float'//lf
    write(66) trim(buffer)
    do k=z_start,z_end,dnz
        do j=y_start,y_end,dny
            do i=x_start,x_end,dnx
                write(66) real(phi(i,j,k))
           enddo
        enddo
    enddo
endif

return
end




! output for fibers 
subroutine generate_fiber_output(nstep,nfib)

use commondata

implicit none 
integer,intent(in) :: nstep,nfib

integer :: m, f, j
integer :: beads_per_fiber
integer :: first_bead, last_bead
integer :: total_connectivity
character(len=100) :: fname
character(len=8) :: numfile


! check fibers topology first

if (nfib <= 0 ) then
    write(*,*) 'Invalid number of fibers'
    stop 
endif
if (mod(NBEADS,nfib) /= 0) then 
    write(*,*) 'Invalid number of beads per fiber'
    write(*,*) 'NBEADS =',NBEADS
    write(*,*) 'NFIBERS =',nfib
    stop
endif

beads_per_fiber = NBEADS/nfib

total_connectivity = nfib *(beads_per_fiber + 1) 


! Construct output filename
write(numfile,'(i8.8)') nstep
fname = './output/FIBERS_'//numfile//'.vtk'
open(unit=91, file=trim(fname), status='replace', action='write', form='formatted')


! VTK header
write(91,'(a)') '# vtk DataFile Version 3.0'
write(91,'(a)') 'Flexible fibers'
write(91,'(a)') 'ASCII'
write(91,'(a)') 'DATASET POLYDATA'

! ---------------------------------------------------------
! Write bead positions
!
! fibers(2,:) = x
! fibers(3,:) = y
! fibers(4,:) = z
! ---------------------------------------------------------

write(91,'(a,1x,i0,1x,a)') 'POINTS', NBEADS, 'double'

do m = 1, NBEADS

    write(91,'(3(es24.16,1x))') &
        fibers(2,m), &
        fibers(3,m), &
        fibers(4,m)

enddo

! ---------------------------------------------------------
! Write fiber connectivity
!
! VTK point indices start from ZERO.
! Fortran array indices start from ONE.
! ---------------------------------------------------------

write(91,'(a,1x,i0,1x,i0)') &
    'LINES', nfib, total_connectivity

do f = 1, nfib

    first_bead = (f-1)*beads_per_fiber + 1
    last_bead  = f*beads_per_fiber
    ! Number of points in this fiber
    write(91,'(i0)',advance='no') beads_per_fiber
    ! Point indices: convert from Fortran to VTK indexing
    do j = first_bead, last_bead
        write(91,'(1x,i0)',advance='no') j-1
    enddo

    write(91,*)

enddo
close(91)

write(*,'(a,a)') ' Wrote fiber VTK: ', trim(fname)

end subroutine
