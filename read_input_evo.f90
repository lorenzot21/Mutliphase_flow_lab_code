! read_input.f90
! Reads simulation parameters from  input.dat  using Fortran NAMELIST.
! Call once at the start of main.f90, before any allocation.
!
! NAMELIST groups in input.dat (must appear in this order):
!   &run_control   – restart flag, topology, sub-iterations, fiber IC flags
!   &flow_control  – run length, dump frequency, flow IC flag
!   &flow_physics  – time step, viscosity, ABC forcing, phase-field
!   &fibers_input  – fiber geometry, mechanics, output frequency
!
! Only  param%nx  (grid resolution) requires recompile when changed.
! All other parameters — including nfibers and lambda — are runtime.
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! Evo_version author: L.Terlicher 
! Date: 10/2026
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz


subroutine read_input
use param,      only: lx, dx, dt, pi, nx
use flowparam
use fibers,     only: nfibers, lambda, nbeads, subiterations, &
                      ic_pos, ic_vel, &
                      a_par, arl, dens_part, E_bend, curv_eq, &
                      fiber_dump, fiber_sdump, nstart, &
                      fiber_types_flag, n_fiber_types, &
                      grouping_rule, max_fiber_types, &
                      a_par_type, arl_type, dens_part_type, &
                      E_bend_type, curv_eq_type, Ffric, &
                      bead_types_flag, bead_asym, zeta_lb, &
                      arl_lb, bead_eq, a0,a1,a2
implicit none
! to make iteration loop on fibers_inp
integer :: f, ios


! ── Namelist declarations ─────────────────────────────────────────────────────

! Run / restart control + fiber topology + IC flags
namelist /run_control/ restart, restart_step, &
                       nfibers, lambda, subiterations, &
                       ic_pos, ic_vel

! Flow run control
namelist /flow_control/ tfin, dump, inflow

! Flow physics
namelist /flow_physics/ dt_flow, mu, rho, f1, f2, f3, k0, radius, sigma

! Fiber physical parameters + output frequency
namelist /fibers_input/ fiber_types_flag, bead_types_flag, &
                        a_par, arl, dens_part, E_bend, & 
                        curv_eq, Ffric, fiber_dump, fiber_sdump, nstart

! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! Namelist of separate fiber file - for multifiber and asymmetry
namelist /fiber_types_input/ n_fiber_types, grouping_rule, &
                             a_par_type, arl_type, dens_part_type, &
                             E_bend_type, curv_eq_type, bead_asym, &
                             zeta_lb, arl_lb, bead_eq, a0,a1,a2
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
                             
! ── Read file (groups must appear in this order in input.dat) ─────────────────
open(unit=10, file='input.dat', status='old', action='read', form='formatted')
read(10, nml=run_control)
read(10, nml=flow_control)
read(10, nml=flow_physics)
read(10, nml=fibers_input)
close(10)

! Validate fiber topology before calculating nbeads
if (nfibers < 1) then
   write(*,'(a,i8)') 'ERROR read_input: nfibers must be positive: ', nfibers
   call exit(1)
endif
if (lambda < 1) then
   write(*,'(a,i8)') 'ERROR read_input: lambda must be positive: ', lambda
   call exit(1)
endif
if (subiterations < 1) then
   write(*,'(a,i8)') 'ERROR read_input: subiterations must be positive: ', &
                     subiterations
   call exit(1)
endif

! ── Derived quantities ────────────────────────────────────────────────────────
nbeads = nfibers * lambda          ! total beads across all fibers

pi  = 4.d0 * datan(1.d0)
lx  = 2.d0 * pi                   ! periodic box side = [0, 2π]
dx  = lx / (dble(nx) - 1.d0)
dt  = dt_flow

! added for reading other file zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!if (fiber_types_flag == 1) then  
!open(unit=99, file='fibers_inp.dat', &
!        status='old', action='read', form='formatted')
!   read(99, nml=fiber_types_input)
!   close(99)
!endif

! check on fiber types flag
if (fiber_types_flag /= 0 .and. fiber_types_flag /= 1) then
   write(*,'(a,i6)') &
      'ERROR read_input: invalid fiber_types_flag = ', fiber_types_flag
   call exit(1)
endif

if (bead_types_flag /= 0 .and. bead_types_flag /= 1) then
   write(*,'(a,i6)') &
      'ERROR read_input: invalid bead_types_flag = ', bead_types_flag
   call exit(1)
endif

! to avoid cross - initialization
if (fiber_types_flag == 1 .and. bead_types_flag == 1) then
   write(*,'(a)') &
      'ERROR: fiber_types_flag and bead_types_flag cannot both be 1'
   error stop 1
endif

!added to read with also a status variable that tells if file read or not
if (fiber_types_flag == 1 .or. bead_types_flag == 1) then

   open(unit=99, file='fibers_inp.dat', &
        status='old', action='read', form='formatted', iostat=ios)

   if (ios /= 0) then
      write(*,'(a)') &
         'ERROR read_input: cannot open fibers_inp.dat'
      call exit(1)
   endif

   read(99, nml=fiber_types_input, iostat=ios)

   if (ios /= 0) then
      write(*,'(a)') &
         'ERROR read_input: failed to read fiber_types_input from fibers_inp.dat'
      close(99)
      call exit(1)
   endif

   close(99)

endif
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz


! ── Consistency checks ────────────────────────────────────────────────────────
if (mod(nbeads, lambda) /= 0) then
   write(*,'(a,2i6)') 'ERROR read_input: nbeads not divisible by lambda:', &
        nbeads, lambda
   call exit(1)
endif

if (fiber_types_flag == 1) then 
   ! Validate the number of fiber types
   if (n_fiber_types < 1 .or. n_fiber_types > max_fiber_types) then
      write(*,'(a,i6)') &
         'ERROR read_input: invalid n_fiber_types = ', n_fiber_types
      call exit(1)
   endif

   ! Ensure fibers can be distributed equally among types
   if (mod(nfibers, n_fiber_types) /= 0) then
      write(*,'(a,2i8)') &
         'ERROR read_input: nfibers must be divisible by n_fiber_types:', &
         nfibers, n_fiber_types
      call exit(1)
   endif
   
elseif(bead_types_flag == 1) then 

   if (bead_asym /= 0 .and. bead_asym /= 1) then
   write(*,'(a,i6)') &
      'ERROR read_input: invalid bead_asym = ', bead_asym
   call exit(1)
   endif

endif

if (ic_pos == 1 .and. restart_step <= 0) then
   write(*,'(a)') 'ERROR read_input: ic_pos=1 requires restart_step > 0'
   call exit(1)
endif

if (restart == 1 .and. restart_step <= 0) then
   write(*,'(a)') 'ERROR read_input: restart=1 requires restart_step > 0'
   call exit(1)
endif


! ── Echo to stdout ────────────────────────────────────────────────────────────
 !enable/disable for debug check parameters
write(*,*) "----------------------------------------------"
write(*,*) "███    ███ ██   ██ ██ ████████ ██████   ██████"  
write(*,*) "████  ████ ██   ██ ██    ██         ██ ██     "       
write(*,*) "██ ████ ██ ███████ ██    ██     █████  ███████"  
write(*,*) "██  ██  ██ ██   ██ ██    ██         ██ ██    ██" 
write(*,*) "██      ██ ██   ██ ██    ██    ██████   ██████"         
write(*,*) "----------------------------------------------"
write(*,*) "----Homogeneous isotropic turbulence setup----"
write(*,'(a)') ''
write(*,'(a)') ' ════════════════ Input parameters ════════════════'
write(*,'(a)')       '  ── run / restart ──'
write(*,'(a,i4)')    '  restart       = ', restart
write(*,'(a,i8)')    '  restart_step  = ', restart_step
write(*,'(a,i8)')    '  nfibers       = ', nfibers
write(*,'(a,i6)')    '  lambda        = ', lambda
write(*,'(a,i8)')    '  nbeads        = ', nbeads
write(*,'(a,i4)')    '  subiterations = ', subiterations
write(*,'(a,i4)')    '  ic_pos        = ', ic_pos, &
                     '  (0=random 1=file)'
write(*,'(a,i4)')    '  ic_vel        = ', ic_vel, &
                     '  (0=zero 1=fluid, ignored when ic_pos=1)'
write(*,'(a)')       '  ── flow control ──'
write(*,'(a,i8)')    '  tfin          = ', tfin
write(*,'(a,i8)')    '  dump          = ', dump
write(*,'(a,i4)')    '  inflow        = ', inflow
write(*,'(a)')       '  ── flow physics ──'
write(*,'(a,es12.4)') '  dt            = ', dt
write(*,'(a,es12.4)') '  mu            = ', mu
write(*,'(a,f8.4)')   '  rho           = ', rho
write(*,'(a,3f8.4)')  '  f1 f2 f3      = ', f1, f2, f3
write(*,'(a,f8.4)')   '  k0            = ', k0
write(*,'(a)')        '  ── fibers ──'
write(*,'(a,i4)')    '   fiber_types_flag       = ', fiber_types_flag
write(*,'(a,i4)')    '   bead_types_flag       = ',  bead_types_flag
write(*,'(a,es12.4)') '  a_par         = ', a_par
write(*,'(a,f8.4)')   '  arl           = ', arl
write(*,'(a,f8.4)')   '  dens_part     = ', dens_part
write(*,'(a,es12.4)') '  E_bend        = ', E_bend
write(*,'(a,f8.4)')   '  curv_eq       = ', curv_eq
write(*,'(a,f8.4)')   '  Ffric         = ', Ffric
write(*,'(a,i8)')     '  fiber_dump    = ', fiber_dump
write(*,'(a,i8)')     '  fiber_sdump   = ', fiber_sdump
write(*,'(a)') ' ══════════════════════════════════════════════════'
write(*,'(a)') ''

! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
if (fiber_types_flag == 1) then

   write(*,'(a,i4)') '  n_fiber_types = ', n_fiber_types
   write(*,'(a,i4)') '  grouping_rule = ', grouping_rule

   do f = 1, n_fiber_types

      write(*,*)
      write(*,'(a,i4)') '  Fiber type = ', f
      write(*,'(a,es12.4)') '  a_par      = ', a_par_type(f)
      write(*,'(a,f8.4)')   '  arl        = ', arl_type(f)
      write(*,'(a,f8.4)')   '  dens_part  = ', dens_part_type(f)
      write(*,'(a,es12.4)') '  E_bend     = ', E_bend_type(f)
      write(*,'(a,f8.4)')   '  curv_eq    = ', curv_eq_type(f)

   enddo
   write(*,'(a)') ' ══════════════════════════════════════════════════'
endif


if (bead_types_flag == 1) then

   write(*,'(a,i4)') '  bead_asym = ', bead_asym
   write(*,'(a,es12.4)') '  zeta_lb      = ', zeta_lb
   write(*,'(a,f8.4)')   '  arl_lb        = ', arl_lb
   write(*,'(a,i4)') '   bead_eq = ', bead_eq
   write(*,'(a,3i4)') ' a0, a1, a2 =', a0,a1,a2
   write(*,'(a)') ' ══════════════════════════════════════════════════'
endif

end subroutine read_input
