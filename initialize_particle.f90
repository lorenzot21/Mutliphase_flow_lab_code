! initialize_particle.f90
! Fiber initialization for the triperiodic HIT solver.
!
! Controlled by four flags read from input.dat  &run_control:
!
!   restart  restart_step  ic_pos  ic_vel   behaviour
!   ────────────────────────────────────────────────────────────────────────
!      1        N            *       *       Full restart: read fiber_N.dat
!                                            (positions AND velocities from file)
!      0        *            1       *       Fiber-only restart: read fiber_N.dat
!                                            (positions AND velocities from file;
!                                             flow field may still be fresh)
!      0        *            0       0       Fresh start: random positions,
!                                            zero velocities
!      0        *            0       1       Fresh start: random positions,
!                                            fluid velocity at bead positions
!
! Call order: read_input → initialize_fibers (→ get_fiber_fluid_quantities
!             only if ic_vel=1 and ic_pos=0, i.e. flow field must be loaded first)


! Multifiber_evo mmmmmmmmmmmmmmmmm
! Author : L.Terlicher
! Date : 10/2026
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm


subroutine initialize_fibers
use param
!use flowparam,  only: restart, restart_step
use flowparam,  only: mu, rho, k0, restart, restart_step, radius
use velocity,   only: u, v, w
use fibers
implicit none

integer :: i, j, k, f, b
double precision :: phi_ang, theta_ang, fac, norm, LLL_loc, LLL_loc_lb
double precision :: LLL_loc_type(max_fiber_types)
!added
double precision :: nu_f, tau_p, L_int_f


! ── Allocate global fiber arrays ─────────────────────────────────────────────
allocate(xf(nbeads, 3))
allocate(angles_f(nbeads, 3))
allocate(uf(nbeads, 3))
allocate(omgf(nbeads, 3))
allocate(ff(nbeads, 3))
allocate(vp_o(nbeads, 3, 2))
allocate(angl_o(nbeads, 3, 2))
allocate(uf_at_p(nbeads, 3))
allocate(vort_at_p(nbeads, 3))
allocate(grad_at_p(nbeads, 9))
allocate(kzero(3, nbeads))
allocate(qzero(3, 3, nbeads))
allocate(rzero(3, nbeads))
allocate(szero(3, 3, nbeads))
! added
allocate(X_constraint(nbeads, 3))
allocate(torque_bend(nbeads, 3))
!allocate(coll_partner(nbeads))phi_ang
allocate(coll_partner(nbeads))
allocate(bead_on_interface(nbeads))
allocate(ndetach(nbeads))
allocate(n_iface(nbeads, 3))
allocate(phi_at_p(nbeads))
allocate(dphi_at_p(nbeads, 3))
allocate(f_reaction(nbeads, 3))

! Zero-initialise
xf=0.d0; angles_f=0.d0; uf=0.d0; omgf=0.d0; ff=0.d0
vp_o=0.d0; angl_o=0.d0
uf_at_p=0.d0; vort_at_p=0.d0; grad_at_p=0.d0
kzero=0.d0; qzero=0.d0; rzero=0.d0; szero=0.d0
! added
X_constraint=0.d0; torque_bend=0.d0; coll_partner=0
bead_on_interface=0; n_iface=0.d0; phi_at_p=0.d0; dphi_at_p=0.d0
f_reaction=0.d0; ndetach=0

! ── Derived time step ────────────────────────────────────────────────────────
dt_part = dt / dble(subiterations)

! ── Geometric helpers ────────────────────────────────────────────────────────
Idm = 0.d0
do i = 1, 3
   Idm(i,i) = 1.d0
enddo
epsi = 0.d0
epsi(1,2,3)= 1.d0;  epsi(2,3,1)= 1.d0;  epsi(3,1,2)= 1.d0
epsi(3,2,1)=-1.d0;  epsi(1,3,2)=-1.d0;  epsi(2,1,3)=-1.d0


! calculation of the fiber distribution mmmmmmmmmmmmmmmmmmm

!if (bead_types_flag == 1) then

!   bead_length = fiber_length / 4
!   if (bead_eq == 1) then 
!      do b = 1,lambda
         ! linear equation to define the bead semiradius
!         a_par_bead(i) = a0 + a1*(i/2)     
!         arl_bead(i) = a_par_bead(i) / bead_length
      ! to check the stiffness and other terms that there are
!      enddo 
!   elseif ( bead_eq == 2)
!      do b = 1,lambda 
!         a_par_bead(i) = a0 + a1*(i/2) + a2*((i/2)**2)
!         arl_bead(i) = a_par_bead(i) / bead_length
      !!! to check the rest
!      enddo 
!   else 
   ! return an errror because no other formulations implemented

!endif


! ── Resistance coefficients ─────────────────────────────────────────────────
b_ff  = (a_par / 1.24d0) * sqrt(log(arl))
ec_ff = sqrt(((a_par*arl)**2.d0 - b_ff**2.d0) / (a_par*arl)**2.d0)
LLL_loc = log((1.d0 + ec_ff) / (1.d0 - ec_ff))

X_a = 8.d0/3.d0 * ec_ff**3.d0 / &
      (-2.d0*ec_ff + (1.d0 + ec_ff**2.d0)*LLL_loc)
Y_a = 16.d0/3.d0 * ec_ff**3.d0 / &
      (2.d0*ec_ff + (3.d0*ec_ff**2.d0 - 1.d0)*LLL_loc)
X_c = 4.d0/3.d0 * ec_ff**3.d0 * (1.d0 - ec_ff**2.d0) / &
      (2.d0*ec_ff - (1.d0 - ec_ff**2.d0)*LLL_loc)
Y_c = 4.d0/3.d0 * ec_ff**3.d0 * (2.d0 - ec_ff**2.d0) / &
      (-2.d0*ec_ff + (ec_ff**2.d0 + 1.d0)*LLL_loc)
Y_h = 4.d0/3.d0 * ec_ff**5.d0 / &
      (-2.d0*ec_ff + (1.d0 + ec_ff**2.d0)*LLL_loc)

! to add also the first bead coefs in case of asym
if (bead_types_flag == 1)then
   if (bead_asym == 0) then
   
      ! calculation of last bead semiradius absolute value
      a_par_lb = a_par*zeta_lb 
      
      ! calculation of resistance coefs for last bead
      b_ff_lb  = (a_par_lb / 1.24d0) * sqrt(log(arl_lb))
      ec_ff_lb = sqrt(((a_par_lb*arl_lb)**2.d0 - b_ff_lb**2.d0) / (a_par_lb*arl_lb)**2.d0)
      LLL_loc_lb = log((1.d0 + ec_ff_lb) / (1.d0 - ec_ff_lb))

      X_a_lb = 8.d0/3.d0 * ec_ff_lb**3.d0 / &
            (-2.d0*ec_ff_lb + (1.d0 + ec_ff_lb**2.d0)*LLL_loc_lb)
      Y_a_lb = 16.d0/3.d0 * ec_ff_lb**3.d0 / &
            (2.d0*ec_ff_lb + (3.d0*ec_ff_lb**2.d0 - 1.d0)*LLL_loc_lb)
      X_c_lb = 4.d0/3.d0 * ec_ff_lb**3.d0 * (1.d0 - ec_ff_lb**2.d0) / &
            (2.d0*ec_ff_lb - (1.d0 - ec_ff_lb**2.d0)*LLL_loc_lb)
      Y_c_lb = 4.d0/3.d0 * ec_ff_lb**3.d0 * (2.d0 - ec_ff_lb**2.d0) / &
            (-2.d0*ec_ff_lb + (ec_ff_lb**2.d0 + 1.d0)*LLL_loc_lb)
      Y_h_lb = 4.d0/3.d0 * ec_ff_lb**5.d0 / &
            (-2.d0*ec_ff_lb + (1.d0 + ec_ff_lb**2.d0)*LLL_loc_lb)
   
   elseif (bead_asym /= 0) then 
      write(*,*) ' ERROR : put bead_asym to 0'
      call exit(1)
   endif
endif
 
if (gravity_flag == 1) then
   geff(1) = (1.d0 - 1.d0/dens_part) * grav_x
   geff(2) = (1.d0 - 1.d0/dens_part) * grav_y
   geff(3) = (1.d0 - 1.d0/dens_part) * grav_z
else
   geff = 0.d0
endif

write(*,'(a)') ' --- Flexible fiber initialisation ---'
write(*,'(a,i6,a,i4,a,i8)') &
     '     nfibers=', nfibers, '  lambda=', lambda, '  nbeads=', nbeads
write(*,'(a,f10.4,a,f10.4)') '     a_par=', a_par, '   arl=', arl
write(*,'(a,5f10.4)') '     Xa Ya Xc Yc Yh =', X_a, Y_a, X_c, Y_c, Y_h
if (bead_types_flag == 1)then
if (bead_asym == 0) then
write(*,'(a,f10.4,a,f10.4)') '     a_par_lb =', a_par_lb, '   arl_lb =', arl_lb
write(*,'(a,5f10.4)') 'Xa_lb Ya_lb Xc_lb Yc_lb Yh_lb =', X_a_lb, Y_a_lb, X_c_lb, Y_c_lb, Y_h_lb
endif
endif

! Multiple fiber types routine zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
if (fiber_types_flag == 1) then

   if (mod(nfibers, n_fiber_types) /= 0) then
      write(*,'(a)') 'ERROR: nfibers must be divisible by n_fiber_types.'
      stop 1
   endif
   
   allocate(fiber_type(nfibers))
   type_idx = nfibers / n_fiber_types
   !defining the fiber characteristics for all fiber types
   do i = 1, nfibers
      fiber_type(i) = (i-1)/type_idx + 1
   enddo

   do f=1,n_fiber_types 

      b_ff_type(f)  = (a_par_type(f) / 1.24d0) * sqrt(log(arl_type(f)))
      ec_ff_type(f) = sqrt(((a_par_type(f)*arl_type(f))**2.d0 - b_ff_type(f)**2.d0) / (a_par_type(f)*arl_type(f))**2.d0)
      LLL_loc_type(f) = log((1.d0 + ec_ff_type(f)) / (1.d0 - ec_ff_type(f)))

      ! ── Translational resistance coefficients ─────────────────────
      X_a_type(f) = (8.d0 / 3.d0) * ec_ff_type(f)**3.d0 / &
       ( -2.d0 * ec_ff_type(f) + &
         (1.d0 + ec_ff_type(f)**2.d0) * LLL_loc_type(f) )
      Y_a_type(f) = (16.d0 / 3.d0) * ec_ff_type(f)**3.d0 / &
       ( 2.d0 * ec_ff_type(f) + &
         (3.d0 * ec_ff_type(f)**2.d0 - 1.d0) * LLL_loc_type(f) )
      ! ── Rotational resistance coefficients ────────────────────────
      X_c_type(f) = (4.d0 / 3.d0) * ec_ff_type(f)**3.d0 * &
       (1.d0 - ec_ff_type(f)**2.d0) / &
       ( 2.d0 * ec_ff_type(f) - &
         (1.d0 - ec_ff_type(f)**2.d0) * LLL_loc_type(f) )
      Y_c_type(f) = (4.d0 / 3.d0) * ec_ff_type(f)**3.d0 * &
       (2.d0 - ec_ff_type(f)**2.d0) / &
       ( -2.d0 * ec_ff_type(f) + &
         (ec_ff_type(f)**2.d0 + 1.d0) * LLL_loc_type(f) )
      Y_h_type(f) = (4.d0 / 3.d0) * ec_ff_type(f)**5.d0 / &
       ( -2.d0 * ec_ff_type(f) + &
         (1.d0 + ec_ff_type(f)**2.d0) * LLL_loc_type(f) )

      if (gravity_flag == 1) then
         geff_type(f,1) = (1.d0 - 1.d0/dens_part_type(f)) * grav_x
         geff_type(f,2) = (1.d0 - 1.d0/dens_part_type(f)) * grav_y
         geff_type(f,3) = (1.d0 - 1.d0/dens_part_type(f)) * grav_z
         else
         geff_type(f,:) = 0.d0
      endif

      write(*,'(a)') ' --- Flexible fiber initialisation ---'
      write(*,'(a,i6,a,i4,a,i8)') &
         '     nfibers=', nfibers, '  lambda=', lambda, '  nbeads=', nbeads
      write(*,'(a,f10.4,a,f10.4)') '     a_par=', a_par_type(f), '   arl=', arl_type(f)
      write(*,'(a,5f10.4)') '     Xa Ya Xc Yc Yh =', X_a_type(f), Y_a_type(f), X_c_type(f), Y_c_type(f), Y_h_type(f)                       

   enddo
endif

! added zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
mu_f = mu









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ── Rigid slender-body coefficients (used only when fiber_model=1) ───────────
! The lambda beads double as Nq quadrature points of a rigid rod.  A bead n sits
! at  Xbar + s_quad(n)*p , so the signed arc-lengths must reproduce the flexible
! chain geometry (bead spacing 2*arl*a_par, so L = lambda*2*arl*a_par).  With
!   s_quad(n) = (2n-lambda-1)*arl*a_par
! the spacing is 2*arl*a_par and sum_n s_quad(n) = 0 exactly (CM at the middle).
if ( fiber_model == 1) then
   allocate(s_quad(lambda))
   do i = 1, lambda
      s_quad(i) = dble(2*i - lambda - 1) * arl * a_par
   enddo
   S2_quad = sum(s_quad**2)                       ! least-squares slope normalisation
   Lfib    = dble(lambda) * 2.d0 * arl * a_par     ! rod contour length

   ! Mass = sum of the lambda bead masses (keeps momentum bookkeeping identical to
   ! the flexible/two-way path, where bead_mass = (4/3)*pi*arl*a_par^3*dens_part
   ! with rho_f=1).  I_T is the transverse moment of inertia of a uniform rod.
   Mfib    = dble(lambda) * (4.d0/3.d0)*pi*arl*a_par**3 * dens_part
   I_T_fib = Mfib * Lfib**2 / 12.d0

   ! Slender-body resistance (Batchelor 1970 / Cox 1970).  lam_tot = L/(2a) = total
   ! rod aspect ratio = lambda*arl.
   LLL_loc   = log(2.d0 * dble(lambda) * arl)      ! log(2*lam_tot); reuse temp
   zeta_par  = 2.d0*pi*mu_f*Lfib / (LLL_loc - 0.5d0)
   zeta_perp = 4.d0*pi*mu_f*Lfib / (LLL_loc + 0.5d0)
   zeta_rot  = pi*mu_f*Lfib**3 / (3.d0*(LLL_loc - 0.5d0))
   tau_par   = Mfib / zeta_par
   tau_perp  = Mfib / zeta_perp
   tau_rot   = I_T_fib / zeta_rot

   ! Effective (buoyancy-reduced) gravity.  geff=0 unless gravity_flag=1, so the
   ! validated zero-gravity behaviour is preserved byte-for-byte when it is off.
   if (gravity_flag == 1) then
      geff(1) = (1.d0 - 1.d0/dens_part) * grav_x
      geff(2) = (1.d0 - 1.d0/dens_part) * grav_y
      geff(3) = (1.d0 - 1.d0/dens_part) * grav_z
   else
      geff = 0.d0
   endif
endif

if (fiber_model == 1) then
   write(*,'(a)') ' --- Rigid slender-body fiber coefficients ---'
   write(*,'(a,es12.4,a,es12.4)') '     L =', Lfib, '   M =', Mfib
   write(*,'(a,es12.4)')          '     I_T =', I_T_fib
   write(*,'(a,3es12.4)')         '     zeta_par/perp/rot =', zeta_par, zeta_perp, zeta_rot
   write(*,'(a,3es12.4)')         '     tau_par/perp/rot  =', tau_par, tau_perp, tau_rot
   write(*,'(a,es12.4)')          '     S2_quad =', S2_quad
   if (gravity_flag == 1) &
      write(*,'(a,3es12.4)')      '     geff =', geff(1), geff(2), geff(3)
endif

write(*,'(a)') ' --- Flexible fiber initialisation ---'
write(*,'(a,i6,a,i4,a,i8)') &
     '     nfibers=', nfibers, '  lambda=', lambda, '  nbeads=', nbeads
write(*,'(a,f10.4,a,f10.4)') '     a_par=', a_par, '   arl=', arl
write(*,'(a,5f10.4)') '     Xa Ya Xc Yc Yh =', X_a, Y_a, X_c, Y_c, Y_h

! ── Dimensionless numbers at startup ─────────────────────────────────────────
! tau_p = 2*a^2*dens_ratio / (9*nu)    [bead relaxation time, dens_ratio = rho_p/rho_f]
! tau_f = sqrt(nu/eps)                 [Kolmogorov time; eps from flow stats at runtime]
! Re_L  = u_rms * L_int / nu           [L_int = 2*pi/k0: ABC forcing integral scale]
! Re_lam= u_rms * lambda_T / nu        [lambda_T = sqrt(15*nu/eps)*u_rms]
! St_eta= tau_p / tau_f                [particle Stokes number at Kolmogorov scale]
nu_f    = mu / rho
tau_p   = 2.d0*a_par**2.d0*dens_part / (9.d0*nu_f)
L_int_f = 2.d0*pi / k0
write(*,'(a)') ' --- Dimensionless numbers (outer units) ---'
write(*,'(a,es12.4)') '     nu = mu/rho                              = ', nu_f
write(*,'(a)')        '     tau_p = 2*a^2*dens_ratio/(9*nu)          [dens_ratio = rho_p/rho_f]'
write(*,'(a,es12.4)') '                                              = ', tau_p
write(*,'(a)')        '     tau_f = sqrt(nu/eps)                     [eps from flow stats]'
write(*,'(a,es12.4)') '     L_int = 2*pi/k0  (ABC forcing scale)     = ', L_int_f
write(*,'(a)') '     Re_L      = u_rms * L_int / nu           [u_rms from flow stats]'
write(*,'(a)') '     Re_lambda = u_rms * sqrt(15*nu/eps) / nu  [eps from flow stats]'
write(*,'(a)') '     St_eta    = tau_p / tau_f                 [eps from flow stats]'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









! ═══════════════════════════════════════════════════════════════════════════════
!  Branch on ic_pos / restart flags
! ═══════════════════════════════════════════════════════════════════════════════

if (restart == 1 .or. ic_pos == 1) then

   ! ── Restart from file ───────────────────────────────────────────────────────
   ! Both positions AND velocities are read from the restart snapshot.
   ! ic_vel is ignored.
   write(*,'(a,i8,a)') '     Reading fiber restart from step ', restart_step, &
        '  (./output/fiber_XXXXXXXX.dat)'
   call read_fibers(restart_step)

   ! here there is something to do do characterize fibers (done in flexible fibers routine)

! added zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

else if (nfibers == 1) then

   ! ── Validation IC: single fiber tangent to the drop interface ───────────────
   ! For the two-way-coupling bulge test.  The drop is centred at the box centre
   ! with its interface at radius `radius`.  Lay the lone fiber along z, ~1 cell
   ! OUTSIDE the +x interface point, so its beads sit in the attraction zone
   ! (phi just below 0.5) and pull the interface toward them.  Deterministic so
   ! the test is reproducible (a single random fiber would almost never land on
   ! the interface).
   write(*,'(a)') '     nfibers=1: deterministic placement tangent to drop interface (validation)'
   angles_f(:,1) = 0.d0;  angles_f(:,2) = 0.d0;  angles_f(:,3) = 1.d0   ! along z
   xf(1,1) = lx/2.d0 + radius + dx                       ! ~1 cell outside +x interface
   xf(1,2) = lx/2.d0
   xf(1,3) = lx/2.d0 - dble(lambda-1)*arl*a_par          ! so the chain centres on lx/2
   do j = 1, lambda-1
      xf(1+j,1) = xf(j,1)
      xf(1+j,2) = xf(j,2)
      xf(1+j,3) = xf(j,3) + 2.d0*arl*a_par
   enddo
   uf = 0.d0;  omgf = 0.d0
   !!!!!!!!!!!!!!

else

   ! ── Fresh start: random positions ──────────────────────────────────────────
   write(*,'(a)') '     ic_pos=0: random fiber placement'

   ! multifiber mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
   if (fiber_types_flag == 1) then

      ! Equal-sized consecutive groups of fiber types.
      if (mod(nfibers, n_fiber_types) /= 0) then
         write(*,'(a)') &
            'ERROR: nfibers must be divisible by n_fiber_types.'
         call exit(1)
      endif
      ! actually non stricly necessary

      type_idx = nfibers / n_fiber_types

      ! Determine the type of the current fiber.
      do i = 1, nfibers

         f = fiber_type(i)
         ! so with this we get that till type_idx the fibers will belong to 1st type and so on with the others

         ! Random orientation — uniform on the sphere - NB here changed from previous case
         !in fact calls are alternated with angles definition not to overwrite the angle value
         call random_number(fac)
         phi_ang = 2.d0*pi*fac
         call random_number(fac)
         theta_ang = acos(2.d0*fac - 1.d0)
         angles_f(1+(i-1)*lambda,1) = sin(theta_ang)*cos(phi_ang)
         angles_f(1+(i-1)*lambda,2) = sin(theta_ang)*sin(phi_ang)
         angles_f(1+(i-1)*lambda,3) = cos(theta_ang)
         norm = sqrt(sum(angles_f(1+(i-1)*lambda,:)**2.d0))
         if (norm > 1.d-14) then
            angles_f(1+(i-1)*lambda,:) = angles_f(1+(i-1)*lambda,:)/norm
         endif
         ! Assign the same orientation to every bead in this fiber.
         do j = 1, lambda-1
            angles_f(1+(i-1)*lambda+j,:) = angles_f(1+(i-1)*lambda,:)
         enddo

         ! Random position of the first bead.
         call random_number(xf(1+(i-1)*lambda,:))
         xf(1+(i-1)*lambda,:) = xf(1+(i-1)*lambda,:) * lx

         ! Build the chain using the geometry of its assigned type.
         do j = 1, lambda-1
         do k = 1, 3

            xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j-1,k) + &
               2.d0*arl_type(f)*a_par_type(f)*angles_f(1+(i-1)*lambda+j-1,k)

         enddo
            ! Periodic wrapping.
         do k = 1, 3
               if (xf(1+(i-1)*lambda+j,k) < 0.d0) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) + lx
               if (xf(1+(i-1)*lambda+j,k) > lx) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) - lx
         enddo
         enddo
      enddo
   
   !standard case executed only when the flag is not 1 
   elseif (bead_types_flag == 1) then 

      do i = 1, nfibers

         ! Random orientation — uniform on sphere
         call random_number(fac);  phi_ang   = 2.d0*pi*fac
         call random_number(fac);  theta_ang = acos(2.d0*fac - 1.d0)
         angles_f(1+(i-1)*lambda, 1) = sin(theta_ang)*cos(phi_ang)
         angles_f(1+(i-1)*lambda, 2) = sin(theta_ang)*sin(phi_ang)
         angles_f(1+(i-1)*lambda, 3) = cos(theta_ang)
         norm = sqrt(sum(angles_f(1+(i-1)*lambda,:)**2.d0))
         if (norm > 1.d-14) angles_f(1+(i-1)*lambda,:) = angles_f(1+(i-1)*lambda,:)/norm
         do j = 1, lambda-1
            angles_f(1+(i-1)*lambda+j, :) = angles_f(1+(i-1)*lambda, :)
         enddo

         ! Random position for the first bead
         call random_number(xf(1+(i-1)*lambda, :))
         xf(1+(i-1)*lambda, :) = xf(1+(i-1)*lambda, :) * lx

         ! Build chain — bead-to-bead step = 2 * arl * a_par along orientation
         do j = 1, lambda-1 
            do k = 1, 3
               ! to give the last bead a different dimension
               if (j == lambda - 1 ) then
                  xf(1+(i-1)*lambda+j, k) = xf(1+(i-1)*lambda+j-1, k) &
                     + 2.d0*arl_lb*a_par_lb * angles_f(1+(i-1)*lambda+j-1, k)   
               else
                  xf(1+(i-1)*lambda+j, k) = xf(1+(i-1)*lambda+j-1, k) &
                     + 2.d0*arl*a_par * angles_f(1+(i-1)*lambda+j-1, k)
               endif 
            enddo
            do k = 1, 3
               if (xf(1+(i-1)*lambda+j,k) < 0.d0) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) + lx
               if (xf(1+(i-1)*lambda+j,k) > lx) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) - lx
            enddo
         enddo
      enddo

   else
   
      do i = 1, nfibers

         ! Random orientation — uniform on sphere
         call random_number(fac);  phi_ang   = 2.d0*pi*fac
         call random_number(fac);  theta_ang = acos(2.d0*fac - 1.d0)
         angles_f(1+(i-1)*lambda, 1) = sin(theta_ang)*cos(phi_ang)
         angles_f(1+(i-1)*lambda, 2) = sin(theta_ang)*sin(phi_ang)
         angles_f(1+(i-1)*lambda, 3) = cos(theta_ang)
         norm = sqrt(sum(angles_f(1+(i-1)*lambda,:)**2.d0))
         if (norm > 1.d-14) angles_f(1+(i-1)*lambda,:) = angles_f(1+(i-1)*lambda,:)/norm
         do j = 1, lambda-1
            angles_f(1+(i-1)*lambda+j, :) = angles_f(1+(i-1)*lambda, :)
         enddo

         ! Random position for the first bead
         call random_number(xf(1+(i-1)*lambda, :))
         xf(1+(i-1)*lambda, :) = xf(1+(i-1)*lambda, :) * lx

         ! Build chain — bead-to-bead step = 2 * arl * a_par along orientation
         do j = 1, lambda-1 
            do k = 1, 3
               xf(1+(i-1)*lambda+j, k) = xf(1+(i-1)*lambda+j-1, k) &
                  + 2.d0*arl*a_par * angles_f(1+(i-1)*lambda+j-1, k)
            enddo
            do k = 1, 3
               if (xf(1+(i-1)*lambda+j,k) < 0.d0) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) + lx
               if (xf(1+(i-1)*lambda+j,k) > lx) &
                  xf(1+(i-1)*lambda+j,k) = xf(1+(i-1)*lambda+j,k) - lx
            enddo
         enddo
      enddo
   endif

   ! ── Velocity initial condition ──────────────────────────────────────────────
   if (ic_vel == 0) then
      write(*,'(a)') '     ic_vel=0: zero initial bead velocities'
      uf   = 0.d0
      omgf = 0.d0

   else  ! ic_vel == 1
      write(*,'(a)') '     ic_vel=1: interpolating fluid velocity to bead positions'
      call get_fiber_fluid_quantities   ! fills uf_at_p
      uf   = uf_at_p                   ! initialise bead velocity = fluid velocity
      omgf = 0.d0                      ! angular velocity starts at zero

   endif

endif

! ── Seed AB time-history from initial velocity ───────────────────────────────
! This ensures the integrator starts cleanly regardless of IC path.
vp_o(:,:,1) = uf
vp_o(:,:,2) = uf
do i = 1, nbeads
   angl_o(i,1,1) = angles_f(i,2)*omgf(i,3) - omgf(i,2)*angles_f(i,3)
   angl_o(i,2,1) = -(angles_f(i,1)*omgf(i,3) - omgf(i,1)*angles_f(i,3))
   angl_o(i,3,1) = angles_f(i,1)*omgf(i,2) - omgf(i,1)*angles_f(i,2)
   angl_o(i,:,2) = angl_o(i,:,1)
enddo

if (fiber_types_flag == 1) then 

   do f = 1,n_fiber_types
      
      write(*,'(a,i8,f10.4)') ' Fiber length of type ', f, &
           dble(lambda) * 2.d0 * arl_type(f) * a_par_type(f)
   enddo
   write(*,'(a)') ' Initializzation complete'

else

   write(*,'(a,f10.4)') '     Fiber length (outer units) =', &
         dble(lambda)*2.d0*arl*a_par
   write(*,'(a)') '     Initialization complete.'

endif


end subroutine initialize_fibers
