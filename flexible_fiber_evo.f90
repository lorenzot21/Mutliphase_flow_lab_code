! flexible_fiber.f90
! Flexible fiber solver for triperiodic HIT using the Gears bead-chain model.
!
! Physics reference:
!   Delmotte et al. 2015, "A general formulation of Bead Models applied to
!   flexible fibers and active filaments at low Reynolds number"
!   Original channel-flow implementation: D. Di Giusto (Jan 2021), Univ. of Udine
!
! Adaptations for HIT:
!   - Fully periodic domain (no walls, no wall-collision check)
!   - Pure Fortran + OpenMP (no MPI, no CUDA/OpenACC)
!   - Global fiber memory via the `fibers` module (no wind_* scatter/gather)
!   - Outer-unit scaling  (re = 1); dt_part*re = dt_part throughout
!   - One-way coupling (Stokes drag, Rep = 0)
!
! Entry points:
!   get_fiber_fluid_quantities   -- trilinear interp of u, vorticity, grad-u
!   advance_fibers(step)         -- main OpenMP-parallel time-step loop

! Friction routine zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! author : L.Terlicher
! date : 10/2026
! NB friction routine implemented just for the inter-fiber collisions
!    with intra-fiber collision not addressed at the moment

! Multifiber routine mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! author : L.Terlicher
! date : 10/2026


! ============================================================================
!  get_fiber_fluid_quantities
!  Interpolate fluid velocity, vorticity, and velocity-gradient tensor at
!  every bead position.  Uses the same trilinear stencil as particles.f90.
! ============================================================================
subroutine get_fiber_fluid_quantities
use param
use velocity
use fibers
implicit none

integer :: m, ic, jc, kc, icp, jcp, kcp
integer :: cx, cy, cz, ix, iy, iz, ipp, jpp, kpp, imm, jmm, kmm
double precision :: xt, yt, zt, dxi
double precision :: c1,c2,c3,c4,c5,c6,c7,c8
double precision :: xb, yb, zb
! gradient arrays at 8 cell corners (cx+1, cy+1, cz+1)
double precision :: du_dx(2,2,2), du_dy(2,2,2), du_dz(2,2,2)
double precision :: dv_dx(2,2,2), dv_dy(2,2,2), dv_dz(2,2,2)
double precision :: dw_dx(2,2,2), dw_dy(2,2,2), dw_dz(2,2,2)

dxi = 1.0d0 / dx

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP   SHARED(nbeads, xf, uf_at_p, vort_at_p, grad_at_p, u, v, w, dx, lx, nx, dxi) &
!$OMP   PRIVATE(m, ic, jc, kc, icp, jcp, kcp, xt, yt, zt, &
!$OMP           c1,c2,c3,c4,c5,c6,c7,c8, xb, yb, zb, &
!$OMP           cx, cy, cz, ix, iy, iz, ipp, jpp, kpp, imm, jmm, kmm, &
!$OMP           du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)
do m = 1, nbeads
   xb = xf(m,1);  yb = xf(m,2);  zb = xf(m,3)

   ! Cell index (0-based) and fractional position
   ic = floor(xb * dxi);  xt = xb*dxi - dble(ic)
   jc = floor(yb * dxi);  yt = yb*dxi - dble(jc)
   kc = floor(zb * dxi);  zt = zb*dxi - dble(kc)

   ! Shift to 1-based Fortran indexing, clamp into [1,nx]
   ic = mod(ic, nx) + 1
   jc = mod(jc, nx) + 1
   kc = mod(kc, nx) + 1

   ! Periodic +1 neighbour
   icp = ic + 1;  if (icp > nx) icp = 1
   jcp = jc + 1;  if (jcp > nx) jcp = 1
   kcp = kc + 1;  if (kcp > nx) kcp = 1

   ! Trilinear weights
   c1 = (1.d0-xt)*(1.d0-yt)*(1.d0-zt)
   c2 =       xt *(1.d0-yt)*(1.d0-zt)
   c3 = (1.d0-xt)*      yt *(1.d0-zt)
   c4 =       xt *      yt *(1.d0-zt)
   c5 = (1.d0-xt)*(1.d0-yt)*      zt
   c6 =       xt *(1.d0-yt)*      zt
   c7 = (1.d0-xt)*      yt *      zt
   c8 =       xt *      yt *      zt

   ! --- Fluid velocity ---
   uf_at_p(m,1) = c1*u(ic,jc,kc)  + c2*u(icp,jc,kc)  + c3*u(ic,jcp,kc)  + c4*u(icp,jcp,kc) &
                + c5*u(ic,jc,kcp) + c6*u(icp,jc,kcp) + c7*u(ic,jcp,kcp) + c8*u(icp,jcp,kcp)
   uf_at_p(m,2) = c1*v(ic,jc,kc)  + c2*v(icp,jc,kc)  + c3*v(ic,jcp,kc)  + c4*v(icp,jcp,kc) &
                + c5*v(ic,jc,kcp) + c6*v(icp,jc,kcp) + c7*v(ic,jcp,kcp) + c8*v(icp,jcp,kcp)
   uf_at_p(m,3) = c1*w(ic,jc,kc)  + c2*w(icp,jc,kc)  + c3*w(ic,jcp,kc)  + c4*w(icp,jcp,kc) &
                + c5*w(ic,jc,kcp) + c6*w(icp,jc,kcp) + c7*w(ic,jcp,kcp) + c8*w(icp,jcp,kcp)

   ! --- Velocity gradient at the 8 cell corners (centred FD, O(dx^2)) ---
   do cz = 0, 1
   do cy = 0, 1
   do cx = 0, 1
      ix = ic + cx;  if (ix > nx) ix = 1
      iy = jc + cy;  if (iy > nx) iy = 1
      iz = kc + cz;  if (iz > nx) iz = 1
      ipp = ix + 1;  if (ipp > nx) ipp = 1;  imm = ix - 1;  if (imm < 1) imm = nx
      jpp = iy + 1;  if (jpp > nx) jpp = 1;  jmm = iy - 1;  if (jmm < 1) jmm = nx
      kpp = iz + 1;  if (kpp > nx) kpp = 1;  kmm = iz - 1;  if (kmm < 1) kmm = nx
      du_dx(cx+1,cy+1,cz+1) = 0.5d0*dxi*(u(ipp,iy,iz) - u(imm,iy,iz))
      du_dy(cx+1,cy+1,cz+1) = 0.5d0*dxi*(u(ix,jpp,iz) - u(ix,jmm,iz))
      du_dz(cx+1,cy+1,cz+1) = 0.5d0*dxi*(u(ix,iy,kpp) - u(ix,iy,kmm))
      dv_dx(cx+1,cy+1,cz+1) = 0.5d0*dxi*(v(ipp,iy,iz) - v(imm,iy,iz))
      dv_dy(cx+1,cy+1,cz+1) = 0.5d0*dxi*(v(ix,jpp,iz) - v(ix,jmm,iz))
      dv_dz(cx+1,cy+1,cz+1) = 0.5d0*dxi*(v(ix,iy,kpp) - v(ix,iy,kmm))
      dw_dx(cx+1,cy+1,cz+1) = 0.5d0*dxi*(w(ipp,iy,iz) - w(imm,iy,iz))
      dw_dy(cx+1,cy+1,cz+1) = 0.5d0*dxi*(w(ix,jpp,iz) - w(ix,jmm,iz))
      dw_dz(cx+1,cy+1,cz+1) = 0.5d0*dxi*(w(ix,iy,kpp) - w(ix,iy,kmm))
   enddo; enddo; enddo

   ! Trilinearly interpolate gradient components; layout: grad(1..9) = du/dx..dw/dz
   grad_at_p(m,1)=c1*du_dx(1,1,1)+c2*du_dx(2,1,1)+c3*du_dx(1,2,1)+c4*du_dx(2,2,1) &
                 +c5*du_dx(1,1,2)+c6*du_dx(2,1,2)+c7*du_dx(1,2,2)+c8*du_dx(2,2,2)
   grad_at_p(m,2)=c1*du_dy(1,1,1)+c2*du_dy(2,1,1)+c3*du_dy(1,2,1)+c4*du_dy(2,2,1) &
                 +c5*du_dy(1,1,2)+c6*du_dy(2,1,2)+c7*du_dy(1,2,2)+c8*du_dy(2,2,2)
   grad_at_p(m,3)=c1*du_dz(1,1,1)+c2*du_dz(2,1,1)+c3*du_dz(1,2,1)+c4*du_dz(2,2,1) &
                 +c5*du_dz(1,1,2)+c6*du_dz(2,1,2)+c7*du_dz(1,2,2)+c8*du_dz(2,2,2)
   grad_at_p(m,4)=c1*dv_dx(1,1,1)+c2*dv_dx(2,1,1)+c3*dv_dx(1,2,1)+c4*dv_dx(2,2,1) &
                 +c5*dv_dx(1,1,2)+c6*dv_dx(2,1,2)+c7*dv_dx(1,2,2)+c8*dv_dx(2,2,2)
   grad_at_p(m,5)=c1*dv_dy(1,1,1)+c2*dv_dy(2,1,1)+c3*dv_dy(1,2,1)+c4*dv_dy(2,2,1) &
                 +c5*dv_dy(1,1,2)+c6*dv_dy(2,1,2)+c7*dv_dy(1,2,2)+c8*dv_dy(2,2,2)
   grad_at_p(m,6)=c1*dv_dz(1,1,1)+c2*dv_dz(2,1,1)+c3*dv_dz(1,2,1)+c4*dv_dz(2,2,1) &
                 +c5*dv_dz(1,1,2)+c6*dv_dz(2,1,2)+c7*dv_dz(1,2,2)+c8*dv_dz(2,2,2)
   grad_at_p(m,7)=c1*dw_dx(1,1,1)+c2*dw_dx(2,1,1)+c3*dw_dx(1,2,1)+c4*dw_dx(2,2,1) &
                 +c5*dw_dx(1,1,2)+c6*dw_dx(2,1,2)+c7*dw_dx(1,2,2)+c8*dw_dx(2,2,2)
   grad_at_p(m,8)=c1*dw_dy(1,1,1)+c2*dw_dy(2,1,1)+c3*dw_dy(1,2,1)+c4*dw_dy(2,2,1) &
                 +c5*dw_dy(1,1,2)+c6*dw_dy(2,1,2)+c7*dw_dy(1,2,2)+c8*dw_dy(2,2,2)
   grad_at_p(m,9)=c1*dw_dz(1,1,1)+c2*dw_dz(2,1,1)+c3*dw_dz(1,2,1)+c4*dw_dz(2,2,1) &
                 +c5*dw_dz(1,1,2)+c6*dw_dz(2,1,2)+c7*dw_dz(1,2,2)+c8*dw_dz(2,2,2)

   ! Vorticity from the interpolated gradient
   ! vort(1)=dw/dy-dv/dz, vort(2)=du/dz-dw/dx, vort(3)=dv/dx-du/dy
   vort_at_p(m,1) = grad_at_p(m,8) - grad_at_p(m,6)
   vort_at_p(m,2) = grad_at_p(m,3) - grad_at_p(m,7)
   vort_at_p(m,3) = grad_at_p(m,4) - grad_at_p(m,2)
enddo
!$OMP END PARALLEL DO

end subroutine get_fiber_fluid_quantities


! ============================================================================
!  spread_fiber_reaction
!  Two-way coupling: spread the accumulated per-bead reaction impulses
!  (f_reaction) onto the Eulerian force field fxfib/fyfib/fzfib, which is added
!  to the NSE RHS in main.f90.  Uses the SAME trilinear stencil as
!  get_fiber_fluid_quantities, so the spread is the exact transpose of the
!  velocity interpolation => momentum is conserved (Σ f_reaction == Σ f?fib·dt·dx³).
!  Serial: many beads deposit into shared cells, so it is not parallelised.
! ============================================================================
subroutine spread_fiber_reaction(step)
use param,    only: dx, nx, dt
use velocity, only: fxfib, fyfib, fzfib
use fibers,   only: nbeads, xf, f_reaction
implicit none
integer, intent(in) :: step
integer          :: m, ic, jc, kc, icp, jcp, kcp
double precision :: xb, yb, zb, xt, yt, zt, dxi, sc, rx, ry, rz
double precision :: c1,c2,c3,c4,c5,c6,c7,c8
double precision :: pin(3), pres(3)

fxfib = 0.d0;  fyfib = 0.d0;  fzfib = 0.d0
dxi = 1.d0/dx
! reaction impulse [mass·vel] accumulated over the fluid step -> force density
! [accel] for the NSE RHS:  Δu = dt·f_fib = spread(impulse)/dx³  ⇒  /(dt·dx³)
sc  = 1.d0 / (dt * dx*dx*dx)

do m = 1, nbeads
   xb = xf(m,1);  yb = xf(m,2);  zb = xf(m,3)
   ic = floor(xb*dxi);  xt = xb*dxi - dble(ic)
   jc = floor(yb*dxi);  yt = yb*dxi - dble(jc)
   kc = floor(zb*dxi);  zt = zb*dxi - dble(kc)
   ic = mod(ic,nx)+1;  jc = mod(jc,nx)+1;  kc = mod(kc,nx)+1
   icp = ic+1;  if (icp>nx) icp=1
   jcp = jc+1;  if (jcp>nx) jcp=1
   kcp = kc+1;  if (kcp>nx) kcp=1
   c1=(1.d0-xt)*(1.d0-yt)*(1.d0-zt);  c2=xt*(1.d0-yt)*(1.d0-zt)
   c3=(1.d0-xt)*yt*(1.d0-zt);          c4=xt*yt*(1.d0-zt)
   c5=(1.d0-xt)*(1.d0-yt)*zt;          c6=xt*(1.d0-yt)*zt
   c7=(1.d0-xt)*yt*zt;                  c8=xt*yt*zt

   rx = f_reaction(m,1)*sc;  ry = f_reaction(m,2)*sc;  rz = f_reaction(m,3)*sc

   fxfib(ic ,jc ,kc ) = fxfib(ic ,jc ,kc ) + c1*rx
   fxfib(icp,jc ,kc ) = fxfib(icp,jc ,kc ) + c2*rx
   fxfib(ic ,jcp,kc ) = fxfib(ic ,jcp,kc ) + c3*rx
   fxfib(icp,jcp,kc ) = fxfib(icp,jcp,kc ) + c4*rx
   fxfib(ic ,jc ,kcp) = fxfib(ic ,jc ,kcp) + c5*rx
   fxfib(icp,jc ,kcp) = fxfib(icp,jc ,kcp) + c6*rx
   fxfib(ic ,jcp,kcp) = fxfib(ic ,jcp,kcp) + c7*rx
   fxfib(icp,jcp,kcp) = fxfib(icp,jcp,kcp) + c8*rx

   fyfib(ic ,jc ,kc ) = fyfib(ic ,jc ,kc ) + c1*ry
   fyfib(icp,jc ,kc ) = fyfib(icp,jc ,kc ) + c2*ry
   fyfib(ic ,jcp,kc ) = fyfib(ic ,jcp,kc ) + c3*ry
   fyfib(icp,jcp,kc ) = fyfib(icp,jcp,kc ) + c4*ry
   fyfib(ic ,jc ,kcp) = fyfib(ic ,jc ,kcp) + c5*ry
   fyfib(icp,jc ,kcp) = fyfib(icp,jc ,kcp) + c6*ry
   fyfib(ic ,jcp,kcp) = fyfib(ic ,jcp,kcp) + c7*ry
   fyfib(icp,jcp,kcp) = fyfib(icp,jcp,kcp) + c8*ry

   fzfib(ic ,jc ,kc ) = fzfib(ic ,jc ,kc ) + c1*rz
   fzfib(icp,jc ,kc ) = fzfib(icp,jc ,kc ) + c2*rz
   fzfib(ic ,jcp,kc ) = fzfib(ic ,jcp,kc ) + c3*rz
   fzfib(icp,jcp,kc ) = fzfib(icp,jcp,kc ) + c4*rz
   fzfib(ic ,jc ,kcp) = fzfib(ic ,jc ,kcp) + c5*rz
   fzfib(icp,jc ,kcp) = fzfib(icp,jc ,kcp) + c6*rz
   fzfib(ic ,jcp,kcp) = fzfib(ic ,jcp,kcp) + c7*rz
   fzfib(icp,jcp,kcp) = fzfib(icp,jcp,kcp) + c8*rz
enddo

! ── Validation diagnostic (every 100 steps): momentum conservation ───────────
! Spread is the exact interp transpose, so the residual must be ~machine zero.
if (mod(step,100) == 0) then
   pin(1) = sum(f_reaction(:,1));  pin(2) = sum(f_reaction(:,2));  pin(3) = sum(f_reaction(:,3))
   pres(1) = sum(fxfib)*dt*dx*dx*dx - pin(1)
   pres(2) = sum(fyfib)*dt*dx*dx*dx - pin(2)
   pres(3) = sum(fzfib)*dt*dx*dx*dx - pin(3)
   write(*,'(a,i8,a,3es11.3,a,3es9.1)') ' [2way] step ', step, &
        '  fiber->fluid impulse =', pin, '   conservation residual =', pres
endif
end subroutine spread_fiber_reaction


! ============================================================================
!  advance_fibers
!  Main time-step loop over all fibers.  Called once per fluid time step
!  (possibly with subiterations).  Parallelised over fibers with OMP DO.
! ============================================================================
subroutine advance_fibers(step)
use param
use fibers
use energy_budget_mod, only: eb_ncoll   ! aggregate energy-budget collision counter
implicit none

integer, intent(in) :: step
integer :: i, subst, j, k, ibead, jb, mm
double precision :: bead_mass, bead_mass_lb
integer :: cid, nbr_cid, ixc, iyc, izc, ix, iy, iz, dox, doy, doz
integer :: rr, ss, fid_rr, fid_ss, rid_rr, rid_ss, tmpcid, local_cnt
double precision :: dxxx, dyyy, dzzz, r_cut, cell_size
double precision :: cell_size_x, cell_size_y, cell_size_z
integer :: n_cell_x, n_cell_y, n_cell_z, cell_count
integer, allocatable :: head(:), next_b(:), alreadycoll(:)
integer :: total_coll
integer, dimension(nfibers) :: per_fib_coll
double precision, dimension(nfibers*lambda,3) :: up_add, omgp_add
double precision, dimension(lambda,3) :: fvxx0_add, fvxx1_add, fomg0_add, fomg1_add
integer :: flag_coll_this


! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! for friction
double precision, dimension(3,lambda,nfibers) :: ft_add,torque_add, fttrq_add  ! to be same dimension as X_0
double precision, dimension(3) :: ft1_add,ft0_add,torque_t0_add,torque_t1_add
double precision, dimension(3) :: fttrq0_add,fttrq1_add
integer :: stat_cnt,flag_f

! to understand if this is tooo big ??????????????????????????????????????????????
!integer, dimension(nfibers*lambda,nfibers*lambda) :: stat_flag
! NB this last one useful only in case of static friction!
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!for multifiber
double precision :: bead_mass_type(max_fiber_types)
integer :: kf,jf,fr,fs
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 


! Gears working arrays: one slice per fiber in each thread
double precision, allocatable :: ORIENTm_0(:,:,:,:)  ! (3,3,lambda,nfibers)
double precision, allocatable :: X_0(:,:,:)           ! (3,lambda,nfibers)

! ── Rigid slender-body model: dispatch to the dedicated integrator ───────────
! The rigid rod reuses all the same bead arrays (beads = quadrature points), so
! everything downstream (I/O, restart, coupling spread, viz) is unchanged.
if (fiber_model == 1) then
   call advance_fibers_rigid(step)
   return
endif

allocate(ORIENTm_0(3,3,lambda,nfibers))
allocate(X_0(3,lambda,nfibers))

! Cell list parameters
r_cut      = r_cut_coef * arl * a_par   ! r_cut_coef set in input.dat (&fibers_input)
! Cell size = interaction cutoff: the smallest cell for which the 3x3x3
! neighbour stencil is guaranteed to contain every partner within r_cut.
! (The old value 100*arl*a_par gave only ~3 cells/dim, so the stencil spanned
!  the whole box and the search degenerated to all-pairs, O(N^2).)
!
! CORRECTNESS: the actual cell size must be >= r_cut, else a contact pair within
! r_cut can straddle a 2-cell gap and be MISSED -> interpenetration -> a large
! catch-up collision impulse that breaks inextensibility.  Use FLOOR so that
! cell_size = lx/n_cell >= r_cut (ceiling would give lx/n_cell <= r_cut, which is
! exactly the wrong side of the bound).  max(3,..) keeps the wrapped stencil
! covering the whole box in the degenerate few-cell case.
cell_size  = r_cut
n_cell_x   = max(3, int(floor(lx / cell_size)))
n_cell_y   = n_cell_x
n_cell_z   = n_cell_x
cell_count = n_cell_x * n_cell_y * n_cell_z
cell_size_x = lx / n_cell_x
cell_size_y = lx / n_cell_y
cell_size_z = lx / n_cell_z

allocate(head(cell_count))
allocate(next_b(nbeads))
allocate(alreadycoll(nbeads))

! Two-way coupling: bead inertial mass (sphere of radius a_par, density dens_part
! with rho_f=1 — consistent with the Stokes drag rate in calculate_KQ_f).  Reset
! the reaction accumulator for this fluid step.

! mmmmmmmmmmmmmmmmmmmmmmmmmmm
! calculate mass
if (fiber_types_flag == 1) then
   !do kf = 1,nfibers 
      !jf = floor(min((kf-1)/type_idx + 1, n_fiber_types))
      !bead_mass_type(jf) = (4.d0/3.d0)*pi*arl_type(jf)* &
      !                     a_par_type(jf)**3.d0*dens_part_type(jf)
   !enddo   
   do jf = 1, n_fiber_types
      bead_mass_type(jf) = (4.d0/3.d0)*pi * &
           arl_type(jf)*a_par_type(jf)**3.d0 * dens_part_type(jf)
   enddo
elseif (bead_types_flag == 1) then

   bead_mass_lb = (4.d0/3.d0)*pi*arl_lb*a_par_lb**3.d0*dens_part
   bead_mass = (4.d0/3.d0)*pi*arl*a_par**3.d0*dens_part

else
   bead_mass = (4.d0/3.d0)*pi*arl*a_par**3.d0*dens_part
endif
if (coupling_flag == 1) f_reaction = 0.d0


! ======= starting main routine ================
do subst = 1, subiterations

   total_coll    = 0
   per_fib_coll  = 0
   up_add        = 0.0d0
   omgp_add      = 0.0d0
   X_constraint  = 0.0d0   ! constraint impulse – filled in constraint_flexfibre_f
   torque_bend   = 0.0d0   ! bending torque    – filled in calculate_RS_f
   coll_partner  = 0        ! collision partner – filled in collision CRITICAL block
   ! initializing the frcition terms
   if (subst == 1) then
      ft_add        = 0.0d0
      fttrq_add     = 0.0d0
   endif
   ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
   ! initialize the static counter
   !stat_flag = 0
   ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

   ! ── Gears dynamics (OMP parallel over fibers) ────────────────────────────
   !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1) &
   !$OMP   SHARED(nfibers, lambda, step, subst, nstart, &
   !$OMP          kzero, qzero, rzero, szero, ORIENTm_0, X_0, &
   !$OMP          xf, angles_f, uf, omgf, ff, &
   !$OMP          uf_at_p, vort_at_p, grad_at_p, &
   !$OMP          vp_o, angl_o, lx, dt_part, a_par, arl, &
   !$OMP          dens_part, E_bend, curv_eq, &
   !$OMP          a_par_type, arl_type, dens_part_type, E_bend_type, &
   !$OMP          curv_eq_type, Idm, bead_mass_type, fiber_type, &
   !$OMP          geff_type, fiber_types_flag,type_idx, &
   !$OMP          X_a, Y_a, X_c, Y_c, Y_h, X_a_type, Y_a_type, &
   !$OMP          X_c_type, Y_c_type, Y_h_type, &
   !$OMP          bead_on_interface, n_iface, phi_at_p, &
   !$OMP          coupling_flag, f_reaction, bead_mass, &
   !$OMP          gravity_flag, geff, ft_add,torque_add, &
   !$OMP          fttrq_add,ft1_add,ft0_add, &
   !$OMP          fttrq0_add,fttrq1_add ) &
   !$OMP   PRIVATE(jb, mm, jf )
   do i = 1, nfibers

      call tensor_KQRS_f( &
           kzero(:, 1+(i-1)*lambda:i*lambda), &
           qzero(:,:,1+(i-1)*lambda:i*lambda), &
           rzero(:, 1+(i-1)*lambda:i*lambda), &
           szero(:,:,1+(i-1)*lambda:i*lambda), &
           ORIENTm_0(:,:,:,i), &
           i, step)

      ! Two-way coupling: accumulate the drag reaction this fiber imparts to the
      ! fluid.  kzero (=K) is the drag-relaxed velocity; uf is still the pre-update
      ! value (vv), so mass*(uf - kzero) is the reaction impulse onto the fluid.
      ! Only the external (drag) force is spread; the internal Gears constraint and
      ! bending forces exchange no momentum with the fluid.
      !
      ! With gravity on, kzero also absorbed the +dt_part*geff body-force impulse
      ! (see calculate_KQ_f).  Gravity is NOT a hydrodynamic reaction, so add it
      ! back here: the fluid must feel only the drag part.  Physically this makes a
      ! fiber settling at terminal velocity (kzero=uf) transmit exactly its weight
      ! (bead_mass*geff) to the fluid.  geff=0 unless gravity_flag=1 -> no-op off.
      
      
      ! mmmmmmmm added multifiber force mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

      if (fiber_types_flag == 1) then 

         jf = fiber_type(i)

         if (coupling_flag == 1) then
            do jb = 1, lambda
               mm = 1 + (i-1)*lambda + (jb-1)
               f_reaction(mm,1) = f_reaction(mm,1) + bead_mass_type(jf)*(uf(mm,1) - kzero(1,mm))  ! to modify k0 missing so NB NB 
               f_reaction(mm,2) = f_reaction(mm,2) + bead_mass_type(jf)*(uf(mm,2) - kzero(2,mm))
               f_reaction(mm,3) = f_reaction(mm,3) + bead_mass_type(jf)*(uf(mm,3) - kzero(3,mm))
            enddo
         if (gravity_flag == 1) then
               do jb = 1, lambda
                  mm = 1 + (i-1)*lambda + (jb-1)
                  f_reaction(mm,1) = f_reaction(mm,1) + bead_mass_type(jf)*geff_type(jf,1)*dt_part
                  f_reaction(mm,2) = f_reaction(mm,2) + bead_mass_type(jf)*geff_type(jf,2)*dt_part
                  f_reaction(mm,3) = f_reaction(mm,3) + bead_mass_type(jf)*geff_type(jf,3)*dt_part
               enddo
         endif
         endif
      
      elseif (bead_types_flag == 1) then
         if (coupling_flag == 1) then
            do jb = 1, lambda
               if (jb == lambda) then
                  mm = 1 + (i-1)*lambda + (jb-1)
                  f_reaction(mm,1) = f_reaction(mm,1) + bead_mass_lb*(uf(mm,1) - kzero(1,mm))
                  f_reaction(mm,2) = f_reaction(mm,2) + bead_mass_lb*(uf(mm,2) - kzero(2,mm))
                  f_reaction(mm,3) = f_reaction(mm,3) + bead_mass_lb*(uf(mm,3) - kzero(3,mm))
               else
                  mm = 1 + (i-1)*lambda + (jb-1)
                  f_reaction(mm,1) = f_reaction(mm,1) + bead_mass*(uf(mm,1) - kzero(1,mm))
                  f_reaction(mm,2) = f_reaction(mm,2) + bead_mass*(uf(mm,2) - kzero(2,mm))
                  f_reaction(mm,3) = f_reaction(mm,3) + bead_mass*(uf(mm,3) - kzero(3,mm))
               endif
            enddo
         if (gravity_flag == 1) then
               do jb = 1, lambda
                  if (jb == lambda) then
                     mm = 1 + (i-1)*lambda + (jb-1)
                     f_reaction(mm,1) = f_reaction(mm,1) + bead_mass_lb*geff(1)*dt_part
                     f_reaction(mm,2) = f_reaction(mm,2) + bead_mass_lb*geff(2)*dt_part
                     f_reaction(mm,3) = f_reaction(mm,3) + bead_mass_lb*geff(3)*dt_part
                  else
                     mm = 1 + (i-1)*lambda + (jb-1)
                     f_reaction(mm,1) = f_reaction(mm,1) + bead_mass*geff(1)*dt_part
                     f_reaction(mm,2) = f_reaction(mm,2) + bead_mass*geff(2)*dt_part
                     f_reaction(mm,3) = f_reaction(mm,3) + bead_mass*geff(3)*dt_part
                  endif
               enddo
         endif
         endif

      else
         if (coupling_flag == 1) then
            do jb = 1, lambda
               mm = 1 + (i-1)*lambda + (jb-1)
               f_reaction(mm,1) = f_reaction(mm,1) + bead_mass*(uf(mm,1) - kzero(1,mm))
               f_reaction(mm,2) = f_reaction(mm,2) + bead_mass*(uf(mm,2) - kzero(2,mm))
               f_reaction(mm,3) = f_reaction(mm,3) + bead_mass*(uf(mm,3) - kzero(3,mm))
            enddo
         if (gravity_flag == 1) then
               do jb = 1, lambda
                  mm = 1 + (i-1)*lambda + (jb-1)
                  f_reaction(mm,1) = f_reaction(mm,1) + bead_mass*geff(1)*dt_part
                  f_reaction(mm,2) = f_reaction(mm,2) + bead_mass*geff(2)*dt_part
                  f_reaction(mm,3) = f_reaction(mm,3) + bead_mass*geff(3)*dt_part
               enddo
         endif
         endif
      endif

      ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      ! not done for asymmetric fibers ------- start from here to modify -----------
      ! to understand how to write constraints since they are from ind = 1 to lambda - 1 
      call constraint_flexfibre_f( &
           X_0(:,:,i), &
           kzero(:, 1+(i-1)*lambda:i*lambda), &
           qzero(:,:,1+(i-1)*lambda:i*lambda), &
           rzero(:, 1+(i-1)*lambda:i*lambda), &
           szero(:,:,1+(i-1)*lambda:i*lambda), &
           ORIENTm_0(:,:,:,i), i, &
           ft_add(:,:,i), fttrq_add(:,:,i)) ! added for friction

      call dynamics_flexible_fiber_f( &
           X_0(:,:,i), &
           kzero(:, 1+(i-1)*lambda:i*lambda), &
           qzero(:,:,1+(i-1)*lambda:i*lambda), &
           rzero(:, 1+(i-1)*lambda:i*lambda), &
           szero(:,:,1+(i-1)*lambda:i*lambda), &
           ORIENTm_0(:,:,:,i), &
           i, ft_add(:,:,i), fttrq_add(:,:,i))

      ! Interface holding is now applied as a drag-target attraction inside
      ! tensor_KQRS_f (re-projected through the constraint solve), which holds the
      ! fiber on the drop without breaking inextensibility.  The old post-solve
      ! velocity overwrite (apply_interface_constraints_fiber) is superseded and
      ! no longer called — it could not hold the bead without risking a crash.

      call integrate_flexible_fiber_f(i, step, subst)

      call quality_check_f( &
           xf(1+(i-1)*lambda:i*lambda,:), &
           angles_f(1+(i-1)*lambda:i*lambda,:), &
           uf(1+(i-1)*lambda:i*lambda,:), &
           omgf(1+(i-1)*lambda:i*lambda,:), i)

   enddo
   !$OMP END PARALLEL DO

   ! ── Fiber-fiber collision (skipped when collision_flag=0) ────────────────
   if (collision_flag == 1) then

   ! ── Build cell list for fiber-fiber collision detection ──────────────────
   !$OMP SINGLE
   head        = -1
   next_b      = -1
   alreadycoll = 0
   !initializzation also of friction
   flag_f      = 0
   do i = 1, nfibers
      do j = 1, lambda
         ibead = j + (i-1)*lambda
         ixc = int(floor(xf(ibead,1)/cell_size_x)) + 1
         iyc = int(floor(xf(ibead,2)/cell_size_y)) + 1
         izc = int(floor(xf(ibead,3)/cell_size_z)) + 1
         ixc = max(1, min(n_cell_x, ixc))
         iyc = max(1, min(n_cell_y, iyc))
         izc = max(1, min(n_cell_z, izc))
         cid = ixc + (iyc-1)*n_cell_x + (izc-1)*n_cell_x*n_cell_y
         next_b(ibead) = head(cid)
         head(cid)     = ibead
      enddo
   enddo
   !$OMP END SINGLE

  ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  ! added to model friction and ft reinizializzation
  ft_add = 0.0d0
  torque_add = 0.0d0
  fttrq_add = 0.0d0
  !initialize multifiber to avoid problems
  fs = 0
  fr = 0 ! NB added also in OMP private
  ! so we reinitialize our vector before the next contact calculation righttt????


   ! ── Fiber-fiber collision detection (OMP DO over cells) ──────────────────
   !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC,64) &
   !$OMP   SHARED(n_cell_x, n_cell_y, n_cell_z, cell_count, &
   !$OMP          head, next_b, alreadycoll, xf, angles_f, uf, omgf, &
   !$OMP          lx, r_cut, cell_size_x, cell_size_y, cell_size_z, &
   !$OMP          lambda, a_par, arl, Idm, up_add, omgp_add, &
   !$OMP          a_par_type, arl_type,fiber_type,  &
   !$OMP          per_fib_coll, total_coll, coll_partner, &
   !$OMP          ft_add,torque_add, fiber_types_flag, &
   !$OMP          fttrq_add, type_idx) &
   !$OMP   PRIVATE(cid, nbr_cid, ixc, iyc, izc, ix, iy, iz, dox, doy, doz,fr,fs, &
   !$OMP           rr, ss, fid_rr, fid_ss, rid_rr, rid_ss, tmpcid, local_cnt, &
   !$OMP           dxxx, dyyy, dzzz, fvxx0_add, fvxx1_add, fomg0_add, fomg1_add &
   !$OMP           ft1_add, ft0_add, fttrq0_add,fttrq1_add, &
   !$OMP           torque_t0_add,torque_t1_add, flag_f, stat_cnt)
   do cid = 1, cell_count
      izc   = (cid-1) / (n_cell_x*n_cell_y) + 1
      tmpcid = (cid-1) - (izc-1)*n_cell_x*n_cell_y
      iyc   = tmpcid / n_cell_x + 1
      ixc   = mod(tmpcid, n_cell_x) + 1

      do doz = -1, 1
         iz = izc + doz
         if (iz < 1 .or. iz > n_cell_z) cycle !to check if this is ok or wrong
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do doy = -1, 1
            iy = iyc + doy
            if (iy < 1) iy = iy + n_cell_y
            if (iy > n_cell_y) iy = iy - n_cell_y
            do dox = -1, 1
               ix = ixc + dox
               if (ix < 1) ix = ix + n_cell_x
               if (ix > n_cell_x) ix = ix - n_cell_x
               nbr_cid = ix + (iy-1)*n_cell_x + (iz-1)*n_cell_x*n_cell_y
               if (nbr_cid < cid) cycle

               rr = head(cid)
               do while (rr /= -1)
                  ss = head(nbr_cid)
                  do while (ss /= -1)
                     fid_rr = (rr-1)/lambda + 1
                     fid_ss = (ss-1)/lambda + 1
                     rid_rr = rr - (fid_rr-1)*lambda
                     rid_ss = ss - (fid_ss-1)*lambda
                     
                     ! mmm added for multifiber case - to get grp index of the two mmmmmmmm
                     
                     if (fiber_types_flag == 1) then
                     fr = fiber_type(fid_rr)
                     fs = fiber_type(fid_ss)
                     endif

                     ! mmmmmmmmmmm

                     if (nbr_cid == cid .and. ss <= rr) then
                        ! half-neighbour rule, skip
                     elseif ((abs(ss-rr) <= 2) .and. fid_rr == fid_ss) then
                        ! adjacent beads on same fiber, skip
                     else
                        dxxx = xf(rr,1) - xf(ss,1)
                        if (dxxx > 0.5d0*lx) dxxx = dxxx - lx
                        if (dxxx <= -0.5d0*lx) dxxx = dxxx + lx
                        dyyy = xf(rr,2) - xf(ss,2)
                        if (dyyy > 0.5d0*lx) dyyy = dyyy - lx
                        if (dyyy <= -0.5d0*lx) dyyy = dyyy + lx
                        dzzz = xf(rr,3) - xf(ss,3)
                        if (dzzz > 0.5d0*lx) dzzz = dzzz - lx
                        if (dzzz <= -0.5d0*lx) dzzz = dzzz + lx

                        if (dxxx*dxxx+dyyy*dyyy+dzzz*dzzz < r_cut*r_cut) then
                           local_cnt    = 0
                           !reinitialize also the friction collision 
                           flag_f = 0
                           fvxx0_add = 0.d0; fvxx1_add = 0.d0
                           fomg0_add = 0.d0; fomg1_add = 0.d0
                           !initializzation also of the friction variables
                           ft0_add = 0.d0; ft1_add = 0.d0
                           torque_t0_add = 0.d0
                           torque_t1_add = 0.d0
                           fttrq0_add = 0.d0; fttrq1_add = 0.d0
                           stat_cnt = 0

                           call fibre_collision_check_f( &
                                xf(rr,:), xf(ss,:), &
                                angles_f(rr,:), angles_f(ss,:), local_cnt, &
                                fid_rr, fid_ss, rid_rr, rid_ss, &
                                uf(1+(fid_rr-1)*lambda:fid_rr*lambda,:), &
                                uf(1+(fid_ss-1)*lambda:fid_ss*lambda,:), &
                                angles_f(1+(fid_rr-1)*lambda:fid_rr*lambda,:), &
                                angles_f(1+(fid_ss-1)*lambda:fid_ss*lambda,:), &
                                omgf(1+(fid_rr-1)*lambda:fid_rr*lambda,:), &
                                omgf(1+(fid_ss-1)*lambda:fid_ss*lambda,:), &
                                fvxx0_add, fvxx1_add, fomg0_add, fomg1_add, &
                                alreadycoll(rr), alreadycoll(ss), & 
                                ft0_add, ft1_add,torque_t0_add,torque_t1_add, &
                                stat_cnt,flag_f, fttrq0_add,fttrq1_add, fr,fs ) ! last terms for friction and multifiber

                           ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
                           ! adding also tangential forces on contact point 
                  
                           if (flag_f .gt. 0) then

                              !$OMP CRITICAL 
                     
                              ft_add(:,rid_rr,fid_rr) = ft_add(:,rid_rr,fid_rr) + ft0_add
                              ft_add(:,rid_ss,fid_ss) = ft_add(:,rid_ss,fid_ss) + ft1_add

                              ! torques NB maybe they should be non-dim. in function of the leverage of X_0
                              torque_add(:,rid_rr,fid_rr) = torque_add(:,rid_rr,fid_rr) + torque_t0_add
                              torque_add(:,rid_ss,fid_ss) = torque_add(:,rid_ss,fid_ss) + torque_t1_add
                              ! NBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB da vedere se giustooooooooooooooo
                     
                              ! corrected vector for the force eq. torque terms
                              fttrq_add(:,rid_rr,fid_rr) = fttrq_add(:,rid_rr,fid_rr) + fttrq0_add
                              fttrq_add(:,rid_ss,fid_ss) = fttrq_add(:,rid_ss,fid_ss) + fttrq1_add

                              ! vector and criteria to store all the static contacts                      
                              !if (stat_cnt .gt. 0) then 
                                 !stat_flag(rr, ss) = 1 
                                 !stat_flag(ss,rr) = 1
                              !endif
                              !$OMP END CRITICAL
                           endif

                  
                           ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

                           if (local_cnt > 0) then
                              !$OMP CRITICAL
                              ! One collision per bead per step.  The gate inside
                              ! fibre_collision_check_f is only a pre-filter; re-
                              ! check under the lock so a bead already impulsed by
                              ! another thread this step is not hit again.  Summing
                              ! a second impulse against the same pre-collision
                              ! velocity would inject energy (both sides are
                              ! discarded together, so momentum stays conserved).
                              ! For the flexible chain this is the COMPLETE policy
                              ! (every bead is a genuine distinct contact sphere, so
                              ! forbidding a bead from firing twice removes the whole
                              ! double-count; nearest-only-per-pair would wrongly
                              ! drop real multi-point contacts between curved fibers).
                              ! Validated in validation/collision/ (COLLISION_VALIDATION.md
                              ! sec 6.2, tests R10/R12).  Contrast the rigid path, where
                              ! one-per-bead is only the momentum-exact interim.
                              if (alreadycoll(rr)==0 .and. alreadycoll(ss)==0) then
                                 up_add(1+(fid_rr-1)*lambda:fid_rr*lambda,:) = &
                                    up_add(1+(fid_rr-1)*lambda:fid_rr*lambda,:) + fvxx0_add
                                 up_add(1+(fid_ss-1)*lambda:fid_ss*lambda,:) = &
                                    up_add(1+(fid_ss-1)*lambda:fid_ss*lambda,:) + fvxx1_add
                                 omgp_add(1+(fid_rr-1)*lambda:fid_rr*lambda,:) = &
                                    omgp_add(1+(fid_rr-1)*lambda:fid_rr*lambda,:) + fomg0_add
                                 omgp_add(1+(fid_ss-1)*lambda:fid_ss*lambda,:) = &
                                    omgp_add(1+(fid_ss-1)*lambda:fid_ss*lambda,:) + fomg1_add
                                 ! Record which fiber each colliding bead touched (col 38)
                                 coll_partner(rr) = fid_ss
                                 coll_partner(ss) = fid_rr
                                 alreadycoll(rr) = 1
                                 alreadycoll(ss) = 1
                                 total_coll = total_coll + local_cnt
                                 per_fib_coll(min(fid_rr,fid_ss)) = &
                                    per_fib_coll(min(fid_rr,fid_ss)) + local_cnt
                              endif
                              !$OMP END CRITICAL
                           endif
                        endif
                     endif
                     ss = next_b(ss)
                  enddo
                  rr = next_b(rr)
               enddo
            enddo
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO

   write(*,'(a,i6,a,i6,a,i8)') ' step=', step, ' sub=', subst, &
        '  total fiber-fiber collisions=', total_coll

   !write (*, '(a,f10.4)') 'total friction force added', ft_add

   ! Aggregate energy-budget diagnostic: accumulate collisions across sub-steps
   ! and steps (serial region — no OMP sharing); reset when a budget line is
   ! written in energy_budget_write.
   eb_ncoll = eb_ncoll + total_coll

   ! ── Apply collision velocity corrections ─────────────────────────────────
   !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
   !$OMP   SHARED(nfibers, lambda, up_add, omgp_add, uf, omgf, vp_o, angl_o, angles_f) &
   !$OMP   PRIVATE(i, j, k, flag_coll_this)
   do i = 1, nfibers
      flag_coll_this = 0
      do j = 1, lambda
         k = j + (i-1)*lambda
         if (any(up_add(k,:) /= 0.d0) .or. any(omgp_add(k,:) /= 0.d0)) &
            flag_coll_this = 1
      enddo
      if (flag_coll_this == 1) then
         uf(1+(i-1)*lambda:i*lambda,:) = &
            uf(1+(i-1)*lambda:i*lambda,:) + up_add(1+(i-1)*lambda:i*lambda,:)
         vp_o(1+(i-1)*lambda:i*lambda,:,1) = uf(1+(i-1)*lambda:i*lambda,:)
         vp_o(1+(i-1)*lambda:i*lambda,:,2) = vp_o(1+(i-1)*lambda:i*lambda,:,1)
         omgf(1+(i-1)*lambda:i*lambda,:) = &
            omgf(1+(i-1)*lambda:i*lambda,:) + omgp_add(1+(i-1)*lambda:i*lambda,:)
         do j = 1, lambda
            k = j + (i-1)*lambda
            angl_o(k,1,1) = angles_f(k,2)*omgf(k,3) - omgf(k,2)*angles_f(k,3)
            angl_o(k,2,1) = -(angles_f(k,1)*omgf(k,3) - omgf(k,1)*angles_f(k,3))
            angl_o(k,3,1) = angles_f(k,1)*omgf(k,2) - omgf(k,1)*angles_f(k,2)
            angl_o(k,:,2) = angl_o(k,:,1)
         enddo
      endif
   enddo
   !$OMP END PARALLEL DO

   end if   ! collision_flag

enddo  ! subiterations

if (nfibers == 1) call jeffery_analytics_f(xf(1:lambda,:), step, dble(step)*dt_part)

deallocate(ORIENTm_0, X_0, head, next_b, alreadycoll)
end subroutine advance_fibers


! ============================================================================
!  advance_fibers_rigid
!  Rigid slender-body rod integrator (selected by fiber_model=1).  Each fiber
!  is a RIGID rod with 6 DOF: centre-of-mass position Xbar, velocity Vbar,
!  orientation unit vector p, and transverse angular velocity omegaT (⟂ p).
!
!  The lambda beads are reused as Nq quadrature points sitting at
!      x_bead(n) = Xbar + s_quad(n)*p ,     sum_n s_quad(n) = 0,
!  so every bead-based facility downstream (fluid interpolation uf_at_p,
!  two-way coupling spread, 40-column I/O, restart, visualisation) works
!  unchanged.  The rigid DOF are RECONSTRUCTED from the bead arrays each call
!  (no new persistent state → restart is automatic):
!      p      = angles_f(anchor)                      (all beads share it)
!      Xbar   = xf(anchor) - s_quad(1)*p
!      Vbar   = mean_n uf(bead_n)                      (exact: sum s_quad = 0)
!      omegaT = p × [ (sum_n s_quad(n) uf(bead_n)) / S2_quad ]
!
!  Physics (per user pseudocode / Batchelor 1970, Jeffery 1922):
!    - sample the frozen flow along the rod → Ubar (mean) and Gvec=Gsum/S2
!      (least-squares axial velocity slope, i.e. dU/ds along the rod);
!    - resistance tensor  A = (1/tau_perp)(I-pp) + (1/tau_par)pp;
!    - translation, backward Euler:  Vbar = (I+dt A)^-1 (Vbar + dt A Ubar + dt geff);
!    - rotation, scalar-implicit:  omegaT=(omegaT+(dt/tau_rot)(p×Gvec))/(1+dt/tau_rot),
!      then reproject ⟂ p;
!    - orientation by exponential map:  p = p cosθ + (n̂×p) sinθ,  θ=|omegaT|dt;
!    - position:  Xbar += dt Vbar  (periodic wrap on the rebuilt beads).
!
!  Sub-iterates subiterations×dt_part = dt with the flow sample frozen over the
!  fluid step (same cadence as the flexible path).  Fiber-fiber collisions are
!  not applied in rigid mode (would require rod-level impulse resolution).
! ============================================================================
subroutine advance_fibers_rigid(step)
use param,  only: pi, lx, dt
use fibers
implicit none
integer, intent(in) :: step
integer :: i, n, mm, anchor, subst
double precision :: p(3), Ubar(3), Gvec(3), Vbar(3), omegaT(3), Wvec(3), Xbar(3)
double precision :: pp(3,3), AA(3,3), res(3,3), invres(3,3), rhs(3), temp(3)
double precision :: pxG(3), nhat(3), pnew(3), cx(3)
double precision :: rate, th, cth, sth, wmag, bead_mass, invS2
double precision, parameter :: tiny = 1.d-14
double precision, allocatable :: uf_start(:,:)   ! bead velocities at step start

! ── fiber-fiber collision (rigid two-body impulse) work variables ────────────
integer :: cid, nbr_cid, ixc, iyc, izc, ix, iy, iz, dox, doy, doz
integer :: rr, ss, fid_rr, fid_ss, rid_rr, rid_ss, tmpcid, local_cnt
integer :: n_cell_x, n_cell_y, n_cell_z, cell_count, total_coll
integer, allocatable :: head(:), next_b(:), alreadycoll(:)
double precision :: dxxx, dyyy, dzzz, r_cut, cell_size
double precision :: cell_size_x, cell_size_y, cell_size_z
double precision :: dVb0(3), dVb1(3), dwT0(3), dwT1(3)
double precision, allocatable :: dVbar_add(:,:), domgT_add(:,:)
integer, dimension(nfibers) :: per_fib_coll

! Bead inertial mass (= Mfib/lambda); same definition as the flexible path so
! the two-way momentum bookkeeping is identical.
bead_mass = (4.d0/3.d0)*pi*arl*a_par**3.d0*dens_part

! Two-way coupling: remember the pre-step bead velocities for the drag reaction.
if (coupling_flag == 1) then
   f_reaction = 0.d0
   allocate(uf_start(nbeads,3))
   uf_start = uf
endif

if (S2_quad > tiny) then
   invS2 = 1.d0/S2_quad
else
   invS2 = 0.d0        ! lambda=1 degenerate rod: no rotation / slope
endif

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1) &
!$OMP   SHARED(nfibers, lambda, subiterations, dt_part, lx, &
!$OMP          xf, angles_f, uf, omgf, uf_at_p, vp_o, angl_o, &
!$OMP          s_quad, invS2, Idm, geff, gravity_flag, &
!$OMP          tau_par, tau_perp, tau_rot) &
!$OMP   PRIVATE(i, n, mm, anchor, subst, p, Ubar, Gvec, Vbar, omegaT, Wvec, Xbar, &
!$OMP           pp, AA, res, invres, rhs, temp, pxG, nhat, pnew, cx, &
!$OMP           rate, th, cth, sth, wmag)
do i = 1, nfibers
   anchor = 1 + (i-1)*lambda

   ! ── Reconstruct the rigid DOF from the bead arrays ────────────────────────
   p = angles_f(anchor,:)
   wmag = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
   if (wmag > tiny) p = p/wmag
   Xbar = xf(anchor,:) - s_quad(1)*p          ! CM (unwrapped; rewrapped at end)

   Ubar = 0.d0;  Gvec = 0.d0;  Vbar = 0.d0;  Wvec = 0.d0
   do n = 1, lambda
      mm   = anchor + (n-1)
      Ubar = Ubar + uf_at_p(mm,:)                 ! mean flow along rod
      Gvec = Gvec + s_quad(n)*uf_at_p(mm,:)       ! axial slope numerator
      Vbar = Vbar + uf(mm,:)                       ! mean bead velocity = Vbar
      Wvec = Wvec + s_quad(n)*uf(mm,:)            ! transverse slope = S2*(omegaT×p)
   enddo
   Ubar = Ubar / dble(lambda)
   Vbar = Vbar / dble(lambda)
   Gvec = Gvec * invS2
   Wvec = Wvec * invS2
   ! omegaT = p × Wvec = p × (omegaT × p) → guaranteed perpendicular to p
   omegaT(1) = p(2)*Wvec(3) - p(3)*Wvec(2)
   omegaT(2) = p(3)*Wvec(1) - p(1)*Wvec(3)
   omegaT(3) = p(1)*Wvec(2) - p(2)*Wvec(1)

   rate = dt_part / tau_rot

   ! ── Sub-iterate the rigid dynamics (flow sample frozen) ───────────────────
   do subst = 1, subiterations

      ! resistance tensor  A = (1/tau_perp)(I - pp) + (1/tau_par) pp
      pp(1,1)=p(1)*p(1); pp(1,2)=p(1)*p(2); pp(1,3)=p(1)*p(3)
      pp(2,1)=p(2)*p(1); pp(2,2)=p(2)*p(2); pp(2,3)=p(2)*p(3)
      pp(3,1)=p(3)*p(1); pp(3,2)=p(3)*p(2); pp(3,3)=p(3)*p(3)
      AA = (1.d0/tau_perp)*(Idm - pp) + (1.d0/tau_par)*pp

      ! translation: backward Euler on drag (+ buoyancy-reduced gravity)
      res = Idm + dt_part*AA
      call INVERSE_MATRIX_ff(res, invres)
      call vec_product_f(AA, Ubar, temp)
      rhs = Vbar + dt_part*temp
      if (gravity_flag == 1) rhs = rhs + dt_part*geff
      call vec_product_f(invres, rhs, Vbar)

      ! rotation: scalar-implicit relaxation of omegaT toward p × Gvec
      pxG(1) = p(2)*Gvec(3) - p(3)*Gvec(2)
      pxG(2) = p(3)*Gvec(1) - p(1)*Gvec(3)
      pxG(3) = p(1)*Gvec(2) - p(2)*Gvec(1)
      omegaT = (omegaT + rate*pxG) / (1.d0 + rate)
      wmag   = omegaT(1)*p(1) + omegaT(2)*p(2) + omegaT(3)*p(3)
      omegaT = omegaT - wmag*p               ! reproject perpendicular to p

      ! orientation: exponential map (rotate p about omegaT by θ=|omegaT|dt)
      wmag = sqrt(omegaT(1)**2 + omegaT(2)**2 + omegaT(3)**2)
      if (wmag > tiny) then
         nhat = omegaT/wmag
         th   = wmag*dt_part
         cth  = cos(th);  sth = sin(th)
         cx(1) = nhat(2)*p(3) - nhat(3)*p(2)
         cx(2) = nhat(3)*p(1) - nhat(1)*p(3)
         cx(3) = nhat(1)*p(2) - nhat(2)*p(1)
         pnew  = p*cth + cx*sth
         wmag  = sqrt(pnew(1)**2 + pnew(2)**2 + pnew(3)**2)
         if (wmag > tiny) p = pnew/wmag
      endif

      ! position: advance CM
      Xbar = Xbar + dt_part*Vbar

   enddo

   ! ── Write the updated rigid state back into the bead arrays ───────────────
   ! v(n) = Vbar + omegaT × (s_quad(n) p) = Vbar + s_quad(n) (omegaT × p)
   cx(1) = omegaT(2)*p(3) - omegaT(3)*p(2)
   cx(2) = omegaT(3)*p(1) - omegaT(1)*p(3)
   cx(3) = omegaT(1)*p(2) - omegaT(2)*p(1)
   do n = 1, lambda
      mm = anchor + (n-1)
      xf(mm,:) = Xbar + s_quad(n)*p
      xf(mm,1) = xf(mm,1) - lx*floor(xf(mm,1)/lx)   ! periodic wrap
      xf(mm,2) = xf(mm,2) - lx*floor(xf(mm,2)/lx)
      xf(mm,3) = xf(mm,3) - lx*floor(xf(mm,3)/lx)
      angles_f(mm,:) = p
      uf(mm,:)  = Vbar + s_quad(n)*cx
      omgf(mm,:) = omegaT
      ! AB history is unused by the rigid integrator; keep it consistent so a
      ! restart or a later switch to the flexible model starts cleanly.
      vp_o(mm,:,1) = uf(mm,:);  vp_o(mm,:,2) = uf(mm,:)
      angl_o(mm,:,1) = cx;      angl_o(mm,:,2) = cx
   enddo

enddo
!$OMP END PARALLEL DO

! ── Two-way coupling: spread the rod's drag reaction, excluding gravity ──────
! Reaction impulse on the fluid = -(drag impulse on rod) = bead_mass*(u_old-u_new)
! + weight impulse.  The rotational part sums to zero over each rod, so the NET
! transmitted momentum is the translational drag; a rod settling at terminal
! velocity (u_new=u_old) transmits exactly its weight bead_mass*geff.  geff=0
! unless gravity_flag=1, so this reduces to the validated scheme when off.
if (coupling_flag == 1) then
   do mm = 1, nbeads
      f_reaction(mm,1) = bead_mass*(uf_start(mm,1) - uf(mm,1))
      f_reaction(mm,2) = bead_mass*(uf_start(mm,2) - uf(mm,2))
      f_reaction(mm,3) = bead_mass*(uf_start(mm,3) - uf(mm,3))
   enddo
   if (gravity_flag == 1) then
      do mm = 1, nbeads
         f_reaction(mm,1) = f_reaction(mm,1) + bead_mass*geff(1)*dt
         f_reaction(mm,2) = f_reaction(mm,2) + bead_mass*geff(2)*dt
         f_reaction(mm,3) = f_reaction(mm,3) + bead_mass*geff(3)*dt
      enddo
   endif
   deallocate(uf_start)
endif

! ── Fiber-fiber collision (rigid two-body impulse) ───────────────────────────
! Resolved ONCE per fluid step, AFTER the drag update and AFTER f_reaction is
! set: a collision is internal fiber-fiber momentum exchange and must not
! contaminate the two-way fluid reaction (which carries drag only).  The rod
! moves << a_par per step, so a single post-drag pass is sufficient — there is
! no need to re-resolve inside the drag sub-iterations (contrast the flexible
! path, which resolves every sub-iteration).  Detection reuses the flexible
! cell-list + segment closest-approach geometry; only the impulse response is
! rigid (see rigid_collision_pair_f).  nfibers>1 is required for a partner.
if (collision_flag == 1 .and. nfibers > 1) then

   ! Cell list: cell = interaction cutoff so the 3x3x3 stencil covers r_cut
   ! (identical sizing to the flexible path; FLOOR keeps cell_size >= r_cut).
   r_cut       = r_cut_coef * arl * a_par
   cell_size   = r_cut
   n_cell_x    = max(3, int(floor(lx / cell_size)))
   n_cell_y    = n_cell_x
   n_cell_z    = n_cell_x
   cell_count  = n_cell_x * n_cell_y * n_cell_z
   cell_size_x = lx / n_cell_x
   cell_size_y = lx / n_cell_y
   cell_size_z = lx / n_cell_z

   allocate(head(cell_count), next_b(nbeads), alreadycoll(nbeads))
   allocate(dVbar_add(nfibers,3), domgT_add(nfibers,3))
   dVbar_add = 0.d0;  domgT_add = 0.d0;  alreadycoll = 0
   total_coll = 0;    per_fib_coll = 0;    coll_partner = 0

   ! Build the (serial) cell list over all beads.
   head = -1;  next_b = -1
   do i = 1, nfibers
      do n = 1, lambda
         mm  = n + (i-1)*lambda
         ixc = max(1, min(n_cell_x, int(floor(xf(mm,1)/cell_size_x)) + 1))
         iyc = max(1, min(n_cell_y, int(floor(xf(mm,2)/cell_size_y)) + 1))
         izc = max(1, min(n_cell_z, int(floor(xf(mm,3)/cell_size_z)) + 1))
         cid = ixc + (iyc-1)*n_cell_x + (izc-1)*n_cell_x*n_cell_y
         next_b(mm) = head(cid)
         head(cid)  = mm
      enddo
   enddo

   ! Detect contacts between beads on DIFFERENT rods; accumulate the resulting
   ! rigid-body increment per rod (a straight rod cannot self-collide, so all
   ! same-fiber pairs are skipped — no in-fiber solver for the rigid model).
   !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC,64) &
   !$OMP   SHARED(n_cell_x, n_cell_y, n_cell_z, cell_count, head, next_b, &
   !$OMP          xf, angles_f, uf, omgf, lx, r_cut, lambda, alreadycoll, &
   !$OMP          dVbar_add, domgT_add, per_fib_coll, total_coll, coll_partner) &
   !$OMP   PRIVATE(cid, nbr_cid, ixc, iyc, izc, ix, iy, iz, dox, doy, doz, &
   !$OMP           rr, ss, fid_rr, fid_ss, rid_rr, rid_ss, tmpcid, local_cnt, &
   !$OMP           dxxx, dyyy, dzzz, dVb0, dVb1, dwT0, dwT1)
   do cid = 1, cell_count
      izc    = (cid-1) / (n_cell_x*n_cell_y) + 1
      tmpcid = (cid-1) - (izc-1)*n_cell_x*n_cell_y
      iyc    = tmpcid / n_cell_x + 1
      ixc    = mod(tmpcid, n_cell_x) + 1

      do doz = -1, 1
         iz = izc + doz
         if (iz < 1 .or. iz > n_cell_z) cycle
         do doy = -1, 1
            iy = iyc + doy
            if (iy < 1) iy = iy + n_cell_y
            if (iy > n_cell_y) iy = iy - n_cell_y
            do dox = -1, 1
               ix = ixc + dox
               if (ix < 1) ix = ix + n_cell_x
               if (ix > n_cell_x) ix = ix - n_cell_x
               nbr_cid = ix + (iy-1)*n_cell_x + (iz-1)*n_cell_x*n_cell_y
               if (nbr_cid < cid) cycle

               rr = head(cid)
               do while (rr /= -1)
                  ss = head(nbr_cid)
                  do while (ss /= -1)
                     fid_rr = (rr-1)/lambda + 1
                     fid_ss = (ss-1)/lambda + 1
                     rid_rr = rr - (fid_rr-1)*lambda
                     rid_ss = ss - (fid_ss-1)*lambda
                     if (nbr_cid == cid .and. ss <= rr) then
                        ! half-neighbour rule, skip
                     elseif (fid_rr == fid_ss) then
                        ! same rigid rod: cannot self-collide, skip
                     else
                        dxxx = xf(rr,1) - xf(ss,1)
                        if (dxxx > 0.5d0*lx) dxxx = dxxx - lx
                        if (dxxx <= -0.5d0*lx) dxxx = dxxx + lx
                        dyyy = xf(rr,2) - xf(ss,2)
                        if (dyyy > 0.5d0*lx) dyyy = dyyy - lx
                        if (dyyy <= -0.5d0*lx) dyyy = dyyy + lx
                        dzzz = xf(rr,3) - xf(ss,3)
                        if (dzzz > 0.5d0*lx) dzzz = dzzz - lx
                        if (dzzz <= -0.5d0*lx) dzzz = dzzz + lx

                        ! One collision per bead per step: the alreadycoll pre-check
                        ! skips beads already impulsed this step (rigid_collision_pair_f
                        ! itself has no such gate, unlike the flexible path).
                        if (dxxx*dxxx+dyyy*dyyy+dzzz*dzzz < r_cut*r_cut .and. &
                            alreadycoll(rr)==0 .and. alreadycoll(ss)==0) then
                           local_cnt = 0
                           call rigid_collision_pair_f( &
                                xf(rr,:), xf(ss,:), angles_f(rr,:), angles_f(ss,:), &
                                rid_rr, rid_ss, uf(rr,:), uf(ss,:), omgf(rr,:), omgf(ss,:), &
                                local_cnt, dVb0, dVb1, dwT0, dwT1)

                           if (local_cnt > 0) then
                              !$OMP CRITICAL
                              ! Re-check under the lock so a bead impulsed by another
                              ! thread this step is not hit again (summing a second
                              ! impulse vs. the same pre-collision velocity injects
                              ! energy; discarding both sides keeps momentum conserved).
                              if (alreadycoll(rr)==0 .and. alreadycoll(ss)==0) then
                                 dVbar_add(fid_rr,:) = dVbar_add(fid_rr,:) + dVb0
                                 dVbar_add(fid_ss,:) = dVbar_add(fid_ss,:) + dVb1
                                 domgT_add(fid_rr,:) = domgT_add(fid_rr,:) + dwT0
                                 domgT_add(fid_ss,:) = domgT_add(fid_ss,:) + dwT1
                                 coll_partner(rr) = fid_ss
                                 coll_partner(ss) = fid_rr
                                 alreadycoll(rr) = 1
                                 alreadycoll(ss) = 1
                                 total_coll = total_coll + local_cnt
                                 per_fib_coll(min(fid_rr,fid_ss)) = &
                                    per_fib_coll(min(fid_rr,fid_ss)) + local_cnt
                              endif
                              !$OMP END CRITICAL
                           endif
                        endif
                     endif
                     ss = next_b(ss)
                  enddo
                  rr = next_b(rr)
               enddo
            enddo
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO

   if (total_coll > 0) write(*,'(a,i7,a,i8)') ' step=', step, &
        '  rigid fiber-fiber collisions=', total_coll

   ! Apply the accumulated impulse to each rod: bump (Vbar, omegaT), reproject
   ! omegaT perpendicular to p, then rewrite the bead velocity field.  The rigid
   ! velocity field is linear in s_quad, so mapping a rod's summed (Vbar, omegaT)
   ! increment back onto its beads is exact -- but that is a KINEMATIC statement
   ! only.
   !
   ! CONTACT-RESOLUTION POLICY (validated in validation/collision/, see
   ! COLLISION_VALIDATION.md sec 6.2; tests R10-R14).  Summing SIMULTANEOUS
   ! impulses against a shared pre-collision velocity injects energy because
   ! distinct contacts on one rod are coupled through its 6 DOF (R10: +740%).
   ! The validated policy is:
   !   flexible  -> one collision per bead        (each bead a real contact sphere)
   !   rigid     -> nearest contact per rod pair  when contacts are mostly binary
   !                crossings (dilute; R11 = machine-zero energy), ESCALATING to
   !                sequential impulses / projected Gauss-Seidel when concentration
   !                makes multi-rod contact common (R13/R14 = machine-zero; nearest-
   !                per-pair alone still injects +57.7% across pairs, R14).
   ! What is LIVE below is the momentum-exact INTERIM: the one-collision-per-bead
   ! guard (alreadycoll) stops a single bead being impulsed twice, cutting the
   ! artifact +740% -> -5.8% (R12) with momentum still exact.  A rod contacting
   ! several partners at DISTINCT beads still sums those impulses against the
   ! pre-collision velocity -- a first-order treatment, acceptable because the rod
   ! moves << a_par per step.  Replacing this accumulate-then-apply block with the
   ! sequential (Gauss-Seidel) sweep for the concentrated production regime is the
   ! remaining GPU-box step (must be compiled+bit-verified under nvfortran; the
   ! reference implementation is rigid_integration_sequential() in the harness).
   ! The AB history (vp_o/angl_o) is unused by the rigid integrator but kept
   ! consistent for restart / a later model switch.
   !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1) &
   !$OMP   SHARED(nfibers, lambda, dVbar_add, domgT_add, &
   !$OMP          angles_f, uf, omgf, vp_o, angl_o, s_quad) &
   !$OMP   PRIVATE(i, n, mm, anchor, p, cx, temp, wmag)
   do i = 1, nfibers
      if (all(dVbar_add(i,:) == 0.d0) .and. all(domgT_add(i,:) == 0.d0)) cycle
      anchor = 1 + (i-1)*lambda
      p = angles_f(anchor,:)
      ! reproject the angular-velocity increment perpendicular to p (the rod
      ! carries no tracked axial spin), then form its per-bead velocity slope
      ! cx = (Δomega_T) × p.
      temp = domgT_add(i,:)
      wmag = temp(1)*p(1) + temp(2)*p(2) + temp(3)*p(3)
      temp = temp - wmag*p
      cx(1) = temp(2)*p(3) - temp(3)*p(2)
      cx(2) = temp(3)*p(1) - temp(1)*p(3)
      cx(3) = temp(1)*p(2) - temp(2)*p(1)
      do n = 1, lambda
         mm = anchor + (n-1)
         uf(mm,:)   = uf(mm,:)  + dVbar_add(i,:) + s_quad(n)*cx
         omgf(mm,:) = omgf(mm,:) + temp
         vp_o(mm,:,1) = uf(mm,:);  vp_o(mm,:,2) = uf(mm,:)
         angl_o(mm,1,1) = omgf(mm,2)*angles_f(mm,3) - omgf(mm,3)*angles_f(mm,2)
         angl_o(mm,2,1) = omgf(mm,3)*angles_f(mm,1) - omgf(mm,1)*angles_f(mm,3)
         angl_o(mm,3,1) = omgf(mm,1)*angles_f(mm,2) - omgf(mm,2)*angles_f(mm,1)
         angl_o(mm,:,2) = angl_o(mm,:,1)
      enddo
   enddo
   !$OMP END PARALLEL DO

   deallocate(head, next_b, alreadycoll, dVbar_add, domgT_add)
endif

! Single-fiber Jeffery-orbit diagnostic (rigid rod is an exact Jeffery orbit).
if (nfibers == 1) call jeffery_analytics_f(xf(1:lambda,:), step, dble(step)*dt)

end subroutine advance_fibers_rigid


! ─────────────────────────── rigid_collision_pair_f ──────────────────────────
! Hard-core collision impulse between two beads on DIFFERENT rigid rods.
!
! Detection is identical to the flexible fibre_collision_check_f: the closest
! approach of the two spheroid axes (each a segment of half-length
! a_par*(arl-1) about its bead centre) gives the contact normal n_c and the
! contact points r_c0/r_c1 relative to each bead centre; a contact exists when
! that axis distance <= 2*a_par and the contact points are approaching.
!
! The RESPONSE is a rigid-body impulse, not the inextensible-chain solve.  For
! the contact point at R = s_quad(rid)*p + r_c relative to the rod centre of
! mass, with equal rod mass Mfib and transverse moment I_T_fib (the slender rod
! carries no tracked axial spin, so the rotational response is projected onto
! the plane perpendicular to p):
!
!   j = -(1+e) (n·v_rel) / K ,
!   K = 1/M_a + 1/M_b + n·[ (1/I_T)(I-pp)(R_a×n) × R_a ]
!                     + n·[ (1/I_T)(I-pp)(R_b×n) × R_b ] ,
!   ΔVbar = ± j n / M ,   Δomega_T = ± (j/I_T)(I-pp)(R×n) .
!
! Returns count=1 and the four increments when a collision fires, else count=0
! and zero increments.  Purely local (no shared-state writes) → callable from
! the OpenMP cell loop without synchronisation.
subroutine rigid_collision_pair_f(xx0, xx1, ang0, ang1, rid0, rid1, &
     u0, u1, w0, w1, count, dVbar0, dVbar1, domgT0, domgT1)
use param,  only: lx
use fibers, only: a_par, arl, s_quad, Mfib, I_T_fib, coef_res
implicit none
double precision, intent(in)  :: xx0(3), xx1(3), ang0(3), ang1(3)
double precision, intent(in)  :: u0(3), u1(3), w0(3), w1(3)
integer,          intent(in)  :: rid0, rid1
integer,          intent(out) :: count
double precision, intent(out) :: dVbar0(3), dVbar1(3), domgT0(3), domgT1(3)

double precision :: dxx(3), dss(3), r_c0(3), r_c1(3), n_c(3)
double precision :: RA(3), RB(3), vcA(3), vcB(3), cr0(3), cr1(3)
double precision :: rnA(3), rnB(3), projA(3), projB(3), wxrA(3), wxrB(3)
double precision :: s0, s1, q0, q1, t0, t1, t2, dist, relvel, seg
double precision :: jimp, Kmass, invM, invIT
double precision, parameter :: tiny = 1.d-14

count  = 0
dVbar0 = 0.d0; dVbar1 = 0.d0; domgT0 = 0.d0; domgT1 = 0.d0

seg = a_par*(arl-1.d0)      ! half-length of the closest-approach segment

! minimum-image bead-centre-to-bead-centre vector
dxx = xx0 - xx1
if (dxx(1) >  0.5d0*lx) dxx(1) = dxx(1) - lx
if (dxx(1) <= -0.5d0*lx) dxx(1) = dxx(1) + lx
if (dxx(2) >  0.5d0*lx) dxx(2) = dxx(2) - lx
if (dxx(2) <= -0.5d0*lx) dxx(2) = dxx(2) + lx
if (dxx(3) >  0.5d0*lx) dxx(3) = dxx(3) - lx
if (dxx(3) <= -0.5d0*lx) dxx(3) = dxx(3) + lx

! closest approach of the two axes (clamped to the segment half-length)
t0 = sum(dxx*ang0);  t1 = sum(dxx*ang1);  t2 = sum(ang0*ang1)
s0 = (t0 - t1*t2) / (seg*min(t2**2.d0-1.d0, -1.d-12))
s1 = (-t1 + t0*t2) / (seg*min(t2**2.d0-1.d0, -1.d-12))
if (abs(s0) > 1.d0 .and. abs(s1) <= 1.d0) then
   s0 = sign(1.d0,s0)
   s1 = (t1 + seg*s0*t2) / seg
   if (abs(s1) > 1.d0) s1 = sign(1.d0,s1)
elseif (abs(s0) <= 1.d0 .and. abs(s1) > 1.d0) then
   s1 = sign(1.d0,s1)
   s0 = (-t0 + seg*s1*t2) / seg
   if (abs(s0) > 1.d0) s0 = sign(1.d0,s0)
elseif (abs(s0) > 1.d0 .and. abs(s1) > 1.d0) then
   s0 = sign(1.d0,s0);  q1 = (t1 + seg*s0*t2)/seg
   s1 = sign(1.d0,s1);  q0 = (-t0 + seg*s1*t2)/seg
   if (abs(q0) <= 1.d0) then; s0 = q0
   elseif (abs(q1) <= 1.d0) then; s1 = q1
   else; s0 = sign(1.d0,q0); s1 = sign(1.d0,q1); endif
endif

dss  = dxx + s0*seg*ang0 - s1*seg*ang1
dist = sqrt(sum(dss**2.d0))
if (dist > 2.d0*a_par .or. dist <= tiny) return

n_c  = dss / dist
r_c0 = a_par*(s0*(arl-1.d0)*ang0 - n_c)     ! contact pt rel bead0 centre
r_c1 = a_par*(s1*(arl-1.d0)*ang1 + n_c)     ! contact pt rel bead1 centre

! contact-point velocities from the (rigid) bead state:
!   v_c = u_bead + omega_T × r_c = Vbar + omega_T × (s_quad*p + r_c)
cr0(1) = w0(2)*r_c0(3) - w0(3)*r_c0(2)
cr0(2) = w0(3)*r_c0(1) - w0(1)*r_c0(3)
cr0(3) = w0(1)*r_c0(2) - w0(2)*r_c0(1)
cr1(1) = w1(2)*r_c1(3) - w1(3)*r_c1(2)
cr1(2) = w1(3)*r_c1(1) - w1(1)*r_c1(3)
cr1(3) = w1(1)*r_c1(2) - w1(2)*r_c1(1)
vcA = u0 + cr0
vcB = u1 + cr1
relvel = n_c(1)*(vcA(1)-vcB(1)) + n_c(2)*(vcA(2)-vcB(2)) + n_c(3)*(vcA(3)-vcB(3))
if (relvel >= 0.d0) return                  ! separating / grazing — no impulse

! contact point relative to each rod's centre of mass
RA = s_quad(rid0)*ang0 + r_c0
RB = s_quad(rid1)*ang1 + r_c1

invM  = 1.d0/Mfib
invIT = 1.d0/I_T_fib

! transverse rotational inverse-mass terms  n·[ (1/I_T)(I-pp)(R×n) × R ]
rnA(1) = RA(2)*n_c(3) - RA(3)*n_c(2)
rnA(2) = RA(3)*n_c(1) - RA(1)*n_c(3)
rnA(3) = RA(1)*n_c(2) - RA(2)*n_c(1)
projA  = rnA - (rnA(1)*ang0(1)+rnA(2)*ang0(2)+rnA(3)*ang0(3))*ang0
wxrA(1) = projA(2)*RA(3) - projA(3)*RA(2)
wxrA(2) = projA(3)*RA(1) - projA(1)*RA(3)
wxrA(3) = projA(1)*RA(2) - projA(2)*RA(1)

rnB(1) = RB(2)*n_c(3) - RB(3)*n_c(2)
rnB(2) = RB(3)*n_c(1) - RB(1)*n_c(3)
rnB(3) = RB(1)*n_c(2) - RB(2)*n_c(1)
projB  = rnB - (rnB(1)*ang1(1)+rnB(2)*ang1(2)+rnB(3)*ang1(3))*ang1
wxrB(1) = projB(2)*RB(3) - projB(3)*RB(2)
wxrB(2) = projB(3)*RB(1) - projB(1)*RB(3)
wxrB(3) = projB(1)*RB(2) - projB(2)*RB(1)

Kmass = 2.d0*invM &
      + invIT*(n_c(1)*wxrA(1)+n_c(2)*wxrA(2)+n_c(3)*wxrA(3)) &
      + invIT*(n_c(1)*wxrB(1)+n_c(2)*wxrB(2)+n_c(3)*wxrB(3))
if (Kmass <= tiny) return

jimp = -(1.d0+coef_res)*relvel / Kmass      ! >0 since relvel<0

dVbar0 =  (jimp*invM)*n_c
dVbar1 = -(jimp*invM)*n_c
domgT0 =  (jimp*invIT)*projA                ! (1/I_T)(I-pp)(R_a×J),  J=jimp n_c
domgT1 = -(jimp*invIT)*projB
count  = 1
end subroutine rigid_collision_pair_f


! ============================================================================
!  get_fiber_fluid_quantities and advance_fibers are the public API.
!  All subroutines below are internal physics helpers.
! ============================================================================


! ─────────────────────────── tensor_KQRS_f ───────────────────────────────────
! Compute K, Q, R, S for all lambda beads of fiber idfiber.
! Side-effect: fills ORIENTm_0(:,:,:) with O(orientation) via calculate_RS_f
!              → calculate_Jdot_f → calculate_O_f (pass-by-reference chain).
subroutine tensor_KQRS_f(K0, Q0, R0, S0, oo, idfiber, step)
use fibers
use fiber_droplet_mod, only: interface_attraction_target
implicit none
double precision, intent(inout) :: K0(3,lambda), Q0(3,3,lambda)
double precision, intent(inout) :: R0(3,lambda), S0(3,3,lambda), oo(3,3,lambda)
integer, intent(in) :: idfiber, step
integer :: ind, mm
double precision :: zero3(3), u_tgt(3)

K0 = 0.d0; Q0 = 0.d0; R0 = 0.d0; S0 = 0.d0; oo = 0.d0
zero3 = 0.d0

do ind = 1, lambda
   mm = 1+(idfiber-1)*lambda+(ind-1)

   ! Drag target = local fluid velocity, plus the interface-holding attraction
   ! for attached beads.  Folding it into the drag here (rather than overwriting
   ! uf after the solve) lets the Gears constraint re-project it, so it holds the
   ! fiber on the drop without ever breaking inextensibility.
   u_tgt = uf_at_p(mm,:)
   call interface_attraction_target(mm, u_tgt)
   !fiber droplet routine and relatives not touched NBBBBBBBBBBBBBBBBBBBBBBBBBB


   if (fiber_types_flag == 1 .or. bead_types_flag == 1)then
         call calculate_KQ_f_type(K0(:,ind), Q0(:,:,ind), &
         uf(mm,:), u_tgt, &
         zero3, ff(mm,:), &
         angles_f(mm,:),idfiber,ind)

      call calculate_RS_f_type(R0(:,ind), S0(:,:,ind), oo(:,:,ind), &
         omgf(1+(idfiber-1)*lambda+(ind-1),:), &
         vort_at_p(1+(idfiber-1)*lambda+(ind-1),:), &
         grad_at_p(1+(idfiber-1)*lambda+(ind-1),:), &
         angles_f(1+(idfiber-1)*lambda+(ind-1),:), ind, idfiber)     
   
   else
      call calculate_KQ_f(K0(:,ind), Q0(:,:,ind), &
         uf(mm,:), u_tgt, &
         zero3, ff(mm,:), &
         angles_f(mm,:))

      call calculate_RS_f(R0(:,ind), S0(:,:,ind), oo(:,:,ind), &
         omgf(1+(idfiber-1)*lambda+(ind-1),:), &
         vort_at_p(1+(idfiber-1)*lambda+(ind-1),:), &
         grad_at_p(1+(idfiber-1)*lambda+(ind-1),:), &
         angles_f(1+(idfiber-1)*lambda+(ind-1),:), ind, idfiber)
   endif
enddo
end subroutine tensor_KQRS_f


! ─────────────────────────── constraint_flexfibre_f ──────────────────────────
! Thomas-algorithm block-tridiagonal solve for the constraint impulse X.
subroutine constraint_flexfibre_f(XX, KK, QQ, RR, SS, oo, idfiber, &
                                 ft_add,fttrq_add)
use fibers
implicit none
double precision, intent(inout) :: XX(3,lambda)
double precision, intent(in) :: KK(3,lambda), QQ(3,3,lambda)
double precision, intent(in) :: RR(3,lambda), SS(3,3,lambda), oo(3,3,lambda)
integer, intent(in) :: idfiber

double precision :: Ktot(3,lambda-1), NN(3,3,lambda), MM(3,3,lambda-1)
double precision :: G0(3,3,lambda-1), invG0(3,3,lambda-1), F0(3,lambda-1)
double precision :: temp(3,3), res(3)
integer :: ind
! for friction 
double precision, dimension(3,lambda),intent(in) :: ft_add, fttrq_add
double precision,dimension(3, lambda) :: Ftot,fttrqo, FT, FTTRQ
!multifiber
integer :: f


Ktot = 0.d0; NN = 0.d0
!for multifiber
if (fiber_types_flag == 1) then
f = min((idfiber-1)/type_idx + 1, n_fiber_types)
endif

! Known-term assembly
if (fiber_types_flag == 1) then     ! multifiber

   do ind = 1, lambda-1

      Ktot(1,ind) = KK(1,ind) - KK(1,ind+1) &
         - arl_type(f)*a_par_type(f)*(oo(1,1,ind)*RR(1,ind)+oo(1,2,ind)*RR(2,ind)+oo(1,3,ind)*RR(3,ind) &
                  + oo(1,1,ind+1)*RR(1,ind+1)+oo(1,2,ind+1)*RR(2,ind+1)+oo(1,3,ind+1)*RR(3,ind+1))
      Ktot(2,ind) = KK(2,ind) - KK(2,ind+1) &
         - arl_type(f)*a_par_type(f)*(oo(2,1,ind)*RR(1,ind)+oo(2,2,ind)*RR(2,ind)+oo(2,3,ind)*RR(3,ind) &
                  + oo(2,1,ind+1)*RR(1,ind+1)+oo(2,2,ind+1)*RR(2,ind+1)+oo(2,3,ind+1)*RR(3,ind+1))
      Ktot(3,ind) = KK(3,ind) - KK(3,ind+1) &
         - arl_type(f)*a_par_type(f)*(oo(3,1,ind)*RR(1,ind)+oo(3,2,ind)*RR(2,ind)+oo(3,3,ind)*RR(3,ind) &
                  + oo(3,1,ind+1)*RR(1,ind+1)+oo(3,2,ind+1)*RR(2,ind+1)+oo(3,3,ind+1)*RR(3,ind+1))

   enddo

elseif (bead_types_flag == 1) then     ! asymmetrical beads

   do ind = 1, lambda-1
      if (ind == lambda -1) then       ! calculating just last asymbead
         Ktot(1,ind) = KK(1,ind) - KK(1,ind+1) &
            - arl*a_par*(oo(1,1,ind)*RR(1,ind)+oo(1,2,ind)*RR(2,ind)+oo(1,3,ind)*RR(3,ind)) &
            - arl_lb*a_par_lb* (oo(1,1,ind+1)*RR(1,ind+1)+oo(1,2,ind+1)*RR(2,ind+1)+oo(1,3,ind+1)*RR(3,ind+1))
         Ktot(2,ind) = KK(2,ind) - KK(2,ind+1) &
            - arl*a_par* (oo(2,1,ind)*RR(1,ind)+oo(2,2,ind)*RR(2,ind)+oo(2,3,ind)*RR(3,ind)) &
            - arl_lb*a_par_lb* (oo(2,1,ind+1)*RR(1,ind+1)+oo(2,2,ind+1)*RR(2,ind+1)+oo(2,3,ind+1)*RR(3,ind+1))
         Ktot(3,ind) = KK(3,ind) - KK(3,ind+1) &
            - arl*a_par* (oo(3,1,ind)*RR(1,ind)+oo(3,2,ind)*RR(2,ind)+oo(3,3,ind)*RR(3,ind)) &
            - arl_lb*a_par_lb* (oo(3,1,ind+1)*RR(1,ind+1)+oo(3,2,ind+1)*RR(2,ind+1)+oo(3,3,ind+1)*RR(3,ind+1))
      else     ! normal calc for the other beads
         Ktot(1,ind) = KK(1,ind) - KK(1,ind+1) &
         - arl*a_par*(oo(1,1,ind)*RR(1,ind)+oo(1,2,ind)*RR(2,ind)+oo(1,3,ind)*RR(3,ind) &
                  + oo(1,1,ind+1)*RR(1,ind+1)+oo(1,2,ind+1)*RR(2,ind+1)+oo(1,3,ind+1)*RR(3,ind+1))
         Ktot(2,ind) = KK(2,ind) - KK(2,ind+1) &
         - arl*a_par*(oo(2,1,ind)*RR(1,ind)+oo(2,2,ind)*RR(2,ind)+oo(2,3,ind)*RR(3,ind) &
                  + oo(2,1,ind+1)*RR(1,ind+1)+oo(2,2,ind+1)*RR(2,ind+1)+oo(2,3,ind+1)*RR(3,ind+1))
         Ktot(3,ind) = KK(3,ind) - KK(3,ind+1) &
         - arl*a_par*(oo(3,1,ind)*RR(1,ind)+oo(3,2,ind)*RR(2,ind)+oo(3,3,ind)*RR(3,ind) &
                  + oo(3,1,ind+1)*RR(1,ind+1)+oo(3,2,ind+1)*RR(2,ind+1)+oo(3,3,ind+1)*RR(3,ind+1))
      endif 
   enddo

else
   do ind = 1, lambda-1
      Ktot(1,ind) = KK(1,ind) - KK(1,ind+1) &
         - arl*a_par*(oo(1,1,ind)*RR(1,ind)+oo(1,2,ind)*RR(2,ind)+oo(1,3,ind)*RR(3,ind) &
                  + oo(1,1,ind+1)*RR(1,ind+1)+oo(1,2,ind+1)*RR(2,ind+1)+oo(1,3,ind+1)*RR(3,ind+1))
      Ktot(2,ind) = KK(2,ind) - KK(2,ind+1) &
         - arl*a_par*(oo(2,1,ind)*RR(1,ind)+oo(2,2,ind)*RR(2,ind)+oo(2,3,ind)*RR(3,ind) &
                  + oo(2,1,ind+1)*RR(1,ind+1)+oo(2,2,ind+1)*RR(2,ind+1)+oo(2,3,ind+1)*RR(3,ind+1))
      Ktot(3,ind) = KK(3,ind) - KK(3,ind+1) &
         - arl*a_par*(oo(3,1,ind)*RR(1,ind)+oo(3,2,ind)*RR(2,ind)+oo(3,3,ind)*RR(3,ind) &
                  + oo(3,1,ind+1)*RR(1,ind+1)+oo(3,2,ind+1)*RR(2,ind+1)+oo(3,3,ind+1)*RR(3,ind+1))
   enddo
endif

! N = oo*S*oo * (arl*a_par)
!multifiber
if (fiber_types_flag == 1) then
   do ind = 1, lambda 
      call calculate_prod_array_f(SS(:,:,ind), oo(:,:,ind), temp)
      call calculate_prod_array_f(oo(:,:,ind), temp, NN(:,:,ind))
      NN(:,:,ind) = NN(:,:,ind) * arl_type(f)*a_par_type(f)
   enddo
!multibead
elseif (bead_types_flag == 1) then 
   do ind = 1, lambda
      if (ind == lambda) then
         call calculate_prod_array_f(SS(:,:,ind), oo(:,:,ind), temp)
         call calculate_prod_array_f(oo(:,:,ind), temp, NN(:,:,ind))
         NN(:,:,ind) = NN(:,:,ind) * arl_lb*a_par_lb
      else
         call calculate_prod_array_f(SS(:,:,ind), oo(:,:,ind), temp)
         call calculate_prod_array_f(oo(:,:,ind), temp, NN(:,:,ind))
         NN(:,:,ind) = NN(:,:,ind) * arl*a_par
      endif
   enddo
!normal
else
   do ind = 1, lambda
      call calculate_prod_array_f(SS(:,:,ind), oo(:,:,ind), temp)
      call calculate_prod_array_f(oo(:,:,ind), temp, NN(:,:,ind))
      NN(:,:,ind) = NN(:,:,ind) * arl*a_par
   enddo
endif

! M = -Q_i - Q_{i+1} + N_i + N_{i+1}
do ind = 1, lambda-1
   MM(:,:,ind) = -QQ(:,:,ind) - QQ(:,:,ind+1) + NN(:,:,ind) + NN(:,:,ind+1)
enddo

! Update N <- N + Q
do ind = 1, lambda
   NN(:,:,ind) = NN(:,:,ind) + QQ(:,:,ind)
enddo

! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! calculate the friction term forcing

do ind=1,lambda-1
  ! same ind
  !-------------  
  ! product of forces with Q matrix 
  call vec_product_f(QQ(:,:,ind),ft_add(:,ind),FT(:,ind))
  ! product of torques with S matrix 
  call vec_product_f(SS(:,:,ind), fttrq_add(:,ind),fttrqo(:,ind))
  ! now we want to express it also in function of the a and lambda O
  call vec_product_f(oo(:,:,ind),fttrqo(:,ind), FTTRQ(:,ind) )

  !next ind 
  call vec_product_f(QQ(:,:,ind+1),ft_add(:,ind+1),FT(:,ind+1))
  call vec_product_f(SS(:,:,ind+1), fttrq_add(:,ind+1),fttrqo(:,ind+1))
  call vec_product_f(oo(:,:,ind+1),fttrqo(:,ind+1), FTTRQ(:,ind+1))


  !multiplying by the a and lambda
  if (fiber_types_flag == 1) then
  FTTRQ(:,ind) = FTTRQ(:,ind) * arl_type(f) * a_par_type(f)
  FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl_type(f) * a_par_type(f)

  elseif (bead_types_flag == 1) then
   if (ind == lambda -1) then
      FTTRQ(:,ind) = FTTRQ(:,ind) * arl * a_par
      FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl_lb * a_par_lb
   else
      FTTRQ(:,ind) = FTTRQ(:,ind) * arl * a_par
      FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl * a_par
   endif

  else
  FTTRQ(:,ind) = FTTRQ(:,ind) * arl * a_par
  FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl * a_par
  endif


  ! summing those thoghether 
  Ftot(1,ind) = -FT(1,ind+1) + FT(1,ind) - FTTRQ(1,ind+1) - FTTRQ(1,ind)
  Ftot(2,ind) = -FT(2,ind+1) + FT(2,ind) - FTTRQ(2,ind+1) - FTTRQ(2,ind)
  Ftot(3,ind) = -FT(3,ind+1) + FT(3,ind) - FTTRQ(3,ind+1) - FTTRQ(3,ind)

  ! adding to k_tot

  Ktot(1,ind) = Ktot(1,ind) + Ftot(1,ind)
  Ktot(2,ind) = Ktot(2,ind) + Ftot(2,ind) 
  Ktot(3,ind) = Ktot(3,ind) + Ftot(3,ind)

  ! to check the callable subroutine eventually
  !call friction_calculation(QQ,SS,RR,oo,fttrq_add,Ktot,idfiber)

enddo

! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

! Backward sweep (Thomas)
F0 = 0.d0; invG0 = 0.d0; G0 = 0.d0
do ind = lambda-1, 1, -1
   if (ind == lambda-1) then
      F0(:,ind)    = Ktot(:,ind)
      G0(:,:,ind)  = MM(:,:,ind)
      call invertmatrix_constraint_f(G0(:,:,ind), invG0(:,:,ind))
   else
      call vec_product_f(invG0(:,:,ind+1), F0(:,ind+1), res)
      call vec_product_f(NN(:,:,ind+1), res, F0(:,ind))
      F0(:,ind) = Ktot(:,ind) - F0(:,ind)
      call calculate_prod_array_f(invG0(:,:,ind+1), NN(:,:,ind+1), temp)
      call calculate_prod_array_f(NN(:,:,ind+1), temp, G0(:,:,ind))
      G0(:,:,ind) = MM(:,:,ind) - G0(:,:,ind)
      call invertmatrix_constraint_f(G0(:,:,ind), invG0(:,:,ind))
   endif
enddo

! Forward substitution
XX = 0.d0
do ind = 1, lambda-1
   if (ind == 1) then
      call vec_product_f(invG0(:,:,ind), F0(:,ind), XX(:,ind+1))
   else
      call vec_product_f(NN(:,:,ind), XX(:,ind), res)
      F0(:,ind) = F0(:,ind) - res
      call vec_product_f(invG0(:,:,ind), F0(:,ind), XX(:,ind+1))
   endif
enddo

! Store constraint impulse for output (cols 32-34 of fiber_*.dat)
! X(:,1) is always zero (free end); X(:,2..lambda) are the inter-bead impulses.
do ind = 1, lambda
   X_constraint(1+(idfiber-1)*lambda+(ind-1), :) = XX(:,ind)
enddo
end subroutine constraint_flexfibre_f


! ─────────────────────────── dynamics_flexible_fiber_f ───────────────────────
! Compute bead velocities from constraint impulses XX and K,Q,R,S tensors.
subroutine dynamics_flexible_fiber_f(XX, KK, QQ, RR, SS, oo, idfiber,ft_add,fttrq_add)
use fibers
implicit none
double precision, intent(in) :: XX(3,lambda), KK(3,lambda)
double precision, intent(in) :: QQ(3,3,lambda), RR(3,lambda)
double precision, intent(in) :: SS(3,3,lambda), oo(3,3,lambda)
integer, intent(in) :: idfiber
integer :: ind
double precision :: velp(3), omegap(3), temp(3)
! for friction 
double precision, dimension(3,lambda),intent(in) :: ft_add, fttrq_add
double precision, dimension(3,lambda) :: FT,TT,trqt

velp = 0.d0; omegap = 0.d0

do ind = 1, lambda

   !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
   ! calculation of the friction term - translation 
   call vec_product_f(QQ(:,:,ind),ft_add(:,ind),FT(:,ind))
   ! rotation term
   call vec_product_f(oo(:,:,ind),fttrq_add(:,ind),trqt(:,ind))
   call vec_product_f(SS(:,:,ind),trqt(:,ind),TT(:,ind))
   !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

   if (ind == 1) then
      call vec_product_f(QQ(:,:,ind), XX(:,ind+1), velp)
   elseif (ind < lambda) then
      call vec_product_f(QQ(:,:,ind), XX(:,ind+1)-XX(:,ind), velp)
   else
      call vec_product_f(QQ(:,:,ind), -XX(:,ind), velp)
   endif
   !velp = velp + KK(:,ind)
   ! updated version
   velp = velp + KK(:,ind) + FT(:,ind)

   if (ind == 1) then
      call vec_product_f(oo(:,:,ind), XX(:,ind+1), temp)
      call vec_product_f(SS(:,:,ind), temp, omegap)
   elseif (ind < lambda) then
      call vec_product_f(oo(:,:,ind), XX(:,ind)+XX(:,ind+1), temp)
      call vec_product_f(SS(:,:,ind), temp, omegap)
   else
      call vec_product_f(oo(:,:,ind), XX(:,ind), temp)
      call vec_product_f(SS(:,:,ind), temp, omegap)
   endif
   !omegap = omegap + RR(:,ind)
   ! correction with friction
   omegap = omegap + RR(:,ind) + TT(:,ind)

   uf(1+(idfiber-1)*lambda+(ind-1),:)  = velp
   omgf(1+(idfiber-1)*lambda+(ind-1),:) = omegap
enddo
end subroutine dynamics_flexible_fiber_f


! ─────────────────────────── integrate_flexible_fiber_f ──────────────────────
! Adams-Bashforth 3rd-order position advance + Euler orientation update.
! Fully periodic in all 3 directions (HIT: no walls).
subroutine integrate_flexible_fiber_f(i, step, substep)
use param
use fibers
implicit none
integer, intent(in) :: i, step, substep
integer :: j, totalstep

if (subiterations > 1) then
   totalstep = substep - 1 + (step - nstart) + nstart
else
   totalstep = step
endif

do j = 0, lambda-1
   associate(xp => xf(1+(i-1)*lambda+j,:), &
             vv  => uf(1+(i-1)*lambda+j,:), &
             v1  => vp_o(1+(i-1)*lambda+j,:,1), &
             v2  => vp_o(1+(i-1)*lambda+j,:,2), &
             ang => angles_f(1+(i-1)*lambda+j,:), &
             omg => omgf(1+(i-1)*lambda+j,:))

      if (totalstep == nstart) then
         xp = xp + dt_part * vv
         v1 = vv
      elseif (totalstep == nstart+1) then
         xp = xp + dt_part * (1.5d0*vv - 0.5d0*v1)
         v2 = v1;  v1 = vv
      else
         xp = xp + dt_part * (23.d0/12.d0*vv - 16.d0/12.d0*v1 + 5.d0/12.d0*v2)
         v2 = v1;  v1 = vv
      endif

      call rotate_fiber_f(ang, omg, totalstep, 1+(i-1)*lambda+j)

      ! Periodic BC – all three directions
      if (xp(1) < 0.d0) xp(1) = xp(1) + lx
      if (xp(1) > lx)   xp(1) = xp(1) - lx
      if (xp(2) < 0.d0) xp(2) = xp(2) + lx
      if (xp(2) > lx)   xp(2) = xp(2) - lx
      if (xp(3) < 0.d0) xp(3) = xp(3) + lx
      if (xp(3) > lx)   xp(3) = xp(3) - lx
   end associate
enddo
end subroutine integrate_flexible_fiber_f


! ─────────────────────────── rotate_fiber_f ──────────────────────────────────
! Adams-Bashforth 3rd-order orientation update (no normalisation drift).
subroutine rotate_fiber_f(angls, omega, totalstep, ind)
use param
use fibers
implicit none
double precision, intent(inout) :: angls(3)
double precision, intent(in)    :: omega(3)
integer, intent(in) :: totalstep, ind
double precision :: orient(3), orient0(3), fac

orient0 = angls
orient(1) =  orient0(2)*omega(3) - omega(2)*orient0(3)
orient(2) = -(orient0(1)*omega(3) - omega(1)*orient0(3))
orient(3) =  orient0(1)*omega(2) - omega(1)*orient0(2)

if (totalstep == nstart) then
   angls         = orient0 - dt_part * orient
   angl_o(ind,:,1) = orient
elseif (totalstep == nstart+1) then
   angls         = orient0 - dt_part * (1.5d0*orient - 0.5d0*angl_o(ind,:,1))
   angl_o(ind,:,2) = angl_o(ind,:,1)
   angl_o(ind,:,1) = orient
else
   angls         = orient0 - dt_part * (23.d0/12.d0*orient &
                   - 16.d0/12.d0*angl_o(ind,:,1) + 5.d0/12.d0*angl_o(ind,:,2))
   angl_o(ind,:,2) = angl_o(ind,:,1)
   angl_o(ind,:,1) = orient
endif

! Renormalise
fac = sqrt(angls(1)**2.d0 + angls(2)**2.d0 + angls(3)**2.d0)
if (fac > 1.d-14) angls = angls / fac
end subroutine rotate_fiber_f


! ─────────────────────────── quality_check_f ─────────────────────────────────
! Check inextensibility error and correct bead positions from the central bead.
! All three directions are periodic (HIT).
subroutine quality_check_f(xxx, angl, upart, omgpart, idfiber)
use param
use fibers
implicit none
double precision, intent(inout) :: xxx(lambda,3), angl(lambda,3)
double precision, intent(in)    :: upart(lambda,3), omgpart(lambda,3)
integer, intent(in) :: idfiber
integer :: j, locks
double precision :: sum_dist, dx12, dy12, dz12, per_seg
integer :: f 

if (fiber_types_flag == 1) then
f = min((idfiber-1)/type_idx + 1, n_fiber_types)
endif

locks = int(lambda/2.0d0) + mod(lambda,2)
sum_dist = 0.d0


if (fiber_types_flag == 1) then
   do j = 1, lambda-1
      dx12 = xxx(j+1,1)-xxx(j,1)
      if (dx12 >  0.5d0*lx) dx12 = dx12 - lx
      if (dx12 <= -0.5d0*lx) dx12 = dx12 + lx
      dy12 = xxx(j+1,2)-xxx(j,2)
      if (dy12 >  0.5d0*lx) dy12 = dy12 - lx
      if (dy12 <= -0.5d0*lx) dy12 = dy12 + lx
      dz12 = xxx(j+1,3)-xxx(j,3)
      if (dz12 >  0.5d0*lx) dz12 = dz12 - lx
      if (dz12 <= -0.5d0*lx) dz12 = dz12 + lx
      sum_dist = sum_dist + sqrt( &
         (arl_type(f)*a_par_type(f)*(angl(j,1)+angl(j+1,1))-dx12)**2.d0 + &
         (arl_type(f)*a_par_type(f)*(angl(j,2)+angl(j+1,2))-dy12)**2.d0 + &
         (arl_type(f)*a_par_type(f)*(angl(j,3)+angl(j+1,3))-dz12)**2.d0) / &
         (2.d0*a_par_type(f)*arl_type(f))
   enddo

   per_seg = sum_dist / dble(lambda - 1)

   ! Per-segment violation (lambda-independent threshold).
   ! The original design targeted lambda=20 (threshold 0.01/19 ≈ 5e-4 per segment).
   ! A fixed absolute threshold of 0.01 becomes 5x tighter for lambda=100.

   per_seg = sum_dist / dble(lambda - 1)

   if (per_seg >= 0.1d0) then
      ! Truly diverged (>10% per segment): abort.
      write(*,'(a,es10.3,a,i6)') 'Tolerance broken! err=', sum_dist, ' fiber=', idfiber
      call exit(1)
   elseif (per_seg >= 5.d-4) then
      ! Mild violation: let the position rebuild below correct it; warn once.
      write(*,'(a,es10.3,a,es10.3,a,i6)') &
         'WARNING constraint err=', sum_dist, ' (per-seg=', per_seg, ') fiber=', idfiber
   endif

   ! Rebuild from centre outward with per-step wrapping.
   ! Without wrapping here, the cumulative displacement from the centre bead
   ! can reach lambda/2 * segment_length (e.g. 50 * 0.5 = 25 [L] for lambda=100),
   ! far exceeding the two levels of single-step wrap that follow.  This puts beads
   ! at wrong absolute positions (> lx), which overflows the collision cell-list index.
   ! Wrapping after each bead is equivalent to wrapping all-at-once for the case
   ! where one segment << lx (always true here: segment=0.5 << lx≈6.28).
   do j = locks+1, lambda
      xxx(j,:) = xxx(j-1,:) + arl_type(f)*a_par_type(f)*(angl(j-1,:)+angl(j,:))
      if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
      if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
      if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
      if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
      if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
      if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
   enddo
   do j = locks-1, 1, -1
      xxx(j,:) = xxx(j+1,:) - arl_type(f)*a_par_type(f)*(angl(j+1,:)+angl(j,:))
      if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
      if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
      if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
      if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
      if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
      if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
   enddo

elseif (bead_types_flag == 1) then
   do j = 1, lambda-1
      if (j == lambda -1) then
         dx12 = xxx(j+1,1)-xxx(j,1)
         if (dx12 >  0.5d0*lx) dx12 = dx12 - lx
         if (dx12 <= -0.5d0*lx) dx12 = dx12 + lx
         dy12 = xxx(j+1,2)-xxx(j,2)
         if (dy12 >  0.5d0*lx) dy12 = dy12 - lx
         if (dy12 <= -0.5d0*lx) dy12 = dy12 + lx
         dz12 = xxx(j+1,3)-xxx(j,3)
         if (dz12 >  0.5d0*lx) dz12 = dz12 - lx
         if (dz12 <= -0.5d0*lx) dz12 = dz12 + lx
         sum_dist = sum_dist + sqrt( &
            ((arl*a_par*angl(j,1)+ arl_lb*a_par_lb*angl(j+1,1))-dx12)**2.d0 + &
            ((arl*a_par*angl(j,2)+ arl_lb*a_par_lb*angl(j+1,2))-dy12)**2.d0 + &
            ((arl*a_par*angl(j,3)+ arl_lb*a_par_lb*angl(j+1,3))-dz12)**2.d0) / &
            (arl*a_par + arl_lb*a_par_lb)     

      else 

         dx12 = xxx(j+1,1)-xxx(j,1)
         if (dx12 >  0.5d0*lx) dx12 = dx12 - lx
         if (dx12 <= -0.5d0*lx) dx12 = dx12 + lx
         dy12 = xxx(j+1,2)-xxx(j,2)
         if (dy12 >  0.5d0*lx) dy12 = dy12 - lx
         if (dy12 <= -0.5d0*lx) dy12 = dy12 + lx
         dz12 = xxx(j+1,3)-xxx(j,3)
         if (dz12 >  0.5d0*lx) dz12 = dz12 - lx
         if (dz12 <= -0.5d0*lx) dz12 = dz12 + lx
         sum_dist = sum_dist + sqrt( &
            (arl*a_par*(angl(j,1)+angl(j+1,1))-dx12)**2.d0 + &
            (arl*a_par*(angl(j,2)+angl(j+1,2))-dy12)**2.d0 + &
            (arl*a_par*(angl(j,3)+angl(j+1,3))-dz12)**2.d0) / (2.d0*a_par*arl)
      endif
   enddo

   per_seg = sum_dist / dble(lambda - 1)

   if (per_seg >= 0.1d0) then
      ! Truly diverged (>10% per segment): abort.
      write(*,'(a,es10.3,a,i6)') 'Tolerance broken! err=', sum_dist, ' fiber=', idfiber
      call exit(1)
   elseif (per_seg >= 5.d-4) then
      ! Mild violation: let the position rebuild below correct it; warn once.
      write(*,'(a,es10.3,a,es10.3,a,i6)') &
         'WARNING constraint err=', sum_dist, ' (per-seg=', per_seg, ') fiber=', idfiber
   endif

   do j = locks+1, lambda
      if ( j == lambda) then 
         xxx(j,:) = xxx(j-1,:) + arl*a_par*angl(j-1,:)+arl_lb*a_par_lb*angl(j,:)
         if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
         if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
         if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
         if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
         if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
         if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
      else
         xxx(j,:) = xxx(j-1,:) + arl*a_par*(angl(j-1,:)+angl(j,:))
         if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
         if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
         if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
         if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
         if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
         if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
      endif
   enddo
   do j = locks-1, 1, -1
      xxx(j,:) = xxx(j+1,:) - arl*a_par*(angl(j+1,:)+angl(j,:))
      if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
      if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
      if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
      if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
      if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
      if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
   enddo

else     ! normal case

   do j = 1, lambda-1
      dx12 = xxx(j+1,1)-xxx(j,1)
      if (dx12 >  0.5d0*lx) dx12 = dx12 - lx
      if (dx12 <= -0.5d0*lx) dx12 = dx12 + lx
      dy12 = xxx(j+1,2)-xxx(j,2)
      if (dy12 >  0.5d0*lx) dy12 = dy12 - lx
      if (dy12 <= -0.5d0*lx) dy12 = dy12 + lx
      dz12 = xxx(j+1,3)-xxx(j,3)
      if (dz12 >  0.5d0*lx) dz12 = dz12 - lx
      if (dz12 <= -0.5d0*lx) dz12 = dz12 + lx
      sum_dist = sum_dist + sqrt( &
         (arl*a_par*(angl(j,1)+angl(j+1,1))-dx12)**2.d0 + &
         (arl*a_par*(angl(j,2)+angl(j+1,2))-dy12)**2.d0 + &
         (arl*a_par*(angl(j,3)+angl(j+1,3))-dz12)**2.d0) / (2.d0*a_par*arl)
   enddo

   per_seg = sum_dist / dble(lambda - 1)

   if (per_seg >= 0.1d0) then
      ! Truly diverged (>10% per segment): abort.
      write(*,'(a,es10.3,a,i6)') 'Tolerance broken! err=', sum_dist, ' fiber=', idfiber
      call exit(1)
   elseif (per_seg >= 5.d-4) then
      ! Mild violation: let the position rebuild below correct it; warn once.
      write(*,'(a,es10.3,a,es10.3,a,i6)') &
         'WARNING constraint err=', sum_dist, ' (per-seg=', per_seg, ') fiber=', idfiber
   endif

   do j = locks+1, lambda
      xxx(j,:) = xxx(j-1,:) + arl*a_par*(angl(j-1,:)+angl(j,:))
      if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
      if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
      if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
      if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
      if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
      if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
   enddo
   do j = locks-1, 1, -1
      xxx(j,:) = xxx(j+1,:) - arl*a_par*(angl(j+1,:)+angl(j,:))
      if (xxx(j,1) < 0.d0) xxx(j,1) = xxx(j,1)+lx
      if (xxx(j,1) > lx)   xxx(j,1) = xxx(j,1)-lx
      if (xxx(j,2) < 0.d0) xxx(j,2) = xxx(j,2)+lx
      if (xxx(j,2) > lx)   xxx(j,2) = xxx(j,2)-lx
      if (xxx(j,3) < 0.d0) xxx(j,3) = xxx(j,3)+lx
      if (xxx(j,3) > lx)   xxx(j,3) = xxx(j,3)-lx
   enddo
endif
end subroutine quality_check_f


! ─────────────────────────── calculate_KQ_f ──────────────────────────────────
! K = unconstrained translational velocity; Q = translational compliance tensor.
! Semi-implicit Stokes drag (Re_p = 0, one-way coupling).
subroutine calculate_KQ_f(KK, QQ, vv, uu, disturbance, fdrag, ang)
use param
use fibers
implicit none
double precision, intent(out)  :: KK(3), QQ(3,3)
double precision, intent(in)   :: vv(3), uu(3), disturbance(3), fdrag(3), ang(3)
double precision :: AA(3,3), res(3,3), invres(3,3), temp(3)
double precision :: OZZ(3,3), line
!added for gravity diag
double precision :: gmag, dgmag

call calculate_ooT_f(OZZ, ang)

! AA = 9ν/(2 a² dens_part) * [Y_a I + (X_a-Y_a) pp]
! from: dv/dt = F/m, F=6πμa[Y_a I+(X_a-Y_a)pp](u-v), m=(4/3)πa³ρ_p
! factor 9ν/2 = 3*(3ν/2); written as 3*mu with an explicit (3/2) factor
AA(1,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,1)+(X_a-Y_a)*OZZ(1,1))
AA(2,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,1)+(X_a-Y_a)*OZZ(2,1))
AA(3,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,1)+(X_a-Y_a)*OZZ(3,1))
AA(1,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,2)+(X_a-Y_a)*OZZ(1,2))
AA(2,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,2)+(X_a-Y_a)*OZZ(2,2))
AA(3,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,2)+(X_a-Y_a)*OZZ(3,2))
AA(1,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,3)+(X_a-Y_a)*OZZ(1,3))
AA(2,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,3)+(X_a-Y_a)*OZZ(2,3))
AA(3,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,3)+(X_a-Y_a)*OZZ(3,3))

res = Idm + dt_part * AA
call INVERSE_MATRIX_ff(res, invres)

call vec_product_f(AA, uu, temp)
temp = vv + dt_part * temp
! Buoyancy-reduced gravity as a body-force impulse on the bead (backward Euler:
! folded into the RHS before the drag inverse).  geff=0 unless gravity_flag=1, so
! this is a no-op — bit-for-bit identical — when gravity is off.
if (gravity_flag == 1) temp = temp + dt_part * geff

call vec_product_f(invres, temp, KK)

line = dt_part / (2.d0*pi*arl*a_par**3.d0*dens_part)
QQ   = line * invres
end subroutine calculate_KQ_f


! ─────────────────────────── calculate_RS_f ──────────────────────────────────
! R = unconstrained angular velocity; S = rotational compliance tensor.
! Side-effect: fills oo with O(orientation) via calculate_Jdot_f.
subroutine calculate_RS_f(RR, SS, oo, omegap, vortph, pgradh, o, ind, idfiber)
use param
use fibers
implicit none
double precision, intent(out)   :: RR(3), SS(3,3)
double precision, intent(inout) :: oo(3,3)
double precision, intent(in)    :: omegap(3), vortph(3), pgradh(9), o(3)
integer, intent(in) :: ind, idfiber

double precision :: JJ(3,3), JJdot(3,3), CC(3,3), HH(3), E(3,3)
double precision :: res(3,3), invres(3,3), temp(3), temp2(3), torque(3)
double precision :: line

call calculate_J_f(JJ, o)

E(1,1) = pgradh(1);          E(1,2) = 0.5d0*(pgradh(2)+pgradh(4))
E(1,3) = 0.5d0*(pgradh(3)+pgradh(7))
E(2,1) = E(1,2);              E(2,2) = pgradh(5)
E(2,3) = 0.5d0*(pgradh(6)+pgradh(8))
E(3,1) = E(1,3);              E(3,2) = E(2,3);  E(3,3) = pgradh(9)

call calculate_H_f(HH, E, o)
call calculate_C_f(CC, o)
call calculate_Jdot_f(JJdot, oo, omegap, o)   ! also fills oo = O(o)
call calculate_torque_f(torque, o, ind, idfiber)
! Store bending torque for output (cols 35-37 of fiber_*.dat)
torque_bend(1+(idfiber-1)*lambda+(ind-1), :) = torque

res = JJ + dt_part * (JJdot + CC)
call INVERSE_MATRIX_ff(res, invres)

call vec_product_f(JJ, omegap, temp2)
call vec_product_f(CC, vortph, temp)
temp2 = temp2 + dt_part*(temp + HH) + torque
call vec_product_f(invres, temp2, RR)

line = dt_part / (pi * a_par**4.d0 * dens_part)
SS   = line * invres
end subroutine calculate_RS_f


! mmmmmmmmmmmmmmmmmmmm added for multifiber case mmmmmmmmmmmmmmmmmm


! ─────────────────────────── calculate_KQ_f multi──────────────────────────────────
! K = unconstrained translational velocity; Q = translational compliance tensor.
! Semi-implicit Stokes drag (Re_p = 0, one-way coupling).
subroutine calculate_KQ_f_type(KK, QQ, vv, uu, disturbance, fdrag, ang,idfiber,ind)
use param
use fibers
implicit none
double precision, intent(out)  :: KK(3), QQ(3,3)
double precision, intent(in)   :: vv(3), uu(3), disturbance(3), fdrag(3), ang(3)
double precision :: AA(3,3), res(3,3)
double precision :: invres(3,3), temp(3)
double precision :: OZZ(3,3), line
!added for gravity diag
double precision :: gmag, dgmag
!added for multifiber
integer,intent(in) :: idfiber
integer,intent(in) :: ind
integer :: f

call calculate_ooT_f(OZZ, ang)
! to get the actuall right term adopting the right fiber values
!NB fibers are grouped in progressive way => 1st group (1-3), 2nd group (4-6)...

if (fiber_types_flag == 1) then

   f = min((idfiber-1)/type_idx + 1, n_fiber_types)

   ! AA = 9ν/(2 a² dens_part) * [Y_a I + (X_a-Y_a) pp]
   ! from: dv/dt = F/m, F=6πμa[Y_a I+(X_a-Y_a)pp](u-v), m=(4/3)πa³ρ_p
   ! factor 9ν/2 = 3*(3ν/2); written as 3*mu with an explicit (3/2) factor
   AA(1,1) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
      )*(Y_a_type(f)*Idm(1,1)+(X_a_type(f)-Y_a_type(f))*OZZ(1,1))
   AA(2,1) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(2,1)+(X_a_type(f)-Y_a_type(f))*OZZ(2,1))
   AA(3,1) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(3,1)+(X_a_type(f)-Y_a_type(f))*OZZ(3,1))
   AA(1,2) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(1,2)+(X_a_type(f)-Y_a_type(f))*OZZ(1,2))
   AA(2,2) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(2,2)+(X_a_type(f)-Y_a_type(f))*OZZ(2,2))
   AA(3,2) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(3,2)+(X_a_type(f)-Y_a_type(f))*OZZ(3,2))
   AA(1,3) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(1,3)+(X_a_type(f)-Y_a_type(f))*OZZ(1,3))
   AA(2,3) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(2,3)+(X_a_type(f)-Y_a_type(f))*OZZ(2,3))
   AA(3,3) = 9.d0*mu_f/(2.d0*arl_type(f)*a_par_type(f)**2.d0*dens_part_type(f) &
         )*(Y_a_type(f)*Idm(3,3)+(X_a_type(f)-Y_a_type(f))*OZZ(3,3))
      
   res = Idm + dt_part * AA

elseif( bead_types_flag == 1) then   !asymmetry

   if (ind == lambda) then
      AA(1,1) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(1,1)+(X_a_lb-Y_a_lb)*OZZ(1,1))
      AA(2,1) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(2,1)+(X_a_lb-Y_a_lb)*OZZ(2,1))
      AA(3,1) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(3,1)+(X_a_lb-Y_a_lb)*OZZ(3,1))
      AA(1,2) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(1,2)+(X_a_lb-Y_a_lb)*OZZ(1,2))
      AA(2,2) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(2,2)+(X_a_lb-Y_a_lb)*OZZ(2,2))
      AA(3,2) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(3,2)+(X_a_lb-Y_a_lb)*OZZ(3,2))
      AA(1,3) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(1,3)+(X_a_lb-Y_a_lb)*OZZ(1,3))
      AA(2,3) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(2,3)+(X_a_lb-Y_a_lb)*OZZ(2,3))
      AA(3,3) = 9.d0*mu_f/(2.d0*arl_lb*a_par_lb**2.d0*dens_part)* &
         (Y_a_lb*Idm(3,3)+(X_a_lb-Y_a_lb)*OZZ(3,3))

   else 

      AA(1,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,1)+(X_a-Y_a)*OZZ(1,1))
      AA(2,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,1)+(X_a-Y_a)*OZZ(2,1))
      AA(3,1) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,1)+(X_a-Y_a)*OZZ(3,1))
      AA(1,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,2)+(X_a-Y_a)*OZZ(1,2))
      AA(2,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,2)+(X_a-Y_a)*OZZ(2,2))
      AA(3,2) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,2)+(X_a-Y_a)*OZZ(3,2))
      AA(1,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(1,3)+(X_a-Y_a)*OZZ(1,3))
      AA(2,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(2,3)+(X_a-Y_a)*OZZ(2,3))
      AA(3,3) = 9.d0*mu_f/(2.d0*arl*a_par**2.d0*dens_part)*(Y_a*Idm(3,3)+(X_a-Y_a)*OZZ(3,3))
   
   endif

   res = Idm + dt_part * AA


endif


call INVERSE_MATRIX_ff(res, invres)

call vec_product_f(AA, uu, temp)
temp = vv + dt_part * temp
! Buoyancy-reduced gravity as a body-force impulse on the bead (backward Euler:
! folded into the RHS before the drag inverse).  geff=0 unless gravity_flag=1, so
! this is a no-op — bit-for-bit identical — when gravity is off.
if (fiber_types_flag == 1) then

   if (gravity_flag == 1) temp = temp + dt_part * geff_type(f,:)

   call vec_product_f(invres, temp, KK)

   line = dt_part / (2.d0*pi*arl_type(f)*a_par_type(f)**3.d0*dens_part_type(f))

elseif (bead_types_flag == 1) then

   if (gravity_flag == 1) temp = temp + dt_part * geff(:)
   call vec_product_f(invres, temp, KK)

   if (ind == lambda) then 
      line = dt_part / (2.d0*pi*arl_lb*a_par_lb**3.d0*dens_part)
   else
      line = dt_part / (2.d0*pi*arl*a_par**3.d0*dens_part) 
   endif
endif

QQ   = line * invres


end subroutine calculate_KQ_f_type


! ─────────────────────────── calculate_RS_f ──────────────────────────────────
! R = unconstrained angular velocity; S = rotational compliance tensor.
! Side-effect: fills oo with O(orientation) via calculate_Jdot_f.
subroutine calculate_RS_f_type(RR, SS, oo, omegap, vortph, pgradh, o, ind, idfiber)
use param
use fibers
implicit none
double precision, intent(out)   :: RR(3), SS(3,3)
double precision, intent(inout) :: oo(3,3)
double precision, intent(in)    :: omegap(3), vortph(3), pgradh(9), o(3)
integer, intent(in) :: ind, idfiber

double precision :: JJ(3,3), JJdot(3,3), CC(3,3), HH(3), E(3,3)
double precision :: res(3,3), invres(3,3), temp(3), temp2(3), torque(3)
double precision :: line
integer :: f

!for multifiber
if (fiber_types_flag == 1) then
f = min((idfiber-1)/type_idx + 1, n_fiber_types)
else 
f = 0
endif

call calculate_J_f_type(JJ, o, f,ind) 
E(1,1) = pgradh(1);          E(1,2) = 0.5d0*(pgradh(2)+pgradh(4))
E(1,3) = 0.5d0*(pgradh(3)+pgradh(7))
E(2,1) = E(1,2);              E(2,2) = pgradh(5)
E(2,3) = 0.5d0*(pgradh(6)+pgradh(8))
E(3,1) = E(1,3);              E(3,2) = E(2,3);  E(3,3) = pgradh(9)

call calculate_H_f_type(HH, E, o, f, ind)
call calculate_C_f_type(CC, o, f, ind)
call calculate_Jdot_f_type(JJdot, oo, omegap, o, f, ind)   ! also fills oo = O(o)
call calculate_torque_f_type(torque, o, ind, idfiber, f)
! Store bending torque for output (cols 35-37 of fiber_*.dat)
torque_bend(1+(idfiber-1)*lambda+(ind-1), :) = torque

res = JJ + dt_part * (JJdot + CC)
call INVERSE_MATRIX_ff(res, invres)

call vec_product_f(JJ, omegap, temp2)
call vec_product_f(CC, vortph, temp)
temp2 = temp2 + dt_part*(temp + HH) + torque
call vec_product_f(invres, temp2, RR)

if (fiber_types_flag == 1 ) then
   line = dt_part / (pi * a_par_type(f)**4.d0 * dens_part_type(f))
elseif (bead_types_flag == 1) then
   if (ind == lambda) then
      line = dt_part / (pi * a_par_lb**4.d0 * dens_part)
   else
      line = dt_part / (pi * a_par**4.d0 * dens_part)
   endif
endif

SS   = line * invres
end subroutine calculate_RS_f_type


! ─────────────────────────── calculate_torque_f ──────────────────────────────
subroutine calculate_torque_f_type(torque, o, ind, idfiber, f)
use param
use fibers
implicit none
double precision, intent(out) :: torque(3)
double precision, intent(in)  :: o(3)
integer, intent(in) :: ind, idfiber,f

double precision :: oback(3), onext(3), ebend(3), ebendnext(3)
double precision :: coeff, alpha1, alpha2, norm
integer :: segno

torque = 0.d0
if (fiber_types_flag == 1)then
   if (E_bend_type(f) == 0.d0) return
else 
   if (E_bend == 0.0d0) return 
endif
! Neighbour orientations (global arrays, single-node)
if (ind == 1) then
   oback = 0.d0
else
   oback = angles_f(1+(idfiber-1)*lambda+(ind-2),:)
endif
if (ind == lambda) then
   onext = 0.d0
else
   onext = angles_f(1+(idfiber-1)*lambda+ind,:)
endif

segno = merge(1, -1, ind <= lambda/2)
ebendnext(1) = onext(2)*o(3) - o(2)*onext(3)
ebendnext(2) = -(onext(1)*o(3) - o(1)*onext(3))
ebendnext(3) = onext(1)*o(2) - o(1)*onext(2)
ebend(1) = o(2)*oback(3) - oback(2)*o(3)
ebend(2) = -(o(1)*oback(3) - oback(1)*o(3))
ebend(3) = o(1)*oback(2) - oback(1)*o(2)
if (fiber_types_flag == 1) then
   ebendnext(3) = ebendnext(3) + segno*curv_eq_type(f)
   ebend(3)     = ebend(3)     - segno*curv_eq_type(f)
else
   ebendnext(3) = ebendnext(3) + segno*curv_eq
   ebend(3)     = ebend(3)     - segno*curv_eq
endif
norm = sqrt(ebend(1)**2.d0+ebend(2)**2.d0+ebend(3)**2.d0)
if (norm > 1.d-9) then; ebend = ebend/norm; else; ebend = 0.d0; endif

norm = sqrt(ebendnext(1)**2.d0+ebendnext(2)**2.d0+ebendnext(3)**2.d0)

if (norm > 1.d-9) then; ebendnext = ebendnext/norm; else; ebendnext = 0.d0; endif

if (fiber_types_flag == 1) then
   coeff = -dt_part * E_bend_type(f) / &
            (8.d0*(arl_type(f)*a_par_type(f))**2.d0*dens_part_type(f))
elseif (bead_types_flag == 1) then
   if (ind == lambda) then
      coeff = -dt_part * E_bend / &
            (8.d0*(arl_lb*a_par_lb)**2.d0*dens_part)
   else
      coeff = -dt_part * E_bend / (8.d0*(arl*a_par)**2.d0*dens_part)
   endif
endif

alpha1 = max(-1.d0, min(1.d0, o(1)*oback(1)+o(2)*oback(2)+o(3)*oback(3)))
alpha1 = acos(alpha1)
alpha2 = max(-1.d0, min(1.d0, onext(1)*o(1)+onext(2)*o(2)+onext(3)*o(3)))
alpha2 = acos(alpha2)

torque = coeff * (alpha2*ebendnext - alpha1*ebend)
end subroutine calculate_torque_f_type


! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm







! ─────────────────────────── calculate_torque_f ──────────────────────────────
subroutine calculate_torque_f(torque, o, ind, idfiber)
use param
use fibers
implicit none
double precision, intent(out) :: torque(3)
double precision, intent(in)  :: o(3)
integer, intent(in) :: ind, idfiber

double precision :: oback(3), onext(3), ebend(3), ebendnext(3)
double precision :: coeff, alpha1, alpha2, norm
integer :: segno

torque = 0.d0
if (E_bend == 0.d0) return

! Neighbour orientations (global arrays, single-node)
if (ind == 1) then
   oback = 0.d0
else
   oback = angles_f(1+(idfiber-1)*lambda+(ind-2),:)
endif
if (ind == lambda) then
   onext = 0.d0
else
   onext = angles_f(1+(idfiber-1)*lambda+ind,:)
endif

segno = merge(1, -1, ind <= lambda/2)

ebendnext(1) = onext(2)*o(3) - o(2)*onext(3)
ebendnext(2) = -(onext(1)*o(3) - o(1)*onext(3))
ebendnext(3) = onext(1)*o(2) - o(1)*onext(2)
ebend(1) = o(2)*oback(3) - oback(2)*o(3)
ebend(2) = -(o(1)*oback(3) - oback(1)*o(3))
ebend(3) = o(1)*oback(2) - oback(1)*o(2)

ebendnext(3) = ebendnext(3) + segno*curv_eq
ebend(3)     = ebend(3)     - segno*curv_eq

norm = sqrt(ebend(1)**2.d0+ebend(2)**2.d0+ebend(3)**2.d0)
if (norm > 1.d-9) then; ebend = ebend/norm; else; ebend = 0.d0; endif
norm = sqrt(ebendnext(1)**2.d0+ebendnext(2)**2.d0+ebendnext(3)**2.d0)
if (norm > 1.d-9) then; ebendnext = ebendnext/norm; else; ebendnext = 0.d0; endif

coeff = -dt_part * E_bend / (8.d0*(arl*a_par)**2.d0*dens_part)

alpha1 = max(-1.d0, min(1.d0, o(1)*oback(1)+o(2)*oback(2)+o(3)*oback(3)))
alpha1 = acos(alpha1)
alpha2 = max(-1.d0, min(1.d0, onext(1)*o(1)+onext(2)*o(2)+onext(3)*o(3)))
alpha2 = acos(alpha2)

torque = coeff * (alpha2*ebendnext - alpha1*ebend)
end subroutine calculate_torque_f


! ─────────────────────────── fiber-fiber collision ───────────────────────────
subroutine fibre_collision_check_f(xx0,xx1,ang0,ang1,count, &
     fid0,fid1,rid0,rid1, fvxx0,fvxx1,fang0,fang1,fomg0,fomg1, &
     fvxx0_add,fvxx1_add,fomg0_add,fomg1_add,alreadycoll0,alreadycoll1, &
     ft0_add, ft1_add, torque_t0_add,torque_t1_add,stat_cnt,flag_f, & ! added for friction
     fttrq0_add,fttrq1_add,fr,fs) ! added for friction)
use param
use fibers
implicit none
double precision, intent(in)  :: xx0(3), xx1(3), ang0(3), ang1(3)
integer, intent(inout) :: count, alreadycoll0, alreadycoll1
integer, intent(in) :: fid0, fid1, rid0, rid1
double precision, intent(in)  :: fvxx0(lambda,3), fvxx1(lambda,3)
double precision, intent(in)  :: fang0(lambda,3), fang1(lambda,3)
double precision, intent(in)  :: fomg0(lambda,3), fomg1(lambda,3)
double precision, intent(out) :: fvxx0_add(lambda,3), fvxx1_add(lambda,3)
double precision, intent(out) :: fomg0_add(lambda,3), fomg1_add(lambda,3)
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! friction variables
double precision, dimension(3) :: relvel_cv,relvel_t, t_c  ! normal vectorial component, tangential vectorial and tangential module of relvel
double precision :: relvel_tm, r_c0m, r_c1m   ! term addressing normal force * tangential friction force
double precision, dimension(3) :: ft0_add, ft1_add, torque_t0_add, torque_t1_add,fttrq0_add,fttrq1_add
integer :: stat_cnt,flag_f
double precision :: thresh
! multifiber mmmmmmmmmm
integer, intent(in) :: fr,fs
double precision :: arl_cr0, a_par_cr0, arl_cr1, a_par_cr1
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

double precision :: dxx(3), dss(3), r_c0(3), r_c1(3), n_c(3)
double precision :: s0, s1, q0, q1, t0, t1, t2, dist, relvel_c
double precision :: vec1(3), vec2(3), vec3(3)
double precision :: RR_c0(3,3), RR_c1(3,3)

fvxx0_add = 0.d0; fvxx1_add = 0.d0
fomg0_add = 0.d0; fomg1_add = 0.d0

! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! initializing also the friction values - NB Ffric initialized and defined in fibers_module
! done two times as the quantities above
ft0_add       = 0.0d0
ft1_add       = 0.0d0
torque_t0_add = 0.0d0
torque_t1_add = 0.0d0
stat_cnt = 0
!threshold from static to dynamic
thresh = 1.0e-4
! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

! for multibead case - assigning arl and a_par the right values in fucntion of bead numbering
if ( bead_types_flag == 1) then 
   if (rid0 == lambda) then
      arl_cr0 = arl_lb 
      a_par_cr0 = a_par_lb
   else
      arl_cr0 = arl
      a_par_cr0 = a_par
   endif

   if (rid1 == lambda) then
      arl_cr1 = arl_lb 
      a_par_cr1 = a_par_lb
   else
      arl_cr1 = arl
      a_par_cr1 = a_par
   endif

else  ! to initialize such values in every condition
   arl_cr0 = arl; a_par_cr0 = a_par
   arl_cr1 = arl; a_par_cr1 = a_par
endif


dxx = xx0 - xx1
if (dxx(1) >  0.5d0*lx) dxx(1) = dxx(1) - lx
if (dxx(1) <= -0.5d0*lx) dxx(1) = dxx(1) + lx
if (dxx(2) >  0.5d0*lx) dxx(2) = dxx(2) - lx
if (dxx(2) <= -0.5d0*lx) dxx(2) = dxx(2) + lx
if (dxx(3) >  0.5d0*lx) dxx(3) = dxx(3) - lx
if (dxx(3) <= -0.5d0*lx) dxx(3) = dxx(3) + lx


if (fiber_types_flag == 1) then 
   !NBBBBBB s0 are points belonging to fiber rr and s1 are points belonging to fiber ss  
   ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
   t0 = sum(dxx*ang0);  t1 = sum(dxx*ang1);  t2 = sum(ang0*ang1)
   s0 = (t0 - t1*t2) / (a_par_type(fr)*(arl_type(fr)-1.d0)*min(t2**2.d0-1.d0, -1.d-12))
   s1 = (-t1 + t0*t2) / (a_par_type(fs)*(arl_type(fs)-1.d0)*min(t2**2.d0-1.d0, -1.d-12))
   if (abs(s0) > 1.d0 .and. abs(s1) <= 1.d0) then
      s0 = sign(1.d0,s0)
      s1 = (t1 + a_par_type(fs)*(arl_type(fs)-1.d0)*s0*t2) / & 
            (a_par_type(fs)*(arl_type(fs)-1.d0))
      if (abs(s1) > 1.d0) s1 = sign(1.d0,s1)
   elseif (abs(s0) <= 1.d0 .and. abs(s1) > 1.d0) then
      s1 = sign(1.d0,s1)
      s0 = (-t0 + a_par_type(fr)*(arl_type(fr)-1.d0)*s1*t2)/ &
            (a_par_type(fr)*(arl_type(fr)-1.d0))
      if (abs(s0) > 1.d0) s0 = sign(1.d0,s0)
   elseif (abs(s0) > 1.d0 .and. abs(s1) > 1.d0) then
      s0 = sign(1.d0,s0)  
      q1 = (t1+a_par_type(fs)*(arl_type(fs)-1.d0)*s0*t2)/&
      (a_par_type(fs)*(arl_type(fs)-1.d0))
      s1 = sign(1.d0,s1)
      q0 = (-t0+a_par_type(fr)*(arl_type(fr)-1.d0)*s1*t2)/ &
      (a_par_type(fr)*(arl_type(fr)-1.d0))
      if (abs(q0) <= 1.d0) then; s0 = q0
      elseif (abs(q1) <= 1.d0) then; s1 = q1
      else; s0 = sign(1.d0,q0); s1 = sign(1.d0,q1); endif
   endif

   dss = dxx + s0*a_par_type(fr)*(arl_type(fr)-1.d0)*ang0 - &
         s1*a_par_type(fs)*(arl_type(fs)-1.d0)*ang1
   dist = sqrt(sum(dss**2.d0))

   !if (dist <= 2.d0*a_par .and. alreadycoll0==0 .and. alreadycoll1==0) then
   if (dist <= (a_par_type(fr) + a_par_type(fs)) .and. alreadycoll0==0 .and. alreadycoll1==0) then
      n_c  = dss / dist
      r_c0 = a_par_type(fr)*(s0*(arl_type(fr)-1.d0)*ang0 - n_c)
      r_c1 = a_par_type(fs)*(s1*(arl_type(fs)-1.d0)*ang1 + n_c)
      call calculate_O_f(r_c0, RR_c0);  call calculate_O_f(r_c1, RR_c1)
      call vec_product_f(RR_c0, fomg0(rid0,:), vec1)
      call vec_product_f(RR_c1, fomg1(rid1,:), vec2)
      vec3    = fvxx0(rid0,:) - fvxx1(rid1,:) - vec1 + vec2
      relvel_c = sum(n_c*vec3)

      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! modification to calculate tangent delta velocity 
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      flag_f = 1
      relvel_cv = relvel_c * n_c
      !tangent and module
      relvel_t = vec3 - relvel_cv
      relvel_tm = sqrt(relvel_t(1)*relvel_t(1) + relvel_t(2)*relvel_t(2)+relvel_t(3)*relvel_t(3))
      if (relvel_tm .le. thresh) then 
         ! to avoid division by 0 and define also a static case flag
         t_c = 0.0d0
         stat_cnt = 1 
      else 
         ! calculate tangent versor
         t_c = relvel_t/relvel_tm
      endif 

      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! calculate the friction force to add as an output
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! in this case when particles in contact we always sum it because we just need them in contact 
      ft0_add = - Ffric * t_c
      ft1_add = Ffric * t_c    ! we have equal and opposite component on the other fiber
      ! friction torque - using the same RR matrix calculated above 
      call vec_product_f(RR_c0,ft0_add,torque_t0_add)  
      call vec_product_f(RR_c1,ft1_add,torque_t1_add)  
      ! calc of the distance module 
      r_c0m = sqrt(r_c0(1)*r_c0(1) + r_c0(2)*r_c0(2) + r_c0(3)*r_c0(3))
      r_c1m = sqrt(r_c1(1)*r_c1(1) + r_c1(2)*r_c1(2) + r_c1(3)*r_c1(3))
      ! rescaled force
      fttrq0_add = ft0_add*r_c0m/( arl_type(fr) * a_par_type(fr)) ! to see if actually arl*apiu(1) means the actual lenght
      fttrq1_add = ft1_add*r_c1m/(arl_type(fs) * a_par_type(fs))
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      if (relvel_c <= 0.d0) then
         count = 1
         if (fid0 == fid1) then
            call in_fibre_collision_solver_f_type(fvxx0,fang0,fomg0,rid0,rid1, &
               r_c0,r_c1,n_c,coef_res,relvel_c,fvxx0_add,fomg0_add,fr)
         else
            call ex_fibre_collision_solver_f_type(fvxx0,fvxx1,fang0,fang1, &
               fomg0,fomg1,rid0,rid1,r_c0,r_c1,n_c,coef_res,relvel_c, &
               fvxx0_add,fvxx1_add,fomg0_add,fomg1_add,fr,fs)
         endif
      endif
   endif

elseif (bead_types_flag == 1) then

   !NBBBBBB s0 are points belonging to fiber rr and s1 are points belonging to fiber ss
   ! we wrote in a way that arl_cr1 and cr0 will automatically input the right dimensions in function of the choosen bead  
   ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
   t0 = sum(dxx*ang0);  t1 = sum(dxx*ang1);  t2 = sum(ang0*ang1)
   s0 = (t0 - t1*t2) / (a_par_cr0*(arl_cr0-1.d0)*min(t2**2.d0-1.d0, -1.d-12))
   s1 = (-t1 + t0*t2) / (a_par_cr1*(arl_cr1-1.d0)*min(t2**2.d0-1.d0, -1.d-12))
   if (abs(s0) > 1.d0 .and. abs(s1) <= 1.d0) then
      s0 = sign(1.d0,s0)
      s1 = (t1 + a_par_cr1*(arl_cr1-1.d0)*s0*t2) / & 
            (a_par_cr1*(arl_cr1-1.d0))
      if (abs(s1) > 1.d0) s1 = sign(1.d0,s1)
   elseif (abs(s0) <= 1.d0 .and. abs(s1) > 1.d0) then
      s1 = sign(1.d0,s1)
      s0 = (-t0 + a_par_cr0*(arl_cr0-1.d0)*s1*t2)/ &
            (a_par_cr0*(arl_cr0-1.d0))
      if (abs(s0) > 1.d0) s0 = sign(1.d0,s0)
   elseif (abs(s0) > 1.d0 .and. abs(s1) > 1.d0) then
      s0 = sign(1.d0,s0)  
      q1 = (t1+a_par_cr1*(arl_cr1-1.d0)*s0*t2)/&
      (a_par_cr1*(arl_cr1-1.d0))
      s1 = sign(1.d0,s1)
      q0 = (-t0+a_par_cr0*(arl_cr0-1.d0)*s1*t2)/ &
      (a_par_cr0*(arl_cr0-1.d0))
      if (abs(q0) <= 1.d0) then; s0 = q0
      elseif (abs(q1) <= 1.d0) then; s1 = q1
      else; s0 = sign(1.d0,q0); s1 = sign(1.d0,q1); endif
   endif

   dss = dxx + s0*a_par_cr0*(arl_cr0-1.d0)*ang0 - &
         s1*a_par_cr1*(arl_cr1-1.d0)*ang1
   dist = sqrt(sum(dss**2.d0))

   !if (dist <= 2.d0*a_par .and. alreadycoll0==0 .and. alreadycoll1==0) then
   if (dist <= (a_par_cr0 + a_par_cr1) .and. alreadycoll0==0 .and. alreadycoll1==0) then
      n_c  = dss / dist
      r_c0 = a_par_cr0*(s0*(arl_cr0-1.d0)*ang0 - n_c)
      r_c1 = a_par_cr1*(s1*(arl_cr1-1.d0)*ang1 + n_c)
      call calculate_O_f(r_c0, RR_c0);  call calculate_O_f(r_c1, RR_c1)
      call vec_product_f(RR_c0, fomg0(rid0,:), vec1)
      call vec_product_f(RR_c1, fomg1(rid1,:), vec2)
      vec3    = fvxx0(rid0,:) - fvxx1(rid1,:) - vec1 + vec2
      relvel_c = sum(n_c*vec3)

      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! modification to calculate tangent delta velocity 
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      flag_f = 1
      relvel_cv = relvel_c * n_c
      !tangent and module
      relvel_t = vec3 - relvel_cv
      relvel_tm = sqrt(relvel_t(1)*relvel_t(1) + relvel_t(2)*relvel_t(2)+relvel_t(3)*relvel_t(3))
      if (relvel_tm .le. thresh) then 
         ! to avoid division by 0 and define also a static case flag
         t_c = 0.0d0
         stat_cnt = 1 
      else 
         ! calculate tangent versor
         t_c = relvel_t/relvel_tm
      endif 

      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! calculate the friction force to add as an output
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! in this case when particles in contact we always sum it because we just need them in contact 
      ft0_add = - Ffric * t_c
      ft1_add = Ffric * t_c    ! we have equal and opposite component on the other fiber
      ! friction torque - using the same RR matrix calculated above 
      call vec_product_f(RR_c0,ft0_add,torque_t0_add)  
      call vec_product_f(RR_c1,ft1_add,torque_t1_add)  
      ! calc of the distance module 
      r_c0m = sqrt(r_c0(1)*r_c0(1) + r_c0(2)*r_c0(2) + r_c0(3)*r_c0(3))
      r_c1m = sqrt(r_c1(1)*r_c1(1) + r_c1(2)*r_c1(2) + r_c1(3)*r_c1(3))
      ! rescaled force
      fttrq0_add = ft0_add*r_c0m/( arl_cr0 * a_par_cr0) ! to see if actually arl*apiu(1) means the actual lenght
      fttrq1_add = ft1_add*r_c1m/(arl_cr1 * a_par_cr1)
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      if (relvel_c <= 0.d0) then
         count = 1
         if (fid0 == fid1) then
            call in_fibre_collision_solver_f_type(fvxx0,fang0,fomg0,rid0,rid1, &
               r_c0,r_c1,n_c,coef_res,relvel_c,fvxx0_add,fomg0_add,fr)
         else
            call ex_fibre_collision_solver_f_type(fvxx0,fvxx1,fang0,fang1, &
               fomg0,fomg1,rid0,rid1,r_c0,r_c1,n_c,coef_res,relvel_c, &
               fvxx0_add,fvxx1_add,fomg0_add,fomg1_add,fr,fs)
         endif
      endif
   endif
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
else
   t0 = sum(dxx*ang0);  t1 = sum(dxx*ang1);  t2 = sum(ang0*ang1)
   s0 = (t0 - t1*t2) / (a_par*(arl-1.d0)*min(t2**2.d0-1.d0, -1.d-12))
   s1 = (-t1 + t0*t2) / (a_par*(arl-1.d0)*min(t2**2.d0-1.d0, -1.d-12))

   if (abs(s0) > 1.d0 .and. abs(s1) <= 1.d0) then
      s0 = sign(1.d0,s0)
      s1 = (t1 + a_par*(arl-1.d0)*s0*t2) / (a_par*(arl-1.d0))
      if (abs(s1) > 1.d0) s1 = sign(1.d0,s1)
   elseif (abs(s0) <= 1.d0 .and. abs(s1) > 1.d0) then
      s1 = sign(1.d0,s1)
      s0 = (-t0 + a_par*(arl-1.d0)*s1*t2) / (a_par*(arl-1.d0))
      if (abs(s0) > 1.d0) s0 = sign(1.d0,s0)
   elseif (abs(s0) > 1.d0 .and. abs(s1) > 1.d0) then
      s0 = sign(1.d0,s0);  q1 = (t1+a_par*(arl-1.d0)*s0*t2)/(a_par*(arl-1.d0))
      s1 = sign(1.d0,s1);  q0 = (-t0+a_par*(arl-1.d0)*s1*t2)/(a_par*(arl-1.d0))
      if (abs(q0) <= 1.d0) then; s0 = q0
      elseif (abs(q1) <= 1.d0) then; s1 = q1
      else; s0 = sign(1.d0,q0); s1 = sign(1.d0,q1); endif
   endif

   dss = dxx + s0*a_par*(arl-1.d0)*ang0 - s1*a_par*(arl-1.d0)*ang1
   dist = sqrt(sum(dss**2.d0))

   if (dist <= 2.d0*a_par .and. alreadycoll0==0 .and. alreadycoll1==0) then
      n_c  = dss / dist
      r_c0 = a_par*(s0*(arl-1.d0)*ang0 - n_c)
      r_c1 = a_par*(s1*(arl-1.d0)*ang1 + n_c)
      call calculate_O_f(r_c0, RR_c0);  call calculate_O_f(r_c1, RR_c1)
      call vec_product_f(RR_c0, fomg0(rid0,:), vec1)
      call vec_product_f(RR_c1, fomg1(rid1,:), vec2)
      vec3    = fvxx0(rid0,:) - fvxx1(rid1,:) - vec1 + vec2
      relvel_c = sum(n_c*vec3)

      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! modification to calculate tangent delta velocity 
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      flag_f = 1

      relvel_cv = relvel_c * n_c
      !tangent
      relvel_t = vec3 - relvel_cv
      ! tangent module
      relvel_tm = sqrt(relvel_t(1)*relvel_t(1) + relvel_t(2)*relvel_t(2)+relvel_t(3)*relvel_t(3))
      
      if (relvel_tm .le. thresh) then 
         ! to avoid division by 0 and define also a static case flag
         t_c = 0.0d0
         stat_cnt = 1 
      else 
         ! calculate tangent versor
         t_c = relvel_t/relvel_tm
      endif 
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! calculate the friction force to add as an output
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ! in this case when particles in contact we always sum it because we just need them in contact 
      ft0_add = - Ffric * t_c
      ft1_add = Ffric * t_c    ! we have equal and opposite component on the other fiber
      ! friction torque - using the same RR matrix calculated above 
      call vec_product_f(RR_c0,ft0_add,torque_t0_add)  
      call vec_product_f(RR_c1,ft1_add,torque_t1_add)  
      ! NB in the force and moments eq. actually internal forces applied on the center of bead so adim.
      ! To see in function of how X_0 is defined IMPPPPPPPPPPPPPPP
      ! calc of the distance module 
      r_c0m = sqrt(r_c0(1)*r_c0(1) + r_c0(2)*r_c0(2) + r_c0(3)*r_c0(3))
      r_c1m = sqrt(r_c1(1)*r_c1(1) + r_c1(2)*r_c1(2) + r_c1(3)*r_c1(3))
      ! rescaled force
      fttrq0_add = ft0_add*r_c0m/( arl * a_par) ! to see if actually arl*apiu(1) means the actual lenght
      fttrq1_add = ft1_add*r_c1m/(arl* a_par)
      ! zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

      if (relvel_c <= 0.d0) then
         count = 1
         if (fid0 == fid1) then
            call in_fibre_collision_solver_f(fvxx0,fang0,fomg0,rid0,rid1, &
               r_c0,r_c1,n_c,coef_res,relvel_c,fvxx0_add,fomg0_add)
         else
            call ex_fibre_collision_solver_f(fvxx0,fvxx1,fang0,fang1, &
               fomg0,fomg1,rid0,rid1,r_c0,r_c1,n_c,coef_res,relvel_c, &
               fvxx0_add,fvxx1_add,fomg0_add,fomg1_add)
         endif
      endif
   endif
endif
end subroutine fibre_collision_check_f !----------------- modificata ma da vedere start from here to changeeeee NBBBBBBB


! ─── in_fibre_collision_solver_f ──────────────────────────────────────────────
subroutine in_fibre_collision_solver_f(vxx,ang,omg,j_c,k_c,r_cj,r_ck,n_c, &
                                       coef_res,relvel_c,vxx_add,omg_add)
use fibers
implicit none
integer, intent(in) :: j_c, k_c
double precision, intent(in)  :: r_cj(3), r_ck(3), n_c(3), coef_res, relvel_c
double precision, intent(in)  :: vxx(lambda,3), ang(lambda,3), omg(lambda,3)
double precision, intent(out) :: vxx_add(lambda,3), omg_add(lambda,3)

double precision :: OO(lambda,3,3), JJ(lambda,3,3), invJJ(lambda,3,3)
double precision :: DD(3,3,lambda-1), UU(3,3,lambda-2), LL(3,3,lambda-2)
double precision :: RHS_f(3,lambda-1), VV_c(3,lambda-1)
double precision :: vec_array1(3,lambda-1), vec_array2(3,lambda-1)
double precision :: II(3,lambda+1), RR_cj(3,3), RR_ck(3,3)
double precision :: mat1(3,3), mat2(3,3), mat3(3,3), mat4(3,3), mat5(3,3), mat6(3,3)
double precision :: vec1(3), vec2(3)
double precision :: SS_c, FF_c, II_c, tempa, tempb, tempc, tempd
integer :: j

do j = 1, lambda
   call calculate_O_f(ang(j,:), OO(j,:,:))
   call calculate_J_f(JJ(j,:,:), ang(j,:))
   JJ(j,:,:) = 0.5d0*a_par**2.d0*JJ(j,:,:)
   call INVERSE_MATRIX_ff(JJ(j,:,:), invJJ(j,:,:))
   if (j > 1) then
      call calculate_prod_array_f(OO(j-1,:,:),invJJ(j-1,:,:),mat1)
      call calculate_prod_array_f(mat1,OO(j-1,:,:),mat2)
      call calculate_prod_array_f(OO(j,:,:),invJJ(j,:,:),mat3)
      call calculate_prod_array_f(mat3,OO(j,:,:),mat4)
      DD(:,:,j-1) = 2.d0*Idm - (arl*a_par)**2.d0*(mat2+mat4)
      if (j < lambda) then
         UU(:,:,j-1) = -Idm - (arl*a_par)**2.d0*mat4
         LL(:,:,j-1) = UU(:,:,j-1)
      endif
      call vec_product_f(OO(j-1,:,:),omg(j-1,:),vec1)
      call vec_product_f(OO(j,:,:),omg(j,:),vec2)
      RHS_f(:,j-1) = -vxx(j-1,:)+vxx(j,:)+arl*a_par*(vec1+vec2)
   endif
enddo

VV_c = 0.d0
call calculate_O_f(r_cj, RR_cj)
call calculate_prod_array_f(RR_cj,invJJ(j_c,:,:),mat1)
call calculate_prod_array_f(mat1,OO(j_c,:,:),mat2)
if (j_c > 1) then
   mat4 = -Idm-arl*a_par*mat2
   call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1))
endif
if (j_c < lambda) then
   mat5 = Idm-arl*a_par*mat2
   call vec_product_left_f(n_c,mat5,VV_c(:,j_c))
endif
call calculate_prod_array_f(mat1,RR_cj,mat3); mat6 = Idm-mat3

call calculate_O_f(r_ck, RR_ck)
call calculate_prod_array_f(RR_ck,invJJ(k_c,:,:),mat1)
call calculate_prod_array_f(mat1,OO(k_c,:,:),mat2)
if (k_c > 1) then
   mat4 = Idm+arl*a_par*mat2
   call vec_product_left_f(n_c,mat4,VV_c(:,k_c-1))
endif
if (k_c < lambda) then
   mat5 = -Idm+arl*a_par*mat2
   call vec_product_left_f(n_c,mat5,VV_c(:,k_c))
endif
call calculate_prod_array_f(mat1,RR_ck,mat3); mat6 = mat6+Idm-mat3

call vec_product_f(mat6,n_c,vec1); SS_c = sum(n_c*vec1)
FF_c = -(1.d0+coef_res)*relvel_c

call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,vec_array1)
call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,VV_c,vec_array2)

II = 0.d0; tempc = 0.d0; tempd = 0.d0
do j = 1, lambda-1
   tempc = tempc + sum(VV_c(:,j)*vec_array1(:,j))
   tempd = tempd + sum(VV_c(:,j)*vec_array2(:,j))
enddo
do j = 1, lambda-1
   II(:,j+1) = vec_array1(:,j) + tempc/(SS_c-tempd)*vec_array2(:,j)
enddo
tempc = 0.d0
do j = 1, lambda-1; tempc = tempc + sum(VV_c(:,j)*II(:,j+1)); enddo
II_c = (FF_c-tempc)/SS_c

do j = 1, lambda
   vxx_add(j,:) = II(:,j+1) - II(:,j)
   call vec_product_f(OO(j,:,:),II(:,j+1)+II(:,j),vec1)
   call vec_product_f(invJJ(j,:,:),vec1,vec2)
   omg_add(j,:) = arl*a_par*vec2
enddo
vxx_add(j_c,:) = vxx_add(j_c,:) + II_c*n_c
call vec_product_f(RR_cj,n_c,vec1)
call vec_product_f(invJJ(j_c,:,:),vec1,vec2)
omg_add(j_c,:) = omg_add(j_c,:) + II_c*vec2
vxx_add(k_c,:) = vxx_add(k_c,:) - II_c*n_c
call vec_product_f(RR_ck,n_c,vec1)
call vec_product_f(invJJ(k_c,:,:),vec1,vec2)
omg_add(k_c,:) = omg_add(k_c,:) - II_c*vec2
end subroutine in_fibre_collision_solver_f


! ─── ex_fibre_collision_solver_f ──────────────────────────────────────────────
subroutine ex_fibre_collision_solver_f(vxx0,vxx1,ang0,ang1,omg0,omg1, &
     j_c,k_c,r_cj,r_ck,n_c,coef_res,relvel_c,vxx0_add,vxx1_add,omg0_add,omg1_add)
use fibers
implicit none
integer, intent(in) :: j_c, k_c
double precision, intent(in)  :: r_cj(3), r_ck(3), n_c(3), coef_res, relvel_c
double precision, intent(in)  :: vxx0(lambda,3), vxx1(lambda,3)
double precision, intent(in)  :: ang0(lambda,3), ang1(lambda,3)
double precision, intent(in)  :: omg0(lambda,3), omg1(lambda,3)
double precision, intent(out) :: vxx0_add(lambda,3), vxx1_add(lambda,3)
double precision, intent(out) :: omg0_add(lambda,3), omg1_add(lambda,3)

! dim2 = 2*lambda-2 is now a runtime quantity (lambda is a module variable, not parameter)
double precision :: OO0(lambda,3,3), invJJ0(lambda,3,3), JJ0(lambda,3,3)
double precision :: OO1(lambda,3,3), invJJ1(lambda,3,3), JJ1(lambda,3,3)
double precision :: DD(3,3,2*lambda-2), UU(3,3,2*lambda-3), LL(3,3,2*lambda-3)
double precision :: RHS_f(3,2*lambda-2), VV_c(3,2*lambda-2)
double precision :: va1(3,2*lambda-2), va2(3,2*lambda-2)
double precision :: II0(3,lambda+1), II1(3,lambda+1)
double precision :: RR_cj(3,3), RR_ck(3,3)
double precision :: mat1(3,3), mat2(3,3), mat3(3,3), mat4(3,3), mat5(3,3), mat6(3,3)
double precision :: vec1(3), vec2(3)
double precision :: SS_c, FF_c, II_c, tempa, tempb, tempc, tempd
integer :: j, k

do j = 1, lambda
   call calculate_O_f(ang0(j,:),OO0(j,:,:))
   call calculate_J_f(JJ0(j,:,:),ang0(j,:));  JJ0(j,:,:)=0.5d0*a_par**2.d0*JJ0(j,:,:)
   call INVERSE_MATRIX_ff(JJ0(j,:,:),invJJ0(j,:,:))
   if (j > 1) then
      call calculate_prod_array_f(OO0(j-1,:,:),invJJ0(j-1,:,:),mat1)
      call calculate_prod_array_f(mat1,OO0(j-1,:,:),mat2)
      call calculate_prod_array_f(OO0(j,:,:),invJJ0(j,:,:),mat3)
      call calculate_prod_array_f(mat3,OO0(j,:,:),mat4)
      DD(:,:,j-1) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2+mat4)
      if (j < lambda) then; UU(:,:,j-1)=-Idm-(arl*a_par)**2.d0*mat4; LL(:,:,j-1)=UU(:,:,j-1); endif
      call vec_product_f(OO0(j-1,:,:),omg0(j-1,:),vec1)
      call vec_product_f(OO0(j,:,:),omg0(j,:),vec2)
      RHS_f(:,j-1) = -vxx0(j-1,:)+vxx0(j,:)+arl*a_par*(vec1+vec2)
   endif
enddo
do k = 1, lambda
   call calculate_O_f(ang1(k,:),OO1(k,:,:))
   call calculate_J_f(JJ1(k,:,:),ang1(k,:));  JJ1(k,:,:)=0.5d0*a_par**2.d0*JJ1(k,:,:)
   call INVERSE_MATRIX_ff(JJ1(k,:,:),invJJ1(k,:,:))
   if (k > 1) then
      call calculate_prod_array_f(OO1(k-1,:,:),invJJ1(k-1,:,:),mat1)
      call calculate_prod_array_f(mat1,OO1(k-1,:,:),mat2)
      call calculate_prod_array_f(OO1(k,:,:),invJJ1(k,:,:),mat3)
      call calculate_prod_array_f(mat3,OO1(k,:,:),mat4)
      DD(:,:,k+lambda-2) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2+mat4)
      if (k < lambda) then; UU(:,:,k+lambda-2)=-Idm-(arl*a_par)**2.d0*mat4; LL(:,:,k+lambda-2)=UU(:,:,k+lambda-2); endif
      call vec_product_f(OO1(k-1,:,:),omg1(k-1,:),vec1)
      call vec_product_f(OO1(k,:,:),omg1(k,:),vec2)
      RHS_f(:,k+lambda-2) = -vxx1(k-1,:)+vxx1(k,:)+arl*a_par*(vec1+vec2)
   endif
enddo
UU(:,:,lambda-1) = 0.d0;  LL(:,:,lambda-1) = 0.d0

VV_c = 0.d0
call calculate_O_f(r_cj,RR_cj)
call calculate_prod_array_f(RR_cj,invJJ0(j_c,:,:),mat1)
call calculate_prod_array_f(mat1,OO0(j_c,:,:),mat2)
if (j_c > 1) then; mat4=-Idm-arl*a_par*mat2; call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1)); endif
if (j_c < lambda) then; mat5=Idm-arl*a_par*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,j_c)); endif
call calculate_prod_array_f(mat1,RR_cj,mat3); mat6=Idm-mat3

call calculate_O_f(r_ck,RR_ck)
call calculate_prod_array_f(RR_ck,invJJ1(k_c,:,:),mat1)
call calculate_prod_array_f(mat1,OO1(k_c,:,:),mat2)
if (k_c > 1) then; mat4=Idm+arl*a_par*mat2; call vec_product_left_f(n_c,mat4,VV_c(:,k_c+lambda-2)); endif
if (k_c < lambda) then; mat5=-Idm+arl*a_par*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,k_c+lambda-1)); endif
call calculate_prod_array_f(mat1,RR_ck,mat3); mat6=mat6+Idm-mat3

call vec_product_f(mat6,n_c,vec1); SS_c=sum(n_c*vec1)
FF_c = -(1.d0+coef_res)*relvel_c

call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,va1)
call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,VV_c,va2)

II0=0.d0; II1=0.d0; tempc=0.d0; tempd=0.d0
do j=1,2*lambda-2; tempc=tempc+sum(VV_c(:,j)*va1(:,j)); tempd=tempd+sum(VV_c(:,j)*va2(:,j)); enddo
do j=1,lambda-1; II0(:,j+1)=va1(:,j)+tempc/(SS_c-tempd)*va2(:,j); enddo
do k=1,lambda-1; II1(:,k+1)=va1(:,k+lambda-1)+tempc/(SS_c-tempd)*va2(:,k+lambda-1); enddo
tempc=0.d0
do j=1,lambda-1; tempc=tempc+sum(VV_c(:,j)*II0(:,j+1)); enddo
do k=1,lambda-1; tempc=tempc+sum(VV_c(:,k+lambda-1)*II1(:,k+1)); enddo
II_c=(FF_c-tempc)/SS_c

do j=1,lambda
   vxx0_add(j,:)=II0(:,j+1)-II0(:,j)
   call vec_product_f(OO0(j,:,:),II0(:,j+1)+II0(:,j),vec1)
   call vec_product_f(invJJ0(j,:,:),vec1,vec2); omg0_add(j,:)=arl*a_par*vec2
enddo
do k=1,lambda
   vxx1_add(k,:)=II1(:,k+1)-II1(:,k)
   call vec_product_f(OO1(k,:,:),II1(:,k+1)+II1(:,k),vec1)
   call vec_product_f(invJJ1(k,:,:),vec1,vec2); omg1_add(k,:)=arl*a_par*vec2
enddo
vxx0_add(j_c,:)=vxx0_add(j_c,:)+II_c*n_c
call vec_product_f(RR_cj,n_c,vec1); call vec_product_f(invJJ0(j_c,:,:),vec1,vec2)
omg0_add(j_c,:)=omg0_add(j_c,:)+II_c*vec2
vxx1_add(k_c,:)=vxx1_add(k_c,:)-II_c*n_c
call vec_product_f(RR_ck,n_c,vec1); call vec_product_f(invJJ1(k_c,:,:),vec1,vec2)
omg1_add(k_c,:)=omg1_add(k_c,:)-II_c*vec2
end subroutine ex_fibre_collision_solver_f



! multifiber versionssss  mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! ─── in_fibre_collision_solver_f_type──────────────────────────────────────────────
subroutine in_fibre_collision_solver_f_type(vxx,ang,omg,j_c,k_c,r_cj,r_ck,n_c, &
                                       coef_res,relvel_c,vxx_add,omg_add,fr)
use fibers
implicit none
integer, intent(in) :: j_c, k_c, fr
double precision, intent(in)  :: r_cj(3), r_ck(3), n_c(3), coef_res, relvel_c
double precision, intent(in)  :: vxx(lambda,3), ang(lambda,3), omg(lambda,3)
double precision, intent(out) :: vxx_add(lambda,3), omg_add(lambda,3)

double precision :: OO(lambda,3,3), JJ(lambda,3,3), invJJ(lambda,3,3)
double precision :: DD(3,3,lambda-1), UU(3,3,lambda-2), LL(3,3,lambda-2)
double precision :: RHS_f(3,lambda-1), VV_c(3,lambda-1)
double precision :: vec_array1(3,lambda-1), vec_array2(3,lambda-1)
double precision :: II(3,lambda+1), RR_cj(3,3), RR_ck(3,3)
double precision :: mat1(3,3), mat2(3,3), mat3(3,3), mat4(3,3), mat5(3,3), mat6(3,3)
double precision :: vec1(3), vec2(3)
double precision :: SS_c, FF_c, II_c, tempa, tempb, tempc, tempd
integer :: j
!  multibead 
double precision :: arl_cr, a_par_cr, arl_pr, a_par_pr


if (fiber_types_flag == 1) then 

   do j = 1, lambda
      call calculate_O_f(ang(j,:), OO(j,:,:))
      call calculate_J_f_type(JJ(j,:,:), ang(j,:),fr)
      JJ(j,:,:) = 0.5d0*a_par_type(fr)**2.d0*JJ(j,:,:)
      call INVERSE_MATRIX_ff(JJ(j,:,:), invJJ(j,:,:))
      if (j > 1) then
         call calculate_prod_array_f(OO(j-1,:,:),invJJ(j-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO(j-1,:,:),mat2)
         call calculate_prod_array_f(OO(j,:,:),invJJ(j,:,:),mat3)
         call calculate_prod_array_f(mat3,OO(j,:,:),mat4)
         DD(:,:,j-1) = 2.d0*Idm - (arl_type(fr)*a_par_type(fr))**2.d0*(mat2+mat4)
         if (j < lambda) then
            UU(:,:,j-1) = -Idm - (arl_type(fr)*a_par_type(fr))**2.d0*mat4
            LL(:,:,j-1) = UU(:,:,j-1)
         endif
         call vec_product_f(OO(j-1,:,:),omg(j-1,:),vec1)
         call vec_product_f(OO(j,:,:),omg(j,:),vec2)
         RHS_f(:,j-1) = -vxx(j-1,:)+vxx(j,:)+arl_type(fr)*a_par_type(fr)*(vec1+vec2)
      endif
   enddo

   VV_c = 0.d0
   call calculate_O_f(r_cj, RR_cj)
   call calculate_prod_array_f(RR_cj,invJJ(j_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO(j_c,:,:),mat2)
   if (j_c > 1) then
      mat4 = -Idm-arl_type(fr)*a_par_type(fr)*mat2
      call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1))
   endif
   if (j_c < lambda) then
      mat5 = Idm-arl_type(fr)*a_par_type(fr)*mat2
      call vec_product_left_f(n_c,mat5,VV_c(:,j_c))
   endif
   call calculate_prod_array_f(mat1,RR_cj,mat3); mat6 = Idm-mat3

   call calculate_O_f(r_ck, RR_ck)
   call calculate_prod_array_f(RR_ck,invJJ(k_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO(k_c,:,:),mat2)
   if (k_c > 1) then
      mat4 = Idm+arl_type(fr)*a_par_type(fr)*mat2
      call vec_product_left_f(n_c,mat4,VV_c(:,k_c-1))
   endif
   if (k_c < lambda) then
      mat5 = -Idm+arl_type(fr)*a_par_type(fr)*mat2
      call vec_product_left_f(n_c,mat5,VV_c(:,k_c))
   endif
   call calculate_prod_array_f(mat1,RR_ck,mat3); mat6 = mat6+Idm-mat3

   call vec_product_f(mat6,n_c,vec1); SS_c = sum(n_c*vec1)
   FF_c = -(1.d0+coef_res)*relvel_c

   call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,vec_array1)
   call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,VV_c,vec_array2)

   II = 0.d0; tempc = 0.d0; tempd = 0.d0
   do j = 1, lambda-1
      tempc = tempc + sum(VV_c(:,j)*vec_array1(:,j))
      tempd = tempd + sum(VV_c(:,j)*vec_array2(:,j))
   enddo
   do j = 1, lambda-1
      II(:,j+1) = vec_array1(:,j) + tempc/(SS_c-tempd)*vec_array2(:,j)
   enddo
   tempc = 0.d0
   do j = 1, lambda-1; tempc = tempc + sum(VV_c(:,j)*II(:,j+1)); enddo
   II_c = (FF_c-tempc)/SS_c

   do j = 1, lambda
      vxx_add(j,:) = II(:,j+1) - II(:,j)
      call vec_product_f(OO(j,:,:),II(:,j+1)+II(:,j),vec1)
      call vec_product_f(invJJ(j,:,:),vec1,vec2)
      omg_add(j,:) = arl_type(fr)*a_par_type(fr)*vec2
   enddo
   vxx_add(j_c,:) = vxx_add(j_c,:) + II_c*n_c
   call vec_product_f(RR_cj,n_c,vec1)
   call vec_product_f(invJJ(j_c,:,:),vec1,vec2)
   omg_add(j_c,:) = omg_add(j_c,:) + II_c*vec2
   vxx_add(k_c,:) = vxx_add(k_c,:) - II_c*n_c
   call vec_product_f(RR_ck,n_c,vec1)
   call vec_product_f(invJJ(k_c,:,:),vec1,vec2)
   omg_add(k_c,:) = omg_add(k_c,:) - II_c*vec2

elseif ( bead_types_flag == 1) then

   do j = 1, lambda
      ! to iterate with the right a_par and arl values-----
      if (j == lambda) then 
         arl_cr = arl_lb; a_par_cr = a_par_lb
         arl_pr = arl; a_par_pr = a_par
      else
         arl_cr = arl; a_par_cr = a_par
         arl_pr = arl; a_par_pr = a_par
      endif
      !continuing the routine ----------------------------
      call calculate_O_f(ang(j,:), OO(j,:,:))
      call calculate_J_f_type(JJ(j,:,:), ang(j,:),fr,j)
      JJ(j,:,:) = 0.5d0*a_par_cr**2.d0*JJ(j,:,:)
      call INVERSE_MATRIX_ff(JJ(j,:,:), invJJ(j,:,:))
      if (j > 1) then
         call calculate_prod_array_f(OO(j-1,:,:),invJJ(j-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO(j-1,:,:),mat2)
         !NB then mat2 is associated to bead j-1 !!
         call calculate_prod_array_f(OO(j,:,:),invJJ(j,:,:),mat3)
         call calculate_prod_array_f(mat3,OO(j,:,:),mat4)
         !mat4 associated to j-th bead
         if (j == lambda) then
            DD(:,:,j-1) = 2.d0*Idm - (arl*a_par)**2.d0*(mat2) &
                                   - (arl_cr*a_par_cr)**2.d0*(mat4)
         else
            DD(:,:,j-1) = 2.d0*Idm - (arl_cr*a_par_cr)**2.d0*(mat2+mat4)
         endif
         if (j < lambda) then
            UU(:,:,j-1) = -Idm - (arl_cr*a_par_cr)**2.d0*mat4
            LL(:,:,j-1) = UU(:,:,j-1)
         endif
         call vec_product_f(OO(j-1,:,:),omg(j-1,:),vec1)    !NB also here on j-1
         call vec_product_f(OO(j,:,:),omg(j,:),vec2)
         if (j == lambda) then 
            RHS_f(:,j-1) = -vxx(j-1,:)+vxx(j,:)+arl*a_par*(vec1) &
                                               +arl_cr*a_par_cr*(vec2)
         else
            RHS_f(:,j-1) = -vxx(j-1,:)+vxx(j,:)+arl_cr*a_par_cr*(vec1+vec2)
         endif
      endif
   enddo

   VV_c = 0.d0
   call calculate_O_f(r_cj, RR_cj)
   call calculate_prod_array_f(RR_cj,invJJ(j_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO(j_c,:,:),mat2)
   if (j_c > 1) then
      if (j_c == lambda) then
         mat4 = -Idm-arl_lb*a_par_lb*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1))
      else 
         mat4 = -Idm-arl*a_par*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1)) 
      endif
   endif
   if (j_c < lambda) then
      mat5 = Idm-arl*a_par*mat2
      call vec_product_left_f(n_c,mat5,VV_c(:,j_c))
   endif
   call calculate_prod_array_f(mat1,RR_cj,mat3); mat6 = Idm-mat3

   call calculate_O_f(r_ck, RR_ck)
   call calculate_prod_array_f(RR_ck,invJJ(k_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO(k_c,:,:),mat2)
   if (k_c > 1) then
      if (k_c == lambda) then
         mat4 = Idm+arl_lb*a_par_lb*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,k_c-1))
      else
         mat4 = Idm+arl*a_par*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,k_c-1))
      endif
   endif
   if (k_c < lambda) then
      mat5 = -Idm+arl*a_par*mat2
      call vec_product_left_f(n_c,mat5,VV_c(:,k_c))
   endif
   call calculate_prod_array_f(mat1,RR_ck,mat3); mat6 = mat6+Idm-mat3

   call vec_product_f(mat6,n_c,vec1); SS_c = sum(n_c*vec1)
   FF_c = -(1.d0+coef_res)*relvel_c

   call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,vec_array1)
   call block_tridiagonal_solver_f(lambda-1,DD,UU,LL,VV_c,vec_array2)

   II = 0.d0; tempc = 0.d0; tempd = 0.d0
   do j = 1, lambda-1
      tempc = tempc + sum(VV_c(:,j)*vec_array1(:,j))
      tempd = tempd + sum(VV_c(:,j)*vec_array2(:,j))
   enddo
   do j = 1, lambda-1
      II(:,j+1) = vec_array1(:,j) + tempc/(SS_c-tempd)*vec_array2(:,j)
   enddo
   tempc = 0.d0
   do j = 1, lambda-1; tempc = tempc + sum(VV_c(:,j)*II(:,j+1)); enddo
   II_c = (FF_c-tempc)/SS_c

   do j = 1, lambda
      if (j == lambda) then 
         vxx_add(j,:) = II(:,j+1) - II(:,j)
         call vec_product_f(OO(j,:,:),II(:,j+1)+II(:,j),vec1)
         call vec_product_f(invJJ(j,:,:),vec1,vec2)
         omg_add(j,:) = arl_lb*a_par_lb*vec2
      else 
         vxx_add(j,:) = II(:,j+1) - II(:,j)
         call vec_product_f(OO(j,:,:),II(:,j+1)+II(:,j),vec1)
         call vec_product_f(invJJ(j,:,:),vec1,vec2)
         omg_add(j,:) = arl*a_par*vec2
      endif
   enddo
   vxx_add(j_c,:) = vxx_add(j_c,:) + II_c*n_c
   call vec_product_f(RR_cj,n_c,vec1)
   call vec_product_f(invJJ(j_c,:,:),vec1,vec2)
   omg_add(j_c,:) = omg_add(j_c,:) + II_c*vec2
   vxx_add(k_c,:) = vxx_add(k_c,:) - II_c*n_c
   call vec_product_f(RR_ck,n_c,vec1)
   call vec_product_f(invJJ(k_c,:,:),vec1,vec2)
   omg_add(k_c,:) = omg_add(k_c,:) - II_c*vec2
   
endif
end subroutine in_fibre_collision_solver_f_type


! ─── ex_fibre_collision_solver_f_type ──────────────────────────────────────────────
subroutine ex_fibre_collision_solver_f_type(vxx0,vxx1,ang0,ang1,omg0,omg1, &
     j_c,k_c,r_cj,r_ck,n_c,coef_res,relvel_c,vxx0_add,vxx1_add,omg0_add,omg1_add,fr,fs)
use fibers
implicit none
integer, intent(in) :: j_c, k_c, fr, fs
double precision, intent(in)  :: r_cj(3), r_ck(3), n_c(3), coef_res, relvel_c
double precision, intent(in)  :: vxx0(lambda,3), vxx1(lambda,3)
double precision, intent(in)  :: ang0(lambda,3), ang1(lambda,3)
double precision, intent(in)  :: omg0(lambda,3), omg1(lambda,3)
double precision, intent(out) :: vxx0_add(lambda,3), vxx1_add(lambda,3)
double precision, intent(out) :: omg0_add(lambda,3), omg1_add(lambda,3)

! dim2 = 2*lambda-2 is now a runtime quantity (lambda is a module variable, not parameter)
double precision :: OO0(lambda,3,3), invJJ0(lambda,3,3), JJ0(lambda,3,3)
double precision :: OO1(lambda,3,3), invJJ1(lambda,3,3), JJ1(lambda,3,3)
double precision :: DD(3,3,2*lambda-2), UU(3,3,2*lambda-3), LL(3,3,2*lambda-3)
double precision :: RHS_f(3,2*lambda-2), VV_c(3,2*lambda-2)
double precision :: va1(3,2*lambda-2), va2(3,2*lambda-2)
double precision :: II0(3,lambda+1), II1(3,lambda+1)
double precision :: RR_cj(3,3), RR_ck(3,3)
double precision :: mat1(3,3), mat2(3,3), mat3(3,3), mat4(3,3), mat5(3,3), mat6(3,3)
double precision :: vec1(3), vec2(3)
double precision :: SS_c, FF_c, II_c, tempa, tempb, tempc, tempd
integer :: j, k
! for multibead
double precision :: arl_cr, a_par_cr

if (fiber_types_flag == 1) then     ! multifiber case

   do j = 1, lambda
      call calculate_O_f(ang0(j,:),OO0(j,:,:))
      call calculate_J_f_type(JJ0(j,:,:),ang0(j,:),fr)
      JJ0(j,:,:)=0.5d0*a_par_type(fr)**2.d0*JJ0(j,:,:)
      call INVERSE_MATRIX_ff(JJ0(j,:,:),invJJ0(j,:,:))
      if (j > 1) then
         call calculate_prod_array_f(OO0(j-1,:,:),invJJ0(j-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO0(j-1,:,:),mat2)
         call calculate_prod_array_f(OO0(j,:,:),invJJ0(j,:,:),mat3)
         call calculate_prod_array_f(mat3,OO0(j,:,:),mat4)
         DD(:,:,j-1) = 2.d0*Idm-(arl_type(fr)*a_par_type(fr))**2.d0*(mat2+mat4)
         if (j < lambda) then; UU(:,:,j-1)=-Idm-(arl_type(fr)*a_par_type(fr))**2.d0*mat4; LL(:,:,j-1)=UU(:,:,j-1); endif
         call vec_product_f(OO0(j-1,:,:),omg0(j-1,:),vec1)
         call vec_product_f(OO0(j,:,:),omg0(j,:),vec2)
         RHS_f(:,j-1) = -vxx0(j-1,:)+vxx0(j,:)+arl_type(fr)*a_par_type(fr)*(vec1+vec2)
      endif
   enddo
   do k = 1, lambda
      call calculate_O_f(ang1(k,:),OO1(k,:,:))
      call calculate_J_f_type(JJ1(k,:,:),ang1(k,:),fs,k)
      JJ1(k,:,:)=0.5d0*a_par_type(fs)**2.d0*JJ1(k,:,:)
      call INVERSE_MATRIX_ff(JJ1(k,:,:),invJJ1(k,:,:))
      if (k > 1) then
         call calculate_prod_array_f(OO1(k-1,:,:),invJJ1(k-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO1(k-1,:,:),mat2)
         call calculate_prod_array_f(OO1(k,:,:),invJJ1(k,:,:),mat3)
         call calculate_prod_array_f(mat3,OO1(k,:,:),mat4)
         DD(:,:,k+lambda-2) = 2.d0*Idm-(arl_type(fs)*a_par_type(fs))**2.d0*(mat2+mat4)
         if (k < lambda) then; UU(:,:,k+lambda-2)=-Idm-(arl_type(fs)*a_par_type(fs))**2.d0*mat4; LL(:,:,k+lambda-2)=UU(:,:,k+lambda-2); endif
         call vec_product_f(OO1(k-1,:,:),omg1(k-1,:),vec1)
         call vec_product_f(OO1(k,:,:),omg1(k,:),vec2)
         RHS_f(:,k+lambda-2) = -vxx1(k-1,:)+vxx1(k,:)+arl_type(fs)*a_par_type(fs)*(vec1+vec2)
      endif
   enddo
   UU(:,:,lambda-1) = 0.d0;  LL(:,:,lambda-1) = 0.d0

   VV_c = 0.d0
   call calculate_O_f(r_cj,RR_cj)
   call calculate_prod_array_f(RR_cj,invJJ0(j_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO0(j_c,:,:),mat2)
   if (j_c > 1) then; mat4=-Idm-arl_type(fr)*a_par_type(fr)*mat2; call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1)); endif
   if (j_c < lambda) then; mat5=Idm-arl_type(fr)*a_par_type(fr)*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,j_c)); endif
   call calculate_prod_array_f(mat1,RR_cj,mat3); mat6=Idm-mat3

   call calculate_O_f(r_ck,RR_ck)
   call calculate_prod_array_f(RR_ck,invJJ1(k_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO1(k_c,:,:),mat2)
   if (k_c > 1) then; mat4=Idm+arl_type(fs)*a_par_type(fs)*mat2; call vec_product_left_f(n_c,mat4,VV_c(:,k_c+lambda-2)); endif
   if (k_c < lambda) then; mat5=-Idm+arl_type(fs)*a_par_type(fs)*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,k_c+lambda-1)); endif
   call calculate_prod_array_f(mat1,RR_ck,mat3); mat6=mat6+Idm-mat3

   call vec_product_f(mat6,n_c,vec1); SS_c=sum(n_c*vec1)
   FF_c = -(1.d0+coef_res)*relvel_c

   call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,va1)
   call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,VV_c,va2)

   II0=0.d0; II1=0.d0; tempc=0.d0; tempd=0.d0
   do j=1,2*lambda-2; tempc=tempc+sum(VV_c(:,j)*va1(:,j)); tempd=tempd+sum(VV_c(:,j)*va2(:,j)); enddo
   do j=1,lambda-1; II0(:,j+1)=va1(:,j)+tempc/(SS_c-tempd)*va2(:,j); enddo
   do k=1,lambda-1; II1(:,k+1)=va1(:,k+lambda-1)+tempc/(SS_c-tempd)*va2(:,k+lambda-1); enddo
   tempc=0.d0
   do j=1,lambda-1; tempc=tempc+sum(VV_c(:,j)*II0(:,j+1)); enddo
   do k=1,lambda-1; tempc=tempc+sum(VV_c(:,k+lambda-1)*II1(:,k+1)); enddo
   II_c=(FF_c-tempc)/SS_c

   do j=1,lambda
      vxx0_add(j,:)=II0(:,j+1)-II0(:,j)
      call vec_product_f(OO0(j,:,:),II0(:,j+1)+II0(:,j),vec1)
      call vec_product_f(invJJ0(j,:,:),vec1,vec2); omg0_add(j,:)=arl_type(fr)*a_par_type(fr)*vec2
   enddo
   do k=1,lambda
      vxx1_add(k,:)=II1(:,k+1)-II1(:,k)
      call vec_product_f(OO1(k,:,:),II1(:,k+1)+II1(:,k),vec1)
      call vec_product_f(invJJ1(k,:,:),vec1,vec2); omg1_add(k,:)=arl_type(fs)*a_par_type(fs)*vec2
   enddo
   vxx0_add(j_c,:)=vxx0_add(j_c,:)+II_c*n_c
   call vec_product_f(RR_cj,n_c,vec1); call vec_product_f(invJJ0(j_c,:,:),vec1,vec2)
   omg0_add(j_c,:)=omg0_add(j_c,:)+II_c*vec2
   vxx1_add(k_c,:)=vxx1_add(k_c,:)-II_c*n_c
   call vec_product_f(RR_ck,n_c,vec1); call vec_product_f(invJJ1(k_c,:,:),vec1,vec2)
   omg1_add(k_c,:)=omg1_add(k_c,:)-II_c*vec2

elseif(bead_types_flag == 1) then      ! multibead case

   do j = 1, lambda
      !to assign right bead dimensions
      if (j == lambda) then 
         a_par_cr = a_par_lb; arl_cr = arl_lb
      else
         a_par_cr = a_par; arl_cr = arl
      endif

      call calculate_O_f(ang0(j,:),OO0(j,:,:))
      call calculate_J_f_type(JJ0(j,:,:),ang0(j,:),fr,j)
      JJ0(j,:,:)=0.5d0*a_par_cr**2.d0*JJ0(j,:,:)
      call INVERSE_MATRIX_ff(JJ0(j,:,:),invJJ0(j,:,:))
      if (j > 1) then
         call calculate_prod_array_f(OO0(j-1,:,:),invJJ0(j-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO0(j-1,:,:),mat2)
         !Also in this case mat2 refered to j-1-th term
         call calculate_prod_array_f(OO0(j,:,:),invJJ0(j,:,:),mat3)
         call calculate_prod_array_f(mat3,OO0(j,:,:),mat4)
         !mat4 associated to j-th term
         if(j == lambda) then
            DD(:,:,j-1) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2)
                                  -(arl_cr*a_par_cr)**2.d0*(mat4)
         else
            DD(:,:,j-1) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2+mat4)
         endif
         if (j < lambda) then; UU(:,:,j-1)=-Idm-(arl*a_par)**2.d0*mat4; LL(:,:,j-1)=UU(:,:,j-1); endif
         call vec_product_f(OO0(j-1,:,:),omg0(j-1,:),vec1)  !refered to j-1-th
         call vec_product_f(OO0(j,:,:),omg0(j,:),vec2)
         if(j == lambda) then
            RHS_f(:,j-1) = -vxx0(j-1,:)+vxx0(j,:)+arl*a_par*(vec1)
                                                 +arl_cr*a_par_cr*(vec2)
         else
            RHS_f(:,j-1) = -vxx0(j-1,:)+vxx0(j,:)+arl*a_par*(vec1+vec2)
         endif
      endif
   enddo
   do k = 1, lambda

      if (k == lambda) then 
         a_par_cr = a_par_lb; arl_cr = arl_lb
      else
         a_par_cr = a_par; arl_cr = arl
      endif

      call calculate_O_f(ang1(k,:),OO1(k,:,:))
      call calculate_J_f_type(JJ1(k,:,:),ang1(k,:),fs)
      JJ1(k,:,:)=0.5d0*a_par_cr**2.d0*JJ1(k,:,:)
      call INVERSE_MATRIX_ff(JJ1(k,:,:),invJJ1(k,:,:))
      if (k > 1) then
         call calculate_prod_array_f(OO1(k-1,:,:),invJJ1(k-1,:,:),mat1)
         call calculate_prod_array_f(mat1,OO1(k-1,:,:),mat2)      !j-1-th
         call calculate_prod_array_f(OO1(k,:,:),invJJ1(k,:,:),mat3)
         call calculate_prod_array_f(mat3,OO1(k,:,:),mat4)
         if (k == lambda) then
            DD(:,:,k+lambda-2) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2)
                                         -(arl_cr*a_par_cr)**2.d0*(mat4)
         else
            DD(:,:,k+lambda-2) = 2.d0*Idm-(arl*a_par)**2.d0*(mat2+mat4)
         endif
         if (k < lambda) then; UU(:,:,k+lambda-2)=-Idm-(arl_cr*a_par_cr)**2.d0*mat4; LL(:,:,k+lambda-2)=UU(:,:,k+lambda-2); endif
         call vec_product_f(OO1(k-1,:,:),omg1(k-1,:),vec1)     !j-1-th
         call vec_product_f(OO1(k,:,:),omg1(k,:),vec2)
         if (k == lambda) then 
            RHS_f(:,k+lambda-2) = -vxx1(k-1,:)+vxx1(k,:)+arl*a_par*(vec1)
                                                        +arl_cr*a_par_cr*(vec2)
         else
            RHS_f(:,k+lambda-2) = -vxx1(k-1,:)+vxx1(k,:)+arl*a_par*(vec1+vec2)
         endif
      endif
   enddo
   UU(:,:,lambda-1) = 0.d0;  LL(:,:,lambda-1) = 0.d0

   VV_c = 0.d0
   call calculate_O_f(r_cj,RR_cj)
   call calculate_prod_array_f(RR_cj,invJJ0(j_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO0(j_c,:,:),mat2)
   if (j_c > 1) then
      if (j_c == lambda) then  
         mat4=-Idm-arl_lb*a_par_lb*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1))
      else 
         mat4=-Idm-arl*a_par*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,j_c-1))
      endif
   endif
   !NB following line modified coherently knowing that asymmetry just on the last bead NB 
   if (j_c < lambda) then; mat5=Idm-arl*a_par*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,j_c)); endif
   call calculate_prod_array_f(mat1,RR_cj,mat3); mat6=Idm-mat3

   call calculate_O_f(r_ck,RR_ck)
   call calculate_prod_array_f(RR_ck,invJJ1(k_c,:,:),mat1)
   call calculate_prod_array_f(mat1,OO1(k_c,:,:),mat2)
   if (k_c > 1) then
      if (k_c == lambda) then 
         mat4=Idm+arl_lb*a_par_lb*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,k_c+lambda-2))
      else
         mat4=Idm+arl*a_par*mat2
         call vec_product_left_f(n_c,mat4,VV_c(:,k_c+lambda-2))
      endif
   endif
   !NB following line modified coherently knowing that asymmetry just on the last bead NB
   if (k_c < lambda) then; mat5=-Idm+arl*a_par*mat2; call vec_product_left_f(n_c,mat5,VV_c(:,k_c+lambda-1)); endif
   call calculate_prod_array_f(mat1,RR_ck,mat3); mat6=mat6+Idm-mat3

   call vec_product_f(mat6,n_c,vec1); SS_c=sum(n_c*vec1)
   FF_c = -(1.d0+coef_res)*relvel_c

   call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,RHS_f-FF_c/SS_c*VV_c,va1)
   call block_tridiagonal_solver_f(2*lambda-2,DD,UU,LL,VV_c,va2)

   II0=0.d0; II1=0.d0; tempc=0.d0; tempd=0.d0
   do j=1,2*lambda-2; tempc=tempc+sum(VV_c(:,j)*va1(:,j)); tempd=tempd+sum(VV_c(:,j)*va2(:,j)); enddo
   do j=1,lambda-1; II0(:,j+1)=va1(:,j)+tempc/(SS_c-tempd)*va2(:,j); enddo
   do k=1,lambda-1; II1(:,k+1)=va1(:,k+lambda-1)+tempc/(SS_c-tempd)*va2(:,k+lambda-1); enddo
   tempc=0.d0
   do j=1,lambda-1; tempc=tempc+sum(VV_c(:,j)*II0(:,j+1)); enddo
   do k=1,lambda-1; tempc=tempc+sum(VV_c(:,k+lambda-1)*II1(:,k+1)); enddo
   II_c=(FF_c-tempc)/SS_c

   do j=1,lambda
      if (j == lambda) then 
         vxx0_add(j,:)=II0(:,j+1)-II0(:,j)
         call vec_product_f(OO0(j,:,:),II0(:,j+1)+II0(:,j),vec1)
         call vec_product_f(invJJ0(j,:,:),vec1,vec2); omg0_add(j,:)=arl_lb*a_par_lb*vec2
      else 
         vxx0_add(j,:)=II0(:,j+1)-II0(:,j)
         call vec_product_f(OO0(j,:,:),II0(:,j+1)+II0(:,j),vec1)
         call vec_product_f(invJJ0(j,:,:),vec1,vec2); omg0_add(j,:)=arl*a_par*vec2
      endif
   enddo
   do k=1,lambda
      if (k == lambda) then 
         vxx1_add(k,:)=II1(:,k+1)-II1(:,k)
         call vec_product_f(OO1(k,:,:),II1(:,k+1)+II1(:,k),vec1)
         call vec_product_f(invJJ1(k,:,:),vec1,vec2); omg1_add(k,:)=arl_lb*a_par_lb*vec2
      else
         vxx1_add(k,:)=II1(:,k+1)-II1(:,k)
         call vec_product_f(OO1(k,:,:),II1(:,k+1)+II1(:,k),vec1)
         call vec_product_f(invJJ1(k,:,:),vec1,vec2); omg1_add(k,:)=arl*a_par*vec2
      endif
   enddo
   vxx0_add(j_c,:)=vxx0_add(j_c,:)+II_c*n_c
   call vec_product_f(RR_cj,n_c,vec1); call vec_product_f(invJJ0(j_c,:,:),vec1,vec2)
   omg0_add(j_c,:)=omg0_add(j_c,:)+II_c*vec2
   vxx1_add(k_c,:)=vxx1_add(k_c,:)-II_c*n_c
   call vec_product_f(RR_ck,n_c,vec1); call vec_product_f(invJJ1(k_c,:,:),vec1,vec2)
   omg1_add(k_c,:)=omg1_add(k_c,:)-II_c*vec2
endif

end subroutine ex_fibre_collision_solver_f_type
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm


! ─────────────────────────── block_tridiagonal_solver_f ──────────────────────
subroutine block_tridiagonal_solver_f(dim,DD,UU,LL,FF,XX)
implicit none
integer, intent(in) :: dim
double precision, intent(in)  :: DD(3,3,dim), UU(3,3,dim-1), LL(3,3,dim-1)
double precision, intent(in)  :: FF(3,dim)
double precision, intent(out) :: XX(3,dim)
double precision :: D0(3,3,dim), invD0(3,3,dim), F0(3,dim), mat_tmp(3,3), vec_tmp(3)
integer :: j

do j = dim, 1, -1
   if (j == dim) then
      F0(:,j) = FF(:,j);  D0(:,:,j) = DD(:,:,j)
      call invertmatrix_constraint_f(D0(:,:,j), invD0(:,:,j))
   else
      call vec_product_f(invD0(:,:,j+1), F0(:,j+1), vec_tmp)
      call vec_product_f(UU(:,:,j), vec_tmp, F0(:,j))
      F0(:,j) = FF(:,j) - F0(:,j)
      call calculate_prod_array_f(invD0(:,:,j+1), LL(:,:,j), mat_tmp)
      call calculate_prod_array_f(UU(:,:,j), mat_tmp, D0(:,:,j))
      D0(:,:,j) = DD(:,:,j) - D0(:,:,j)
      call invertmatrix_constraint_f(D0(:,:,j), invD0(:,:,j))
   endif
enddo
do j = 1, dim
   if (j == 1) then
      call vec_product_f(invD0(:,:,j), F0(:,j), XX(:,j))
   else
      call vec_product_f(LL(:,:,j-1), XX(:,j-1), vec_tmp)
      F0(:,j) = F0(:,j) - vec_tmp
      call vec_product_f(invD0(:,:,j), F0(:,j), XX(:,j))
   endif
enddo
end subroutine block_tridiagonal_solver_f


! ============================================================================
!  Diagnostic: Jeffery-orbit analytics for a single fiber
! ============================================================================
subroutine jeffery_analytics_f(xxx, step, time)
use param
use fibers
implicit none
double precision, intent(in) :: xxx(lambda,3), time
integer, intent(in) :: step
integer :: i
double precision :: cm(3), am(3), nxt(3)
double precision :: x1, x2, deltax, deltay, deltaz, eed
logical :: file_exists

cm = 0.d0;  am = 0.d0;  eed = 0.d0
cm = xxx(1,:);  nxt = xxx(1,:)
do i = 2, lambda
   nxt = nxt + a_par*arl*(angles_f(i-1,:)+angles_f(i,:))
   cm  = cm  + nxt
enddo
cm = cm / lambda

do i = 1, lambda; am = am + angles_f(i,:); enddo
am = am / lambda

x1 = xxx(lambda,1) + a_par*arl*angles_f(lambda,1)
x2 = xxx(1,1)      - a_par*arl*angles_f(1,1)
deltax = x1-x2
if (deltax >  0.5d0*lx) deltax = deltax - lx
if (deltax <= -0.5d0*lx) deltax = deltax + lx
x1 = xxx(lambda,2)+a_par*arl*angles_f(lambda,2)
x2 = xxx(1,2)     -a_par*arl*angles_f(1,2)
deltay = x1-x2
if (deltay >  0.5d0*lx) deltay = deltay - lx
if (deltay <= -0.5d0*lx) deltay = deltay + lx
deltaz = xxx(lambda,3)+a_par*arl*angles_f(lambda,3) - (xxx(1,3)-a_par*arl*angles_f(1,3))
eed = sqrt(deltax**2.d0+deltay**2.d0+deltaz**2.d0)

if (step == nstart) then
   INQUIRE(FILE='jeffery.dat', EXIST=file_exists)
   if (file_exists) then
      open(unit=66,file='jeffery.dat',status='old'); close(66,status='delete')
   endif
   open(66,file='jeffery.dat',form='formatted',status='new')
   write(66,'(a,i4)') '# Jeffery orbit, fiber ar=', int(lambda*arl)
   write(66,'(8(a12,2x))') 'time','xcm','ycm','zcm','px','py','pz','eed'
   write(66,'(8(es12.5,2x))') time, cm(:), am(:), eed
   close(66)
else
   open(66,file='jeffery.dat',form='formatted',status='old',position='append')
   write(66,'(8(es12.5,2x))') time, cm(:), am(:), eed
   close(66)
endif
end subroutine jeffery_analytics_f


! ============================================================================
!  Elementary linear-algebra helpers (all names suffixed _f to avoid conflicts)
! ============================================================================

subroutine vec_product_f(A, v, res)
implicit none
double precision, intent(in)  :: A(3,3), v(3)
double precision, intent(out) :: res(3)
res(1) = A(1,1)*v(1)+A(1,2)*v(2)+A(1,3)*v(3)
res(2) = A(2,1)*v(1)+A(2,2)*v(2)+A(2,3)*v(3)
res(3) = A(3,1)*v(1)+A(3,2)*v(2)+A(3,3)*v(3)
end subroutine vec_product_f

subroutine vec_product_left_f(v, A, res)
implicit none
double precision, intent(in)  :: v(3), A(3,3)
double precision, intent(out) :: res(3)
res(1) = v(1)*A(1,1)+v(2)*A(2,1)+v(3)*A(3,1)
res(2) = v(1)*A(1,2)+v(2)*A(2,2)+v(3)*A(3,2)
res(3) = v(1)*A(1,3)+v(2)*A(2,3)+v(3)*A(3,3)
end subroutine vec_product_left_f

subroutine calculate_prod_array_f(a1, a2, res)
implicit none
double precision, intent(in)  :: a1(3,3), a2(3,3)
double precision, intent(out) :: res(3,3)
res(1,1)=a1(1,1)*a2(1,1)+a1(1,2)*a2(2,1)+a1(1,3)*a2(3,1)
res(2,1)=a1(2,1)*a2(1,1)+a1(2,2)*a2(2,1)+a1(2,3)*a2(3,1)
res(3,1)=a1(3,1)*a2(1,1)+a1(3,2)*a2(2,1)+a1(3,3)*a2(3,1)
res(1,2)=a1(1,1)*a2(1,2)+a1(1,2)*a2(2,2)+a1(1,3)*a2(3,2)
res(2,2)=a1(2,1)*a2(1,2)+a1(2,2)*a2(2,2)+a1(2,3)*a2(3,2)
res(3,2)=a1(3,1)*a2(1,2)+a1(3,2)*a2(2,2)+a1(3,3)*a2(3,2)
res(1,3)=a1(1,1)*a2(1,3)+a1(1,2)*a2(2,3)+a1(1,3)*a2(3,3)
res(2,3)=a1(2,1)*a2(1,3)+a1(2,2)*a2(2,3)+a1(2,3)*a2(3,3)
res(3,3)=a1(3,1)*a2(1,3)+a1(3,2)*a2(2,3)+a1(3,3)*a2(3,3)
end subroutine calculate_prod_array_f

subroutine calculate_O_f(o, oo)
implicit none
double precision, intent(in)  :: o(3)
double precision, intent(out) :: oo(3,3)
oo(1,1)= 0.d0;  oo(1,2)=-o(3);   oo(1,3)= o(2)
oo(2,1)= o(3);  oo(2,2)= 0.d0;   oo(2,3)=-o(1)
oo(3,1)=-o(2);  oo(3,2)= o(1);   oo(3,3)= 0.d0
end subroutine calculate_O_f

subroutine calculate_ooT_f(A, o)
implicit none
double precision, intent(in)  :: o(3)
double precision, intent(out) :: A(3,3)
A(1,1)=o(1)*o(1); A(1,2)=o(1)*o(2); A(1,3)=o(1)*o(3)
A(2,1)=o(2)*o(1); A(2,2)=o(2)*o(2); A(2,3)=o(2)*o(3)
A(3,1)=o(3)*o(1); A(3,2)=o(3)*o(2); A(3,3)=o(3)*o(3)
end subroutine calculate_ooT_f

subroutine calculate_J_f(JJ, o)
use fibers
implicit none
double precision, intent(in)  :: o(3)
double precision, intent(out) :: JJ(3,3)
double precision :: OZZ(3,3)
call calculate_ooT_f(OZZ, o)
JJ = (4.d0*arl**2.d0+3.d0)/6.d0*(Idm-OZZ) + OZZ
end subroutine calculate_J_f

subroutine calculate_Jdot_f(Jdot, oo, omegap, o)
! Compute d(J)/dt and fill oo = O(o) as a side effect (used in calculate_RS_f)
use fibers
implicit none
double precision, intent(out)   :: Jdot(3,3)
double precision, intent(inout) :: oo(3,3)   ! will be overwritten with O(o)
double precision, intent(in)    :: omegap(3), o(3)
double precision :: res(3), ars(3,3)
integer :: j, k
call calculate_O_f(o, oo)   ! fills oo = O(o)  ← the key side-effect
call vec_product_f(oo, omegap, res)
do j=1,3; do k=1,3; ars(j,k)=res(j)*o(k); enddo; enddo
do j=1,3; do k=1,3; Jdot(j,k)=(4.d0*arl**2.d0-3.d0)/6.d0*(ars(j,k)+ars(k,j)); enddo; enddo
end subroutine calculate_Jdot_f

subroutine calculate_C_f(CC, o)
use fibers
implicit none
double precision, intent(in)  :: o(3)
double precision, intent(out) :: CC(3,3)
double precision :: OZZ(3,3)
call calculate_ooT_f(OZZ, o)
CC = 8.d0*mu_f*arl**2.d0/(dens_part*a_par**2.d0)*(Y_c*Idm - (X_c-Y_c)*OZZ)
end subroutine calculate_C_f

subroutine calculate_H_f(HH, E, o)
use fibers
implicit none
double precision, intent(in)  :: E(3,3), o(3)
double precision, intent(out) :: HH(3)
! HH_i = -Y_h * eps_ijk (E_jl o_l) o_k  (see Delmotte 2015 eq. A15)
HH(1) = E(2,1)*o(1)*o(3)+E(2,2)*o(2)*o(3)+E(2,3)*o(3)*o(3) &
       -E(3,1)*o(1)*o(2)-E(3,2)*o(2)*o(2)-E(3,3)*o(3)*o(2)
HH(2) =-E(1,1)*o(1)*o(3)-E(1,2)*o(2)*o(3)-E(1,3)*o(3)*o(3) &
        +E(3,1)*o(1)*o(1)+E(3,2)*o(2)*o(1)+E(3,3)*o(3)*o(1)
HH(3) = E(1,1)*o(1)*o(2)+E(1,2)*o(2)*o(2)+E(1,3)*o(3)*o(2) &
       -E(2,1)*o(1)*o(1)-E(2,2)*o(2)*o(1)-E(2,3)*o(3)*o(1)
HH = -8.d0*mu_f*arl**2.d0/(dens_part*a_par**2.d0)*Y_h*HH
end subroutine calculate_H_f




! mmmmmm multifiber modifications mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

subroutine calculate_J_f_type(JJ, o, f,ind)
use fibers
implicit none
double precision, intent(in)  :: o(3)
integer, intent(in) :: f,ind
double precision, intent(out) :: JJ(3,3)
double precision :: OZZ(3,3)
call calculate_ooT_f(OZZ, o)

!initializzation 
JJ = 0.d0

if (fiber_types_flag == 1) then
   JJ = (4.d0*arl_type(f)**2.d0+3.d0)/6.d0*(Idm-OZZ) + OZZ

elseif (bead_types_flag == 1)
   if (ind == lambda) then
      JJ = (4.d0*arl_lb**2.d0+3.d0)/6.d0*(Idm-OZZ) + OZZ
   else 
      JJ = (4.d0*arl**2.d0+3.d0)/6.d0*(Idm-OZZ) + OZZ
   endif
endif

end subroutine calculate_J_f_type

subroutine calculate_Jdot_f_type(Jdot, oo, omegap, o, f, ind)
! Compute d(J)/dt and fill oo = O(o) as a side effect (used in calculate_RS_f)
use fibers
implicit none
double precision, intent(out)   :: Jdot(3,3)
double precision, intent(inout) :: oo(3,3)   ! will be overwritten with O(o)
double precision, intent(in)    :: omegap(3), o(3)
integer, intent(in) :: f, ind
double precision :: res(3), ars(3,3)
integer :: j, k
call calculate_O_f(o, oo)   ! fills oo = O(o)  ← the key side-effect
call vec_product_f(oo, omegap, res)
do j=1,3; do k=1,3; ars(j,k)=res(j)*o(k); enddo; enddo

!initializzation 
Jdot = 0.d0

if (fiber_types_flag == 1) then
   do j=1,3; do k=1,3; Jdot(j,k)=(4.d0*arl_type(f)**2.d0-3.d0)/6.d0*(ars(j,k)+ars(k,j)); enddo; enddo

elseif(bead_types_flag == 1) then
   if ( ind == lambda) then
      do j=1,3; do k=1,3; Jdot(j,k)=(4.d0*arl_lb**2.d0-3.d0)/6.d0*(ars(j,k)+ars(k,j)); enddo; enddo
   else
      do j=1,3; do k=1,3; Jdot(j,k)=(4.d0*arl**2.d0-3.d0)/6.d0*(ars(j,k)+ars(k,j)); enddo; enddo
   endif
endif

end subroutine calculate_Jdot_f_type

subroutine calculate_C_f_type(CC, o, f, ind)
use fibers
implicit none
double precision, intent(in)  :: o(3)
integer, intent(in) :: f, ind
double precision, intent(out) :: CC(3,3)
double precision :: OZZ(3,3)
call calculate_ooT_f(OZZ, o)

!initializzation 
CC = 0.d0

if (fiber_types_flag == 1) then
   CC = 8.d0*mu_f*arl_type(f)**2.d0/(dens_part_type(f)*a_par_type(f)**2.d0) &
            *(Y_c_type(f)*Idm - (X_c_type(f)-Y_c_type(f))*OZZ)
elseif(bead_types_flag == 1) then
   if (ind == lambda) then
      CC = 8.d0*mu_f*arl_lb**2.d0/(dens_part*a_par_lb**2.d0) &
            *(Y_c_lb*Idm - (X_c_lb-Y_c_lb)*OZZ)
   else
      CC = 8.d0*mu_f*arl**2.d0/(dens_part*a_par**2.d0)*(Y_c*Idm - (X_c-Y_c)*OZZ)
   endif
endif

end subroutine calculate_C_f_type

subroutine calculate_H_f_type(HH, E, o, f, ind)
use fibers
implicit none
double precision, intent(in)  :: E(3,3), o(3)
integer, intent(in) :: f, ind
double precision, intent(out) :: HH(3)
! HH_i = -Y_h * eps_ijk (E_jl o_l) o_k  (see Delmotte 2015 eq. A15)
HH(1) = E(2,1)*o(1)*o(3)+E(2,2)*o(2)*o(3)+E(2,3)*o(3)*o(3) &
       -E(3,1)*o(1)*o(2)-E(3,2)*o(2)*o(2)-E(3,3)*o(3)*o(2)
HH(2) =-E(1,1)*o(1)*o(3)-E(1,2)*o(2)*o(3)-E(1,3)*o(3)*o(3) &
        +E(3,1)*o(1)*o(1)+E(3,2)*o(2)*o(1)+E(3,3)*o(3)*o(1)
HH(3) = E(1,1)*o(1)*o(2)+E(1,2)*o(2)*o(2)+E(1,3)*o(3)*o(2) &
       -E(2,1)*o(1)*o(1)-E(2,2)*o(2)*o(1)-E(2,3)*o(3)*o(1)

if (fiber_types_flag == 1) then
   HH = -8.d0*mu_f*arl_type(f)**2.d0/(dens_part_type(f)*a_par_type(f)**2.d0) &
         *Y_h_type(f)*HH
elseif (bead_types_flag == 1) then 
   if (ind == lambda) then
      HH = -8.d0*mu_f*arl_lb**2.d0/(dens_part*a_par_lb**2.d0)*Y_h_lb*HH
   else
      HH = -8.d0*mu_f*arl**2.d0/(dens_part*a_par**2.d0)*Y_h*HH
   endif
endif

end subroutine calculate_H_f_type

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm



subroutine INVERSE_MATRIX_ff(in, ou)
implicit none
double precision, intent(in)  :: in(3,3)
double precision, intent(out) :: ou(3,3)
double precision :: det
det = in(1,1)*(in(2,2)*in(3,3)-in(2,3)*in(3,2)) &
    - in(1,2)*(in(2,1)*in(3,3)-in(3,1)*in(2,3)) &
    + in(1,3)*(in(2,1)*in(3,2)-in(3,1)*in(2,2))
if (det /= det .or. det == 0.d0) then; ou = 0.d0; return; endif
ou(1,1) = (in(2,2)*in(3,3)-in(2,3)*in(3,2))/det
ou(1,2) = (in(1,3)*in(3,2)-in(1,2)*in(3,3))/det
ou(1,3) = (in(1,2)*in(2,3)-in(1,3)*in(2,2))/det
ou(2,1) = (in(3,1)*in(2,3)-in(2,1)*in(3,3))/det
ou(2,2) = (in(1,1)*in(3,3)-in(1,3)*in(3,1))/det
ou(2,3) = (in(1,3)*in(2,1)-in(2,3)*in(1,1))/det
ou(3,1) = (in(2,1)*in(3,2)-in(3,1)*in(2,2))/det
ou(3,2) = (in(3,1)*in(1,2)-in(3,2)*in(1,1))/det
ou(3,3) = (in(1,1)*in(2,2)-in(1,2)*in(2,1))/det
end subroutine INVERSE_MATRIX_ff

subroutine invertmatrix_constraint_f(in, ou)
implicit none
double precision, intent(in)  :: in(3,3)
double precision, intent(out) :: ou(3,3)
double precision :: det
det = in(1,1)*(in(2,2)*in(3,3)-in(2,3)*in(3,2)) &
    - in(1,2)*(in(2,1)*in(3,3)-in(3,1)*in(2,3)) &
    + in(1,3)*(in(2,1)*in(3,2)-in(3,1)*in(2,2))
if (det /= det .or. det == 0.d0) then; ou = 0.d0; return; endif
ou(1,1) = (in(2,2)*in(3,3)-in(2,3)*in(3,2))/det
ou(1,2) = (in(1,3)*in(3,2)-in(1,2)*in(3,3))/det
ou(1,3) = (in(1,2)*in(2,3)-in(1,3)*in(2,2))/det
ou(2,1) = (in(3,1)*in(2,3)-in(2,1)*in(3,3))/det
ou(2,2) = (in(1,1)*in(3,3)-in(1,3)*in(3,1))/det
ou(2,3) = (in(1,3)*in(2,1)-in(2,3)*in(1,1))/det
ou(3,1) = (in(2,1)*in(3,2)-in(3,1)*in(2,2))/det
ou(3,2) = (in(3,1)*in(1,2)-in(3,2)*in(1,1))/det
ou(3,3) = (in(1,1)*in(2,2)-in(1,2)*in(2,1))/det
end subroutine invertmatrix_constraint_f


!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!!Author: L. Terlicher
!!Set 2026
!!Copyright: Multiphase-flow lab, University of Udine
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

subroutine friction_calculation(QQ,SS,RR,oo,fttrq_add,Ktot,idfiber)

    use param
    use fibers

implicit none


double precision, dimension(3,3,lambda) :: oo
double precision, dimension(3,3,lambda) :: QQ, SS
double precision, dimension(3,lambda) :: KK, RR, Ktot , X_int1, X_int2, XX
double precision, dimension(3,lambda) :: ft_add, fttrq_add, fttrqo, FT, FTTRQ
double precision,dimension(3, lambda) :: Ftot
integer,intent(in) :: idfiber
integer :: ind,f

f = min((idfiber-1)/type_idx + 1, n_fiber_types)

! same ind
!-------------  
! product of forces with Q matrix 
call vec_product_f(QQ(:,:,ind),ft_add(:,ind),FT(:,ind))
! product of torques with S matrix 
call vec_product_f(SS(:,:,ind), fttrq_add(:,ind),fttrqo(:,ind))
! now we want to express it also in function of the a and lambda O
call vec_product_f(oo(:,:,ind),fttrqo(:,ind), FTTRQ(:,ind) )

!next ind 
call vec_product_f(QQ(:,:,ind+1),ft_add(:,ind+1),FT(:,ind+1))
call vec_product_f(SS(:,:,ind+1), fttrq_add(:,ind+1),fttrqo(:,ind+1))
call vec_product_f(oo(:,:,ind+1),fttrqo(:,ind+1), FTTRQ(:,ind+1))


!multiplying by the a and lambda
if (fiber_types_flag == 1) then
   FTTRQ(:,ind) = FTTRQ(:,ind) * arl_type(f) * a_par_type(f)
   FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl_type(f) * a_par_type(f)
else
   FTTRQ(:,ind) = FTTRQ(:,ind) * arl * a_par
   FTTRQ(:,ind+1) = FTTRQ(:,ind+1) * arl * a_par
endif


! summing those thoghether 
Ftot(1,ind) = -FT(1,ind+1) + FT(1,ind) - FTTRQ(1,ind+1) - FTTRQ(1,ind)
Ftot(2,ind) = -FT(2,ind+1) + FT(2,ind) - FTTRQ(2,ind+1) - FTTRQ(2,ind)
Ftot(3,ind) = -FT(3,ind+1) + FT(3,ind) - FTTRQ(3,ind+1) - FTTRQ(3,ind)

! adding to k_tot

Ktot(1,ind) = Ktot(1,ind) + Ftot(1,ind)
Ktot(2,ind) = Ktot(2,ind) + Ftot(2,ind) 
Ktot(3,ind) = Ktot(3,ind) + Ftot(3,ind)


end subroutine friction_calculation
