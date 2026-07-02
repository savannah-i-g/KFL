!> k26astro_ode_iface — ISO_C_BINDING surface for libk26astro_ode.
!>
!> Exposes DLSODA (from M_odepack) as a plain-C symbol. The C-side
!> wrapper in k26astro_ode.c sets thread-local state (rhs_fn, user,
!> n) before invoking the bind(C) entry; the Fortran rhs_shim calls
!> back into a C trampoline that reads that TLS.
!>
!> Workspace sized dynamically based on `n` using F90 allocatable
!> arrays. DLSODA workspace formula (full Jacobian, JT=2):
!>   LRW = 22 + n * MAX(16, n + 9)
!>   LIW = 20 + n
!>
!> Determinism: -O2 -ffp-contract=off; no I/O in the integration
!> path; LSODA's stiff/non-stiff switching is determined by error
!> estimates only (no thread scheduling).

module k26astro_ode_iface
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use M_odepack, only: dlsoda
    implicit none

    private

    interface
        subroutine k26astro_ode_tls_rhs_trampoline(n, t, y, ydot) &
            bind(C, name="k26astro_ode_tls_rhs_trampoline")
            use, intrinsic :: iso_c_binding, only: c_int, c_double
            integer(c_int), value, intent(in) :: n
            real(c_double), value, intent(in) :: t
            real(c_double), intent(in)        :: y(*)
            real(c_double), intent(out)       :: ydot(*)
        end subroutine k26astro_ode_tls_rhs_trampoline
    end interface

contains

    !> bind(C) entry for DLSODA. C-side wrapper sets TLS before
    !> calling; F-side rhs_shim reads TLS via the C trampoline.
    subroutine k26astro_ode_lsoda_call(n, y, t0, tf, rtol, atol, &
                                        out_y, out_istate) &
        bind(C, name="k26astro_ode_lsoda_call")
        integer(c_int), value, intent(in) :: n
        real(c_double), intent(in)        :: y(*)
        real(c_double), value, intent(in) :: t0, tf, rtol, atol
        real(c_double), intent(out)       :: out_y(*)
        integer(c_int), intent(out)       :: out_istate

        ! Workspace (dynamically sized; freed at subroutine exit).
        real(c_double), allocatable :: rwork(:)
        integer, allocatable        :: iwork(:)
        real(c_double)              :: t_loc, tout_loc
        real(c_double)              :: rtol_arr(1), atol_arr(1)
        integer                     :: neq_arr(1)
        integer                     :: itol, itask, istate, iopt
        integer                     :: lrw, liw, jt, i

        ! Workspace formula from DLSODA manual (full Jacobian path).
        lrw = 22 + n * max(16, n + 9)
        liw = 20 + n

        allocate(rwork(lrw))
        allocate(iwork(liw))
        rwork(:) = 0.0_c_double
        iwork(:) = 0

        ! Copy initial state into the output buffer; DLSODA reads/writes
        ! `y` in place, so we use out_y as the working state vector.
        do i = 1, n
            out_y(i) = y(i)
        end do

        t_loc    = t0
        tout_loc = tf
        rtol_arr(1) = rtol
        atol_arr(1) = atol
        neq_arr(1)  = n
        itol   = 1   ! scalar rtol + scalar atol
        itask  = 1   ! normal integration to tout
        istate = 1   ! first call
        iopt   = 1   ! optional inputs in iwork/rwork
        jt     = 2   ! internally generated full Jacobian

        ! Optional inputs: raise MXSTEP from the LSODA default (500)
        ! to 50000 so stiff benchmarks (Van der Pol μ=1000, Robertson
        ! integrated past 1e4) complete without spurious -1 returns.
        ! The cost is purely time-bounded: LSODA still uses adaptive
        ! step control; MXSTEP is only a wall on runaway integrations.
        iwork(6) = 50000   ! MXSTEP

        call dlsoda(rhs_shim, neq_arr, out_y, t_loc, tout_loc, itol, &
                    rtol_arr, atol_arr, itask, istate, iopt, &
                    rwork, lrw, iwork, liw, jac_dummy, jt)

        out_istate = istate

        deallocate(rwork)
        deallocate(iwork)
    end subroutine k26astro_ode_lsoda_call

    !> RHS shim — matches DLSODA's `external f` signature
    !> (subroutine(neq, t, y, ydot)). Forwards via the C trampoline
    !> that reads the thread-local (rhs_fn, user) pair set by the C
    !> wrapper.
    subroutine rhs_shim(neq, t, y, ydot)
        integer       :: neq(*)
        real(c_double) :: t
        real(c_double) :: y(*)
        real(c_double) :: ydot(*)
        call k26astro_ode_tls_rhs_trampoline(neq(1), t, y, ydot)
    end subroutine rhs_shim

    !> Dummy Jacobian — DLSODA requires the symbol but JT=2
    !> generates the Jacobian internally; this subroutine is never
    !> called.
    subroutine jac_dummy(neq, t, y, ml, mu, pd, nrowpd)
        integer       :: neq(*)
        real(c_double) :: t
        real(c_double) :: y(*)
        integer       :: ml, mu, nrowpd
        real(c_double) :: pd(nrowpd, *)
        ! Intentionally empty.
    end subroutine jac_dummy

end module k26astro_ode_iface
