!> k26astro_quad_iface — ISO_C_BINDING surface for libk26astro_quad.
!>
!> Exposes the QUADPACK DQAGS and DQAGI entries as plain-C symbols
!> (no Fortran name mangling). The C-side wrapper in k26astro_quad.c
!> sets a thread-local (fn, user) pair before invoking the bind(C)
!> entry; the Fortran integrand shim calls back through a C trampoline
!> that reads that thread-local state.
!>
!> Determinism: -O2 -ffp-contract=off; no I/O; workspace arrays are
!> stack-allocated to keep allocator state out of the integration path.

module k26astro_quad_iface
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use quadpack_double, only: dqags, dqagi
    implicit none

    private

    ! Workspace sizing — must match include/k26astro_quad/quad_consts.h
    integer, parameter :: K26_QUAD_LIMIT = 500
    integer, parameter :: K26_QUAD_LENW  = 4 * K26_QUAD_LIMIT

    ! C-side trampoline used by the Fortran integrand shim. Declared
    ! once here, accessed from the contains-block shims below.
    interface
        function k26astro_quad_tls_trampoline(x) bind(C, name="k26astro_quad_tls_trampoline") result(y)
            use, intrinsic :: iso_c_binding, only: c_double
            real(c_double), value, intent(in) :: x
            real(c_double) :: y
        end function k26astro_quad_tls_trampoline
    end interface

contains

    !> DQAGS bind(C) entry. Called by the C-side k26astro_quad_dqags
    !> after the thread-local (fn, user) pair has been set.
    subroutine k26astro_quad_dqags_call(a, b, epsabs, epsrel, &
                                         out_result, out_abserr, &
                                         out_neval, out_ier) &
        bind(C, name="k26astro_quad_dqags_call")
        real(c_double), value, intent(in) :: a, b, epsabs, epsrel
        real(c_double), intent(out)       :: out_result, out_abserr
        integer(c_int), intent(out)       :: out_neval, out_ier

        real(c_double) :: result_loc, abserr_loc
        integer        :: neval_loc, ier_loc, last_loc
        integer        :: iwork(K26_QUAD_LIMIT)
        real(c_double) :: work(K26_QUAD_LENW)

        call dqags(integrand_shim, a, b, epsabs, epsrel, &
                   result_loc, abserr_loc, neval_loc, ier_loc, &
                   K26_QUAD_LIMIT, K26_QUAD_LENW, last_loc, &
                   iwork, work)

        out_result = result_loc
        out_abserr = abserr_loc
        out_neval  = neval_loc
        out_ier    = ier_loc
    end subroutine k26astro_quad_dqags_call

    !> DQAGI bind(C) entry. Called by the C-side k26astro_quad_dqagi
    !> after the thread-local (fn, user) pair has been set.
    subroutine k26astro_quad_dqagi_call(bound, inf, epsabs, epsrel, &
                                         out_result, out_abserr, &
                                         out_neval, out_ier) &
        bind(C, name="k26astro_quad_dqagi_call")
        real(c_double), value, intent(in) :: bound, epsabs, epsrel
        integer(c_int), value, intent(in) :: inf
        real(c_double), intent(out)       :: out_result, out_abserr
        integer(c_int), intent(out)       :: out_neval, out_ier

        real(c_double) :: result_loc, abserr_loc
        integer        :: neval_loc, ier_loc, last_loc, inf_loc
        integer        :: iwork(K26_QUAD_LIMIT)
        real(c_double) :: work(K26_QUAD_LENW)

        inf_loc = inf

        call dqagi(integrand_shim, bound, inf_loc, epsabs, epsrel, &
                   result_loc, abserr_loc, neval_loc, ier_loc, &
                   K26_QUAD_LIMIT, K26_QUAD_LENW, last_loc, &
                   iwork, work)

        out_result = result_loc
        out_abserr = abserr_loc
        out_neval  = neval_loc
        out_ier    = ier_loc
    end subroutine k26astro_quad_dqagi_call

    !> Integrand shim — matches QUADPACK's `func` abstract interface
    !> (real(wp) function f(x); real(wp), intent(in) :: x). Forwards
    !> via the C trampoline that reads the thread-local (fn, user).
    !>
    !> Module-scope so both DQAGS and DQAGI entries can pass it.
    function integrand_shim(x) result(y)
        real(c_double), intent(in) :: x
        real(c_double)             :: y
        y = k26astro_quad_tls_trampoline(x)
    end function integrand_shim

end module k26astro_quad_iface
