!> k26astro_atmos_nrlmsise_iface — ISO_C_BINDING surface for NRLMSISE-00.
!>
!> Exposes the upstream NRL FORTRAN `GTD7` and `METERS` entries as
!> plain-C symbols (no Fortran name mangling). The C-side wrapper in
!> atmos_nrlmsise.c calls k26astro_atmos_nrlmsise_init once at init
!> (sets the METERS-mode global to true so subsequent GTD7 returns
!> kg/m³ + /m³), then calls k26astro_atmos_nrlmsise_gtd7_call per
!> query.
!>
!> Build assumption: the upstream NRLMSISE-00.FOR has no IMPLICIT or
!> REAL*8 declarations and relies on default REAL typing. We compile
!> with `-fdefault-real-8 -fdefault-double-8` so default REAL = 8
!> bytes = c_double. This is the standard gfortran convention for
!> NRL FORTRAN codes (matches the brodo port's documented build).
!>
!> Determinism: the upstream uses COMMON blocks for switch state
!> (CSW) and the METERS global. These are process-global; no
!> thread-safe story. Callers that need concurrent NRLMSISE
!> evaluations must serialise externally. The K26 wrapper documents
!> this in atmos_nrlmsise.c.

module k26astro_atmos_nrlmsise_iface
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
    implicit none

    private

contains

    !> Init shim — sets the NRLMSISE METERS global to .TRUE. so all
    !> subsequent GTD7 calls return SI units (per-m³ + kg/m³). Idempotent
    !> at the FORTRAN side; callable repeatedly without harm.
    subroutine k26astro_atmos_nrlmsise_init_call() &
        bind(C, name="k26astro_atmos_nrlmsise_init_call")
        external :: meters
        logical  :: t
        t = .true.
        call meters(t)
    end subroutine k26astro_atmos_nrlmsise_init_call

    !> GTD7 shim. Wraps the NRL FORTRAN `GTD7(IYD,SEC,ALT,GLAT,GLONG,
    !> STL,F107A,F107,AP,MASS,D,T)` with a bind(C) signature.
    !>
    !> Inputs match the upstream's exactly:
    !>   iyd    : year+day as YYDDD (year ignored by current model)
    !>   sec    : UT seconds
    !>   alt    : altitude (km)
    !>   glat   : geodetic latitude (deg)
    !>   glong  : geodetic longitude (deg)
    !>   stl    : local apparent solar time (hr)
    !>   f107a  : 81-day average F10.7 flux
    !>   f107   : daily F10.7 flux (previous day)
    !>   ap(7)  : magnetic indices (daily Ap + 3hr Ap × 4 + 12-33hr avg
    !>            + 36-57hr avg) — passed by reference to match GTD7
    !>   mass   : mass selector (48 = all species + temperature)
    !>
    !> Outputs (passed-by-reference arrays):
    !>   out_d(9) : densities — He, O, N2, O2, Ar, total_mass, H, N,
    !>              anomalous_O. Units depend on METERS state set by
    !>              k26astro_atmos_nrlmsise_init_call (m^-3 + kg/m^3
    !>              when init has been called).
    !>   out_t(2) : exospheric temperature, temperature at alt (K).
    subroutine k26astro_atmos_nrlmsise_gtd7_call( &
            iyd, sec, alt, glat, glong, stl, &
            f107a, f107, ap, mass, &
            out_d, out_t) &
        bind(C, name="k26astro_atmos_nrlmsise_gtd7_call")
        integer(c_int), value, intent(in)  :: iyd, mass
        real(c_double), value, intent(in)  :: sec, alt, glat, glong, stl
        real(c_double), value, intent(in)  :: f107a, f107
        real(c_double),       intent(in)   :: ap(7)
        real(c_double),       intent(out)  :: out_d(9)
        real(c_double),       intent(out)  :: out_t(2)

        external :: gtd7
        real     :: sec_l, alt_l, glat_l, glong_l, stl_l, f107a_l, f107_l
        real     :: ap_l(7), d_l(9), t_l(2)
        integer  :: i

        ! Copy c_double inputs into default-REAL locals (which gfortran
        ! sizes to 8 bytes under -fdefault-real-8, matching c_double).
        sec_l   = sec
        alt_l   = alt
        glat_l  = glat
        glong_l = glong
        stl_l   = stl
        f107a_l = f107a
        f107_l  = f107
        do i = 1, 7
            ap_l(i) = ap(i)
        end do

        call gtd7(iyd, sec_l, alt_l, glat_l, glong_l, stl_l, &
                  f107a_l, f107_l, ap_l, mass, d_l, t_l)

        do i = 1, 9
            out_d(i) = d_l(i)
        end do
        do i = 1, 2
            out_t(i) = t_l(i)
        end do
    end subroutine k26astro_atmos_nrlmsise_gtd7_call

    !> GTD7D shim — variant of GTD7 where D(6) is the "effective total
    !> mass density for drag" including anomalous oxygen. Used for
    !> high-altitude satellite drag calculations.
    subroutine k26astro_atmos_nrlmsise_gtd7d_call( &
            iyd, sec, alt, glat, glong, stl, &
            f107a, f107, ap, mass, &
            out_d, out_t) &
        bind(C, name="k26astro_atmos_nrlmsise_gtd7d_call")
        integer(c_int), value, intent(in)  :: iyd, mass
        real(c_double), value, intent(in)  :: sec, alt, glat, glong, stl
        real(c_double), value, intent(in)  :: f107a, f107
        real(c_double),       intent(in)   :: ap(7)
        real(c_double),       intent(out)  :: out_d(9)
        real(c_double),       intent(out)  :: out_t(2)

        external :: gtd7d
        real     :: sec_l, alt_l, glat_l, glong_l, stl_l, f107a_l, f107_l
        real     :: ap_l(7), d_l(9), t_l(2)
        integer  :: i

        sec_l   = sec
        alt_l   = alt
        glat_l  = glat
        glong_l = glong
        stl_l   = stl
        f107a_l = f107a
        f107_l  = f107
        do i = 1, 7
            ap_l(i) = ap(i)
        end do

        call gtd7d(iyd, sec_l, alt_l, glat_l, glong_l, stl_l, &
                   f107a_l, f107_l, ap_l, mass, d_l, t_l)

        do i = 1, 9
            out_d(i) = d_l(i)
        end do
        do i = 1, 2
            out_t(i) = t_l(i)
        end do
    end subroutine k26astro_atmos_nrlmsise_gtd7d_call

end module k26astro_atmos_nrlmsise_iface
