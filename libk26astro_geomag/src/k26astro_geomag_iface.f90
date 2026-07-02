!> k26astro_geomag_iface — ISO_C_BINDING wrapper around IGRF-14.
!>
!> The upstream `igrf14syn` is a flat F77 subroutine (no module).
!> This module declares it via an interface block and exposes
!> bind(C) entries that C-side code can call directly.
!>
!> Coordinate convention transforms:
!>   - K26 takes geodetic latitude in radians; IGRF wants
!>     colatitude in degrees = 90 - lat_deg.
!>   - K26 takes longitude in radians (east-positive); IGRF wants
!>     east-longitude in degrees in [0, 360).
!>   - K26 altitude in metres; IGRF wants altitude in km
!>     (geodetic, itype=1).
!>   - K26 epoch is years past J2000.0; the C wrapper converts to
!>     AD-year before calling this entry (passes `date` directly).
!>
!> Output: x = north (nT), y = east (nT), z = vertical (nT, downward
!> positive). The C wrapper converts nT → T for SI consistency.

module k26astro_geomag_iface
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    private

    interface
        subroutine igrf14syn(isv, date, itype, alt, colat, elong, &
                              x, y, z, f)
            integer        :: isv, itype
            real(8)        :: date, alt, colat, elong, x, y, z, f
        end subroutine igrf14syn
    end interface

contains

    !> Main-field evaluation (isv=0). Inputs match the IGRF-14 native
    !> coordinate system: colat in degrees [0,180], elong in degrees
    !> [0,360), alt in km (geodetic with itype=1).
    subroutine k26astro_geomag_field_call(date, colat_deg, elong_deg, &
                                           alt_km, &
                                           out_x_nT, out_y_nT, &
                                           out_z_nT, out_f_nT) &
        bind(C, name="k26astro_geomag_field_call")
        real(c_double), value, intent(in) :: date, colat_deg, elong_deg, alt_km
        real(c_double), intent(out)       :: out_x_nT, out_y_nT, out_z_nT, out_f_nT

        integer :: isv, itype

        isv   = 0
        itype = 1   ! geodetic (WGS84)

        call igrf14syn(isv, date, itype, alt_km, colat_deg, elong_deg, &
                       out_x_nT, out_y_nT, out_z_nT, out_f_nT)
    end subroutine k26astro_geomag_field_call

    !> Secular-variation evaluation (isv=1). Same coordinate convention.
    !> Output components are nT/year. The `out_f_nT` slot is unused by
    !> IGRF in isv=1 mode (set to whatever igrf14syn returns; C-side
    !> wrapper recomputes total intensity from the magnitudes).
    subroutine k26astro_geomag_secvar_call(date, colat_deg, elong_deg, &
                                            alt_km, &
                                            out_x_nT, out_y_nT, &
                                            out_z_nT, out_f_nT) &
        bind(C, name="k26astro_geomag_secvar_call")
        real(c_double), value, intent(in) :: date, colat_deg, elong_deg, alt_km
        real(c_double), intent(out)       :: out_x_nT, out_y_nT, out_z_nT, out_f_nT

        integer :: isv, itype

        isv   = 1
        itype = 1

        call igrf14syn(isv, date, itype, alt_km, colat_deg, elong_deg, &
                       out_x_nT, out_y_nT, out_z_nT, out_f_nT)
    end subroutine k26astro_geomag_secvar_call

end module k26astro_geomag_iface
