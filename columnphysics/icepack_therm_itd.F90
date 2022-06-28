!=======================================================================
!
! Thermo calculations after call to coupler, related to ITD:
! ice thickness redistribution, lateral growth and melting.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First icepack_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then icepack_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Streamlined for efficiency by Elizabeth Hunke
! 2014: Column package created by Elizabeth Hunke
!
      module icepack_therm_itd

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c3, c4, c6, c10
      use icepack_parameters, only: p001, p1, p333, p5, p666, puny, bignum
      use icepack_parameters, only: rhos, rhoi, Lfresh, ice_ref_salinity
      use icepack_parameters, only: phi_init, dsin0_frazil, hs_ssl, salt_loss
      use icepack_parameters, only: rhosi, conserv_check, rhosmin
      use icepack_parameters, only: kitd, ktherm, heat_capacity
      use icepack_parameters, only: z_tracers, solve_zsal, hfrazilmin

      use icepack_tracers, only: ntrcr, nbtrcr
      use icepack_tracers, only: nt_qice, nt_qsno, nt_fbri, nt_sice
      use icepack_tracers, only: nt_apnd, nt_hpnd, nt_aero, nt_mp, nt_isosno, nt_isoice
      use icepack_tracers, only: nt_Tsfc, nt_iage, nt_FY, nt_fsd, nt_rhos
      use icepack_tracers, only: nt_alvl, nt_vlvl
      use icepack_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_snow
      use icepack_tracers, only: tr_iage, tr_FY, tr_lvl, tr_aero, tr_mp, tr_iso, tr_brine, tr_fsd
      use icepack_tracers, only: n_aero, n_mp, n_iso
      use icepack_tracers, only: bio_index

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_fsd, only: fsd_weld_thermo, icepack_cleanup_fsd,  get_subdt_fsd
      use icepack_itd, only: reduce_area, cleanup_itd
      use icepack_itd, only: aggregate_area, shift_ice
      use icepack_itd, only: column_sum, column_conservation_check
      use icepack_isotope, only: isoice_alpha, isotope_frac_method
      use icepack_mushy_physics, only: liquidus_temperature_mush, enthalpy_mush
      use icepack_therm_shared, only: hi_min
      use icepack_zbgc, only: add_new_ice_bgc
      use icepack_zbgc, only: lateral_melt_bgc

      implicit none

      private
      public :: linear_itd, &
                add_new_ice, &
                lateral_melt, &
                icepack_step_therm2

!=======================================================================

      contains

!=======================================================================
!
! ITD scheme that shifts ice among categories
!
! See Lipscomb, W. H.  Remapping the thickness distribution in sea
!     ice models. 2001, J. Geophys. Res., Vol 106, 13989--14000.
!
! Using the thermodynamic "velocities", interpolate to find the
! velocities in thickness space at the category boundaries, and
! compute the new locations of the boundaries.  Then for each
! category, compute the thickness distribution function,  g(h),
! between hL and hR, the left and right boundaries of the category.
! Assume g(h) is a linear polynomial that satisfies two conditions:
!
! (1) The ice area implied by g(h) equals aicen(n).
! (2) The ice volume implied by g(h) equals aicen(n)*hicen(n).
!
! Given g(h), at each boundary compute the ice area and volume lying
! between the original and new boundary locations.  Transfer area
! and volume across each boundary in the appropriate direction, thus
! restoring the original boundaries.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine linear_itd (ncat,        hin_max,     &
                             nilyr,       nslyr,       &
