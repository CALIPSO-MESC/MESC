!> Mock land-surface-model (LSM) driver for the MESC coupling API.
!!
!! Stands in for an LSM (e.g., ORCHIDEE) to exercise [[mesc_coupling_module]]
!! end-to-end:
!!
!! * builds a small synthetic grid of four contrasting land cells
!!   (wet tropics, temperate, boreal, arid subtropics),
!! * initialises MESC with time-invariant fields via [[mesc_coupling_init]]
!!   (a synthetic per-BGC parameter look-up table is generated on-the-fly),
!! * advances all cells with synthetic daily forcing via
!!   [[mesc_coupling_step]],
!! * retrieves the carbon pool state via [[mesc_coupling_get_cpool]],
!! * repeats the run from a mid-run pool snapshot via the `cpool0` argument
!!   (restart path) and checks the two runs agree,
!! * finishes with [[mesc_coupling_finalize]].
!!
!! Checks performed (any failure makes the program exit with a non-zero status):
!!
!! 1. all daily fluxes returned by [[mesc_coupling_step]] are finite, and the
!!    carbon input flux is non-negative;
!! 2. the daily carbon input flux equals the litter forcing passed in
!!    (closed carbon input path);
!! 3. the year-to-date flux accumulators equal the sums of the returned
!!    daily fluxes;
!! 4. carbon pools are finite and non-negative whenever they are retrieved;
!! 5. total soil carbon stays within generous physical bounds;
!! 6. the restart run reproduces the continuous run's pool state exactly.
!!
!! The program is self-contained: it generates its own synthetic parameter
!! table (deleted again at the end) and writes a summary table of annual
!! fluxes and stocks to standard output.
program mock_lsm_driver
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use precision_module, only : dp
  use mic_constant, only : mcpool, xrootorchidee, cnleaf2, cnwood2, cnroot2, &
                           ligleaf2, ligwood2, ligroot2
  use mesc_coupling_module, only : mesc_coupling_state, nday_year, &
                                   mesc_coupling_init, mesc_coupling_step, &
                                   mesc_coupling_get_cpool, &
                                   mesc_coupling_finalize
  implicit none

  ! ---- Grid definition (mimics an ORCHIDEE-forced global run, runcase 4) ----
  integer,  parameter :: ncell = 4  !! number of mock land grid cells
  integer,  parameter :: nlev  = 7  !! number of soil layers (runcase 4 layout)
  integer,  parameter :: npft  = 19 !! number of PFTs (ORCHIDEE)
  integer,  parameter :: nbgc  = 10 !! number of BGC parameter clusters
  integer,  parameter :: nyears = 10  !! years to simulate
  integer,  parameter :: nyear_restart = 5  !! year after which to branch the restart run
  real(dp), parameter :: zse(nlev) = [0.2_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.2_dp, &
                                      0.5_dp, 0.5_dp]
  !! soil layer thicknesses [m] (runcase 4 layout)

  ! Synthetic per-cell time-invariant fields: (1) wet tropics, (2) temperate,
  ! (3) boreal, (4) arid subtropics
  real(dp), parameter :: grid_lon(ncell) = [-60.0_dp, 5.0_dp, 130.0_dp, 25.0_dp]
  real(dp), parameter :: grid_lat(ncell) = [-3.0_dp, 45.0_dp, 62.0_dp, -25.0_dp]
  real(dp), parameter :: grid_area(ncell) = 1.0e5_dp  !! grid cell area [km2]
  integer,  parameter :: grid_pft(ncell) = [2, 6, 7, 10]
    !! dominant PFT per cell (tropical/temperate/boreal forest, grassland)
  real(dp), parameter :: grid_npp(ncell) = [1200.0_dp, 650.0_dp, 350.0_dp, 180.0_dp]
    !! annual NPP per cell [g C m-2 yr-1]
  real(dp), parameter :: grid_clay(ncell) = [0.35_dp, 0.22_dp, 0.15_dp, 0.25_dp]
  real(dp), parameter :: grid_silt(ncell) = [0.30_dp, 0.35_dp, 0.25_dp, 0.30_dp]
  real(dp), parameter :: grid_ph(ncell) = [5.2_dp, 6.5_dp, 4.8_dp, 7.8_dp]
  real(dp), parameter :: grid_bulkd(ncell) = [1050.0_dp, 1300.0_dp, 950.0_dp, 1550.0_dp]
    !! bulk density [kg m-3]
  integer,  parameter :: grid_isoil(ncell) = [2, 6, 8, 4]
  integer,  parameter :: grid_sorder(ncell) = [7, 3, 5, 2]
  integer,  parameter :: grid_bgctype(ncell) = [3, 5, 7, 9]

  ! Synthetic climate seasonality: soil temperature mean/amplitude and day of
  ! year of the temperature maximum per cell
  real(dp), parameter :: tsoil_mean(ncell) = [26.0_dp, 11.0_dp, 1.0_dp, 22.0_dp]  !! [degC]
  real(dp), parameter :: tsoil_amp(ncell) = [1.5_dp, 11.0_dp, 16.0_dp, 6.0_dp]    !! [degC]
  real(dp), parameter :: tday_peak(ncell) = [30.0_dp, 200.0_dp, 210.0_dp, 30.0_dp]
  ! Synthetic moisture seasonality (volumetric water content [m3 m-3])
  real(dp), parameter :: moist_mean(ncell) = [0.35_dp, 0.25_dp, 0.28_dp, 0.10_dp]
  real(dp), parameter :: moist_amp(ncell) = [0.05_dp, 0.08_dp, 0.10_dp, 0.04_dp]
  real(dp), parameter :: mday_peak(ncell) = [120.0_dp, 110.0_dp, 240.0_dp, 60.0_dp]
  ! Litterfall partitioning of NPP [1] and day of year of peak leaf litterfall
  real(dp), parameter :: frac_leaf(ncell) = [0.45_dp, 0.45_dp, 0.45_dp, 0.50_dp]
  real(dp), parameter :: frac_root(ncell) = [0.45_dp, 0.45_dp, 0.45_dp, 0.50_dp]
  real(dp), parameter :: frac_wood(ncell) = [0.10_dp, 0.10_dp, 0.10_dp, 0.00_dp]
  real(dp), parameter :: lday_peak(ncell) = [180.0_dp, 285.0_dp, 250.0_dp, 200.0_dp]
  real(dp), parameter :: lday_spread(ncell) = [60.0_dp, 35.0_dp, 30.0_dp, 25.0_dp]

  real(dp), parameter :: pi = 3.141592653589793_dp
  !> Depth attenuation of the soil temperature seasonal cycle per layer [1]
  real(dp), parameter :: depth_damp(nlev) = [1.00_dp, 0.98_dp, 0.95_dp, 0.90_dp, &
                                             0.85_dp, 0.70_dp, 0.55_dp]
  !> Tolerances for the checks
  real(dp), parameter :: tol_flux = 1.0e-6_dp     !! daily flux closure [g C m-2]
  real(dp), parameter :: tol_accum = 1.0e-6_dp    !! accumulator check [g C m-2]
  real(dp), parameter :: tol_pool = 1.0e-9_dp     !! restart pool check [mg C cm-3]
  real(dp), parameter :: soc_min = 1.0_dp         !! min plausible SOC [g C m-2]
  real(dp), parameter :: soc_max = 1.0e5_dp       !! max plausible SOC [g C m-2]

  ! Per-BGC parameter look-up table generated by this driver. NOTE: must be a
  ! fixed-length variable of at least 140 characters because getparam_global
  ! declares the file name as character(len=140).
  character(len=140), parameter :: fparam = "parameters_global_mock.csv"

  ! Local variables
  type(mesc_coupling_state) :: state
  real(dp) :: pools(nlev, mcpool)
  real(dp) :: run_a_final(ncell, nlev, mcpool)  !! pool state at end of run A
  real(dp) :: snapshot(ncell, nlev, mcpool)     !! pool state after year nyear_restart
  real(dp) :: leaf_weight(ncell, nday_year)     !! normalised leaf litter seasonality [1]
  real(dp) :: ytd_cinput(ncell), ytd_rsoil(ncell), ytd_cleach(ncell)
  real(dp) :: litter_leaf, litter_wood, litter_root
  real(dp) :: tsoil(nlev), moist(nlev), matpot(nlev)
  real(dp) :: flux_cinput, flux_rsoil, flux_cleach
  real(dp) :: soc, diff, maxdiff, tol
  integer  :: year, doy, icell, nfail
  logical  :: lx

  nfail = 0

  write(*, "(a)") "==============================================================="
  write(*, "(a)") " Mock LSM driver for the MESC coupling API"
  write(*, "(a,i0,a,i0,a,i0,a,i0)") &
    " Grid: ", ncell, " cells x ", nlev, " layers, ", npft, " PFTs, nbgc=", nbgc
  write(*, "(a,i0,a)") " Run A: ", nyears, " years (cold start)"
  write(*, "(a,i0,a,i0,a)") &
    " Run B: restart from year ", nyear_restart, " snapshot, ", &
    nyears - nyear_restart, " further years"

  call setup_litter_weights()
  call write_parameter_table(fparam)

  ! ------------------------------------------------------------------
  ! Run A: continuous coupled run from a cold start
  ! ------------------------------------------------------------------
  call mesc_coupling_init(state, ncell, nlev, npft, nbgc, zse, fparam, grid_lon, &
                          grid_lat, grid_area, grid_pft, xrootorchidee, cnleaf2, &
                          cnwood2, cnroot2, ligleaf2, ligwood2, ligroot2, &
                          grid_npp, grid_clay, grid_silt, grid_ph, grid_bulkd, &
                          grid_isoil, grid_sorder, grid_bgctype)

  write(*, "(a)") "---------------------------------------------------------------"
  write(*, "(a)") "Run A: annual fluxes [g C m-2] and total soil carbon [g C m-2]"
  write(*, "(a)") " year cell    Cinput      Rsoil      Cleach         SOC"

  do year = 1, nyears
    ytd_cinput(:) = 0.0_dp
    ytd_rsoil(:) = 0.0_dp
    ytd_cleach(:) = 0.0_dp

    do doy = 1, nday_year
      do icell = 1, ncell
        call daily_forcing(icell, doy, litter_leaf, litter_wood, litter_root, &
                           tsoil, moist, matpot)
        call mesc_coupling_step(state, icell, doy, year, litter_leaf, litter_wood, &
                                litter_root, tsoil, moist, matpot, flux_cinput, &
                                flux_rsoil, flux_cleach)
        call check_fluxes(icell, year, doy, litter_leaf + litter_wood + litter_root, &
                          flux_cinput, flux_rsoil, flux_cleach, nfail)
        ytd_cinput(icell) = ytd_cinput(icell) + flux_cinput
        ytd_rsoil(icell) = ytd_rsoil(icell) + flux_rsoil
        ytd_cleach(icell) = ytd_cleach(icell) + flux_cleach
      end do
    end do

    do icell = 1, ncell
      call mesc_coupling_get_cpool(state, icell, pools)
      call check_pools(icell, year, pools, nfail)

      soc = sum(pools * spread(1000.0_dp * zse, dim=2, ncopies=mcpool))
      if (soc < soc_min .or. soc > soc_max) then
        write(*, "(a,i0,a,i0,a,es12.4)") &
          "FAIL: year ", year, " cell ", icell, " total SOC outside plausible bounds: ", soc
        nfail = nfail + 1
      end if

      ! The year-to-date accumulators must equal the sums of the daily fluxes
      tol = tol_accum * max(1.0_dp, ytd_cinput(icell))
      if (abs(state%micoutput%fluxcinput(icell) - ytd_cinput(icell)) > tol) then
        write(*, "(a,i0,a,i0,a,2(es12.4,a))") "FAIL: year ", year, " cell ", icell, &
          " C input accumulator mismatch: ", state%micoutput%fluxcinput(icell), &
          " vs ", ytd_cinput(icell), " (accumulated dailies)"
        nfail = nfail + 1
      end if
      tol = tol_accum * max(1.0_dp, ytd_rsoil(icell))
      if (abs(state%micoutput%fluxrsoil(icell) - ytd_rsoil(icell)) > tol) then
        write(*, "(a,i0,a,i0,a,2(es12.4,a))") "FAIL: year ", year, " cell ", icell, &
          " respiration accumulator mismatch: ", state%micoutput%fluxrsoil(icell), &
          " vs ", ytd_rsoil(icell), " (accumulated dailies)"
        nfail = nfail + 1
      end if

      write(*, "(i5,1x,i4,4(1x,f11.3))") year, icell, state%micoutput%fluxcinput(icell), &
        state%micoutput%fluxrsoil(icell), ytd_cleach(icell), soc

      if (year == nyear_restart) snapshot(icell, :, :) = pools
      if (year == nyears) run_a_final(icell, :, :) = pools
    end do
  end do

  call mesc_coupling_finalize(state)

  ! ------------------------------------------------------------------
  ! Run B: restart from the year-`nyear_restart` pool snapshot via the
  ! `cpool0` initial-state argument, then continue to year `nyears`.
  ! ------------------------------------------------------------------
  write(*, "(a)") "---------------------------------------------------------------"
  write(*, "(a)") "Run B: restart run"

  call mesc_coupling_init(state, ncell, nlev, npft, nbgc, zse, fparam, grid_lon, &
                          grid_lat, grid_area, grid_pft, xrootorchidee, cnleaf2, &
                          cnwood2, cnroot2, ligleaf2, ligwood2, ligroot2, &
                          grid_npp, grid_clay, grid_silt, grid_ph, grid_bulkd, &
                          grid_isoil, grid_sorder, grid_bgctype, cpool0=snapshot)

  do year = nyear_restart + 1, nyears
    do doy = 1, nday_year
      do icell = 1, ncell
        call daily_forcing(icell, doy, litter_leaf, litter_wood, litter_root, &
                           tsoil, moist, matpot)
        call mesc_coupling_step(state, icell, doy, year, litter_leaf, litter_wood, &
                                litter_root, tsoil, moist, matpot, flux_cinput, &
                                flux_rsoil, flux_cleach)
        call check_fluxes(icell, year, doy, litter_leaf + litter_wood + litter_root, &
                          flux_cinput, flux_rsoil, flux_cleach, nfail)
      end do
    end do
  end do

  maxdiff = 0.0_dp
  do icell = 1, ncell
    call mesc_coupling_get_cpool(state, icell, pools)
    call check_pools(icell, nyears, pools, nfail)
    diff = maxval(abs(pools - run_a_final(icell, :, :)))
    maxdiff = max(maxdiff, diff)
  end do
  write(*, "(a,es12.4,a,es12.4)") &
    "Maximum pool-state difference (run B vs run A): ", maxdiff, &
    " mg C cm-3 (tolerance ", tol_pool
  if (maxdiff > tol_pool) then
    write(*, "(a)") "FAIL: restart run does not reproduce the continuous run"
    nfail = nfail + 1
  end if

  call mesc_coupling_finalize(state)
  call delete_file(fparam, lx)

  ! ------------------------------------------------------------------
  ! Summary
  ! ------------------------------------------------------------------
  write(*, "(a)") "==============================================================="
  if (nfail > 0) then
    write(*, "(a,i0,a)") "MOCK LSM DRIVER: FAIL (", nfail, " check(s) failed)"
    error stop 1
  else
    write(*, "(a)") "MOCK LSM DRIVER: PASS"
  end if

contains

  !> Precompute the normalised leaf-litter seasonality weights for each cell.
  !!
  !! A wrapped Gaussian (peak day `lday_peak`, spread `lday_spread`) is used,
  !! normalised so that the weights sum to one over the year. Root and wood
  !! litterfall are distributed uniformly over the year (see [[daily_forcing]]).
  subroutine setup_litter_weights()
    integer  :: jcell, jday
    real(dp) :: dist

    do jcell = 1, ncell
      do jday = 1, nday_year
        dist = abs(real(jday, dp) - lday_peak(jcell))
        dist = min(dist, real(nday_year, dp) - dist)  ! periodic distance
        leaf_weight(jcell, jday) = exp(-0.5_dp * (dist / lday_spread(jcell))**2)
      end do
      leaf_weight(jcell, :) = leaf_weight(jcell, :) / sum(leaf_weight(jcell, :))
    end do
  end subroutine setup_litter_weights

  !> Generate one day of synthetic forcing for cell `icell`.
  !!
  !! Temperature follows a cosine seasonal cycle, attenuated with depth;
  !! moisture follows a shifted cycle, clamped to [0.05, 0.75] m3 m-3; matric
  !! potential is tied to the moisture anomaly. Leaf litterfall follows the
  !! precomputed seasonal weights; root litterfall (kept strictly positive so
  !! that the litter C:N partitioning in MESC stays well-defined) and wood
  !! litterfall are uniform over the year. Annual totals equal the per-cell
  !! NPP by construction.
  subroutine daily_forcing(jcell, jday, dlit_leaf, dlit_wood, dlit_root, &
                           dtsoil, dmoist, dmatpot)
    integer,  intent(in)  :: jcell     !! grid cell index
    integer,  intent(in)  :: jday      !! day of year (1..365)
    real(dp), intent(out) :: dlit_leaf !! leaf litterfall today [g C m-2 day-1]
    real(dp), intent(out) :: dlit_wood !! wood litterfall today [g C m-2 day-1]
    real(dp), intent(out) :: dlit_root !! root litterfall today [g C m-2 day-1]
    real(dp), intent(out) :: dtsoil(nlev)  !! soil temperature per layer [degC]
    real(dp), intent(out) :: dmoist(nlev)  !! soil water content per layer [m3 m-3]
    real(dp), intent(out) :: dmatpot(nlev) !! soil matric potential per layer [MPa]

    integer  :: ns
    real(dp) :: phase, swc

    dlit_leaf = frac_leaf(jcell) * grid_npp(jcell) * leaf_weight(jcell, jday)
    dlit_wood = frac_wood(jcell) * grid_npp(jcell) / real(nday_year, dp)
    dlit_root = frac_root(jcell) * grid_npp(jcell) / real(nday_year, dp)

    phase = 2.0_dp * pi * (real(jday, dp) - tday_peak(jcell)) / real(nday_year, dp)
    do ns = 1, nlev
      dtsoil(ns) = tsoil_mean(jcell) + tsoil_amp(jcell) * cos(phase) * depth_damp(ns)
    end do

    phase = 2.0_dp * pi * (real(jday, dp) - mday_peak(jcell)) / real(nday_year, dp)
    swc = moist_mean(jcell) + moist_amp(jcell) * cos(phase)
    swc = min(0.75_dp, max(0.05_dp, swc))
    dmoist(:) = swc
    dmatpot(:) = -0.5_dp * (moist_mean(jcell) / swc)**2
  end subroutine daily_forcing

  !> Write a synthetic per-BGC parameter look-up table for [[getparam_global]].
  !!
  !! All scaling factors are set to 1.0, i.e. the default MESC
  !! parameterisation is used for every BGC type. In a production coupled run
  !! the LSM configuration would supply the calibrated table instead.
  subroutine write_parameter_table(fname)
    character(len=*), intent(in) :: fname !! look-up table file to write

    integer :: u, ibgc, n

    open(newunit=u, file=fname, status="replace", action="write")
    write(u, "(a)") "ibgc, xav, xak, xfm, xfs, xtvmic, xtvp, xtvc, xtvac, xkba," &
                    // " xqmaxcoeff, xdiffsoc, xnpp, xrootbeta, xvmaxbeta"
    do ibgc = 1, nbgc
      write(u, "(i0,*(1x,f4.1))") ibgc, (1.0_dp, n = 1, 14)
    end do
    close(u)
  end subroutine write_parameter_table

  !> Delete a file if it exists; `deleted` is true if a file was removed.
  subroutine delete_file(fname, deleted)
    character(len=*), intent(in)  :: fname   !! file to delete
    logical,          intent(out) :: deleted !! whether the file was deleted

    integer :: u

    inquire(file=fname, exist=deleted)
    if (deleted) then
      open(newunit=u, file=fname, status="old")
      close(u, status="delete")
    end if
  end subroutine delete_file

  !> Check the daily fluxes returned by [[mesc_coupling_step]].
  !!
  !! Fails (incrementing `nfail`) if any flux is not finite, if the carbon
  !! input flux is negative, or if the carbon input flux does not match the
  !! litter forcing (to within a small numerical tolerance). The latter
  !! verifies that MESC's internal C-input bookkeeping is closed.
  subroutine check_fluxes(jcell, jyear, jday, forcing_cinput, f_cinput, f_rsoil, &
                          f_cleach, nfail)
    integer,  intent(in)    :: jcell         !! grid cell index
    integer,  intent(in)    :: jyear         !! simulation year
    integer,  intent(in)    :: jday          !! day of year
    real(dp), intent(in)    :: forcing_cinput !! litter C forcing today [g C m-2]
    real(dp), intent(in)    :: f_cinput      !! C input flux returned [g C m-2]
    real(dp), intent(in)    :: f_rsoil       !! respiration flux returned [g C m-2]
    real(dp), intent(in)    :: f_cleach      !! leaching flux returned [g C m-2]
    integer,  intent(inout) :: nfail         !! running failure count

    if (.not. (ieee_is_finite(f_cinput) .and. ieee_is_finite(f_rsoil) .and. &
               ieee_is_finite(f_cleach))) then
      write(*, "(a,i0,a,i0,a,i0,a)") "FAIL: year ", jyear, " doy ", jday, &
        " cell ", jcell, ": non-finite flux returned"
      nfail = nfail + 1
      return
    end if

    if (f_cinput < -tol_flux) then
      write(*, "(a,i0,a,i0,a,i0,a,es12.4)") "FAIL: year ", jyear, " doy ", jday, &
        " cell ", jcell, ": negative C input flux ", f_cinput
      nfail = nfail + 1
    end if

    if (abs(f_cinput - forcing_cinput) > tol_flux * max(1.0_dp, forcing_cinput)) then
      write(*, "(a,i0,a,i0,a,i0,a,2(es12.4,a))") "FAIL: year ", jyear, " doy ", &
        jday, " cell ", jcell, ": C input flux ", f_cinput, &
        " does not match litter forcing ", forcing_cinput
      nfail = nfail + 1
    end if
  end subroutine check_fluxes

  !> Check carbon pools retrieved via [[mesc_coupling_get_cpool]].
  !!
  !! Fails (incrementing `nfail`) if any pool value is not finite or is
  !! negative.
  subroutine check_pools(jcell, jyear, jpools, nfail)
    integer,  intent(in)    :: jcell  !! grid cell index
    integer,  intent(in)    :: jyear  !! simulation year
    real(dp), intent(in)    :: jpools(nlev, mcpool) !! pools [mg C cm-3]
    integer,  intent(inout) :: nfail  !! running failure count

    if (.not. all(ieee_is_finite(jpools))) then
      write(*, "(a,i0,a,i0,a)") "FAIL: year ", jyear, " cell ", jcell, &
        ": non-finite pool state"
      nfail = nfail + 1
    end if

    if (any(jpools < 0.0_dp)) then
      write(*, "(a,i0,a,i0,a,es12.4)") "FAIL: year ", jyear, " cell ", jcell, &
        ": negative pool state, min = ", minval(jpools)
      nfail = nfail + 1
    end if
  end subroutine check_pools

end program mock_lsm_driver
