!> \file modbulkmicrostat.f90
!!  Calculates profiles coming from the bulkmicrophysics


!>
!!  Calculates profiles coming from the bulkmicrophysics3
!>
!! Profiles coming from the bulkmicrophysics. Written to precep.expnr for the
!! rain rates etc., and to qlptend.expnr, nptend.expnr and qtptend.expnr for the
!! tendencies is rain water content, droplet number, and total water content,
!! respectively.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Olivier Geoffroy, KNMI
!!  \author Johan van de Dussen, TU Delft
!!  \author Jan Chylik, IGMK
!!  \author Jisk Attema, NLeSC
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2020 Delft University of Technology, Wageningen University, Utrecht University, KNMI, NLeSC
!
module modbulkmicrostat3
  use modglobal, only : longint
  use modlogging, only: finish

implicit none
private
character(len=*), parameter :: modname = 'modbulkmicrostat3'
PUBLIC  :: initbulkmicrostat3, bulkmicrostat3, exitbulkmicrostat3
save
!NetCDF variables
  integer,parameter                 :: nmvar = 55
  integer,parameter                 :: ntvar = 93
  integer                           :: ncid_mphys = 0
  integer                           :: ncid_tends = 0
  integer                           :: nrec_mphys = 0
  integer                           :: nrec_tends = 0
  ! Bug-fix note: ncinfo's intent(out) argument is an explicit-shape array of
  ! shape (4). Passing a non-contiguous row slice ncmname(k,:) of a (nmvar,4)
  ! array forces the compiler to create an array temporary and copy it back to
  ! the original strided locations after the call. On some platforms (arm64
  ! macOS, gfortran) the copy-back is unreliable, leaving ncmname blank and
  ! causing define_nc to abort with "Bad dimensional information".
  !
  ! Workaround: pass a local contiguous scratch (tmp4_m / tmp4_t) to ncinfo,
  ! then assign it into ncmname(k,:) with an ordinary array assignment. An
  ! array assignment is inlined as an element-wise loop by the compiler and
  ! does not go through the intent(out) copy-back path, so it works for
  ! non-contiguous LHS slices.
  character(80),dimension(nmvar,4)  :: ncmname
  character(80),dimension(ntvar,4)  :: nctname
  character(80),dimension(1,4)      :: tncmname
  character(80),dimension(1,4)      :: tnctname
  character(80)                     :: fname_mphys = 'mphysprofiles.xxx.nc'
  character(80)                     :: fname_tends = 'mphystendencies.xxx.nc'
  real                              :: dtav, timeav
  integer(kind=longint)             :: idtav, itimeav, tnext, tnextwrite
  integer                           :: nsamples
  logical                           :: lmicrostat = .false.

contains

!> Initialization routine, reads namelists and inits variables
subroutine initbulkmicrostat3
    use modmpi,    only  : myid, comm3d, mpierr, D_MPI_BCAST
    use fortran_support, only: nnml_output
    use modglobal, only  : ifnamopt, fname_options, cexpnr, &
         dtav_glob, timeav_glob, ladaptive, dtmax,btime,tres,checknamelisterror,kmax
    use modstat_nc, only : lnetcdf,open_nc,define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav
    use modmicrodata3

    implicit none
    integer      :: ierr
    integer      :: cnt
    character(80) :: tmp4_m(4), tmp4_t(4)  ! contiguous scratch for ncinfo (see header note)
    character(len=*), parameter :: routine = modname//'/initbulkmicrostat3'

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav

    dtav  = dtav_glob
    timeav  = timeav_glob
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMBULKMICROSTAT')
      write(nnml_output,NAMBULKMICROSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lmicrostat,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(timeav    ,1,0,comm3d,mpierr)
    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))

    if (idtav <= 0 .or. itimeav <= 0) then
      call finish(routine, 'dtav/timeav produce non-positive timestep counters (NAMBULKMICROSTAT)')
    end if

    tnext      = idtav   + btime
    tnextwrite = itimeav + btime
    nsamples = int(itimeav / idtav)

    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      call finish(routine, 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)')
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      call finish(routine, 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)')
    end if

    ! TODO: text output
    !if (myid == 0 .and. .not. lwarmstart) then
    !  open (ifoutput,file = 'precep.'//cexpnr ,status = 'replace')
    !  close(ifoutput)
    !end if

    if (lnetcdf) then
      ! initbulkmicrostat3 is called before initgenstat in startup.
      ! Only adopt genstat timing once those counters are valid.
      if (idtav_prof > 0 .and. itimeav_prof > 0) then
        idtav      = idtav_prof
        itimeav    = itimeav_prof
      end if
      tnext      = idtav   + btime
      tnextwrite = itimeav + btime
      nsamples = int(itimeav / idtav)
      if (nsamples <= 0) then
        call finish(routine, 'timeav must be >= dtav for microphysics statistics output')
      end if
      if (myid==0) then

        ! statisitcs output
        if (l_statistics) then
          fname_mphys(15:17) = cexpnr
          call nctiminfo(tncmname(1,:))

          cnt = 0
          call ncinfo(tmp4_m, 'imphys_freeze','change in th due to freezing'    ,'-','tt'); ncmname(cnt + imphys_freeze,:) = tmp4_m
          call ncinfo(tmp4_m, 'imphys_melt'  ,'change in th due to melting'     ,'-','tt'); ncmname(cnt + imphys_melt  ,:) = tmp4_m
          call ncinfo(tmp4_m, 'imphys_cond'  ,'change in th due to condensation','-','tt'); ncmname(cnt + imphys_cond  ,:) = tmp4_m
          call ncinfo(tmp4_m, 'imphys_ev'    ,'change in th due to evaporation' ,'-','tt'); ncmname(cnt + imphys_ev    ,:) = tmp4_m
          call ncinfo(tmp4_m, 'imphys_dep'   ,'change in th due to deposition'  ,'-','tt'); ncmname(cnt + imphys_dep   ,:) = tmp4_m
          call ncinfo(tmp4_m, 'imphys_sub'   ,'change in th due to sublimation' ,'-','tt'); ncmname(cnt + imphys_sub   ,:) = tmp4_m
          cnt = cnt + 6

          ! statistic_sv0_fsum(ncols,k1)
          call ncinfo(tmp4_m, 'n_hr','Average number content of rain'          ,'-','tt'); ncmname(cnt + in_hr,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hr','Average water content of rain'           ,'-','tt'); ncmname(cnt + iq_hr,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_cl','Average number content of cloud droplets','-','tt'); ncmname(cnt + in_cl,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_cl','Average water content of cloud droplets' ,'-','tt'); ncmname(cnt + iq_cl,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_ci','Average number content of ice crystals'  ,'-','tt'); ncmname(cnt + in_ci,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_ci','Average water content of ice crystals'   ,'-','tt'); ncmname(cnt + iq_ci,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_hs','Average number content of snow'          ,'-','tt'); ncmname(cnt + in_hs,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hs','Average water content of snow'           ,'-','tt'); ncmname(cnt + iq_hs,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_hg','Average number content of graupel'       ,'-','tt'); ncmname(cnt + in_hg,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hg','Average water content of graupel'        ,'-','tt'); ncmname(cnt + iq_hg,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_cc','Average number content of ccn'           ,'-','tt'); ncmname(cnt + in_cc,:) = tmp4_m
          cnt = cnt + 11

          ! statistic_svp_fsum(ncols,k1)
          call ncinfo(tmp4_m, 'n_hrp','Average tendency of number content of rain'          ,'-','tt'); ncmname(cnt + in_hr,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hrp','Average tendency of water content of rain'           ,'-','tt'); ncmname(cnt + iq_hr,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_clp','Average tendency of number content of cloud droplets','-','tt'); ncmname(cnt + in_cl,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_clp','Average tendency of water content of cloud droplets' ,'-','tt'); ncmname(cnt + iq_cl,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_cip','Average tendency of number content of ice crystals'  ,'-','tt'); ncmname(cnt + in_ci,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_cip','Average tendency of water content of ice crystals'   ,'-','tt'); ncmname(cnt + iq_ci,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_hsp','Average tendency of number content of snow'          ,'-','tt'); ncmname(cnt + in_hs,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hsp','Average tendency of water content of snow'           ,'-','tt'); ncmname(cnt + iq_hs,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_hgp','Average tendency of number content of graupel'       ,'-','tt'); ncmname(cnt + in_hg,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hgp','Average tendency of water content of graupel'        ,'-','tt'); ncmname(cnt + iq_hg,:) = tmp4_m
          call ncinfo(tmp4_m, 'n_ccp','Average tendency of number content of ccn'           ,'-','tt'); ncmname(cnt + in_cc,:) = tmp4_m
          cnt = cnt + 11

          ! statistic_sv0_count(ncols,k1)
          call ncinfo(tmp4_m, 'q_hr_count','Count of water content of rain above threshold'           ,'-','tt'); ncmname(cnt + 1,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_cl_count','Count of water content of cloud droplets above threshold' ,'-','tt'); ncmname(cnt + 2,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_ci_count','Count of water content of ice crystals above threshold'   ,'-','tt'); ncmname(cnt + 3,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hs_count','Count of water content of snow above threshold'           ,'-','tt'); ncmname(cnt + 4,:) = tmp4_m
          call ncinfo(tmp4_m, 'q_hg_count','Count of water content of graupel above threshold'        ,'-','tt'); ncmname(cnt + 5,:) = tmp4_m
          cnt = cnt + 5

          ! statistic_sv0_csum(ncols,k1)
          call ncinfo(tmp4_m, 'cn_hr','Average number content of rain above threshold'          ,'-','tt'); ncmname(cnt + 1,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hr','Average water content of rain above threshold'           ,'-','tt'); ncmname(cnt + 2,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_cl','Average number content of cloud droplets above threshold','-','tt'); ncmname(cnt + 3,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_cl','Average water content of cloud droplets above threshold' ,'-','tt'); ncmname(cnt + 4,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_ci','Average number content of ice crystals above threshold'  ,'-','tt'); ncmname(cnt + 5,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_ci','Average water content of ice crystals above threshold'   ,'-','tt'); ncmname(cnt + 6,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_hs','Average number content of snow above threshold'          ,'-','tt'); ncmname(cnt + 7,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hs','Average water content of snow above threshold'           ,'-','tt'); ncmname(cnt + 8,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_hg','Average number content of graupel above threshold'       ,'-','tt'); ncmname(cnt + 9,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hg','Average water content of graupel above threshold'        ,'-','tt'); ncmname(cnt +10,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_cc','Average number content of ccn above threshold'           ,'-','tt'); ncmname(cnt +11,:) = tmp4_m
          cnt = cnt + 11

          ! statistic_svp_csum(ncols,k1)
          call ncinfo(tmp4_m, 'cn_hrp','Average tendency of number content of rain above threshold'          ,'-','tt'); ncmname(cnt + 1,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hrp','Average tendency of water content of rain above threshold'           ,'-','tt'); ncmname(cnt + 2,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_clp','Average tendency of number content of cloud droplets above threshold','-','tt'); ncmname(cnt + 3,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_clp','Average tendency of water content of cloud droplets above threshold' ,'-','tt'); ncmname(cnt + 4,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_cip','Average tendency of number content of ice crystals above threshold'  ,'-','tt'); ncmname(cnt + 5,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_cip','Average tendency of water content of ice crystals above threshold'   ,'-','tt'); ncmname(cnt + 6,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_hsp','Average tendency of number content of snow above threshold'          ,'-','tt'); ncmname(cnt + 7,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hsp','Average tendency of water content of snow above threshold'           ,'-','tt'); ncmname(cnt + 8,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_hgp','Average tendency of number content of graupel above threshold'       ,'-','tt'); ncmname(cnt + 9,:) = tmp4_m
          call ncinfo(tmp4_m, 'cq_hgp','Average tendency of water content of graupel above threshold'        ,'-','tt'); ncmname(cnt +10,:) = tmp4_m
          call ncinfo(tmp4_m, 'cn_ccp','Average tendency of number content of ccn above threshold'           ,'-','tt'); ncmname(cnt +11,:) = tmp4_m
          cnt = cnt + 11

          ! -- finish
          if (cnt /= nmvar) then
            call finish(routine, 'wrong number of mphys output variables')
          endif

          ! adding dimensions
          call open_nc(fname_mphys,ncid_mphys,nrec_mphys,n3=kmax)
          if (nrec_mphys == 0) then
            call define_nc(ncid_mphys, 1, tncmname)
            call writestat_dims_nc(ncid_mphys)
          end if
          call define_nc(ncid_mphys, nmvar, ncmname)
        endif ! l_statistics

        ! tendencies output
        if (l_tendencies) then
          fname_tends(17:19) = cexpnr
          call nctiminfo(tnctname(1,:))

          ! all tend_fsum
          cnt = 0
          call ncinfo(tmp4_t, 'dn_cl_nu'      , 'droplet nucleation rate'                                      ,'-','tt'); nctname(cnt + idn_cl_nu     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_inu'     , 'ice nucleation rate'                                          ,'-','tt'); nctname(cnt + idn_ci_inu    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_au'      , 'change in number of cloud droplets due to autoconversion'     ,'-','tt'); nctname(cnt + idn_cl_au     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_au'      , 'change in mass of raindrops due to autoconversion'            ,'-','tt'); nctname(cnt + idq_hr_au     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_au'      , 'change in number of raindrops due to autoconversion'          ,'-','tt'); nctname(cnt + idn_hr_au     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_ac'      , 'change in mass of raindrops due to accretion'                 ,'-','tt'); nctname(cnt + idq_hr_ac     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_ac'      , 'change in number of cloud droplets due to accretion'          ,'-','tt'); nctname(cnt + idn_cl_ac     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_br'      , 'change in number of raindrops due to breakup'                 ,'-','tt'); nctname(cnt + idn_hr_br     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_sc'      , 'change in number of raindrops due to self-collection'         ,'-','tt'); nctname(cnt + idn_hr_sc     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_ev'      , 'change in mass of raindrops due to evaporation'               ,'-','tt'); nctname(cnt + idq_hr_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_ev'      , 'change in number of raindrops due to evaporation'             ,'-','tt'); nctname(cnt + idn_hr_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_dep'     , 'deposition rate for clouds'                                   ,'-','tt'); nctname(cnt + idq_ci_dep    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_dep'     , 'deposition rate for snow'                                     ,'-','tt'); nctname(cnt + idq_hs_dep    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_dep'     , 'deposition rate for graupel'                                  ,'-','tt'); nctname(cnt + idq_hg_dep    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_rime'    , 'riming growth of ice'                                         ,'-','tt'); nctname(cnt + idq_ci_rime   ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_rime_ci' , ' - and impact on n_cl'                                        ,'-','tt'); nctname(cnt + idn_cl_rime_ci,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_rime'    , 'riming growth of snow'                                        ,'-','tt'); nctname(cnt + idq_hs_rime   ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_rime_hs' , ' - and impact on n_cl'                                        ,'-','tt'); nctname(cnt + idn_cl_rime_hs,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_rime'    , 'riming growth for graupel'                                    ,'-','tt'); nctname(cnt + idq_hg_rime   ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_rime_hg' , ' - and impact on n_cl'                                        ,'-','tt'); nctname(cnt + idn_cl_rime_hg,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hshr_rime'  , 'riming growth for snow with rain'                             ,'-','tt'); nctname(cnt + idq_hshr_rime ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_rime_hs' , ' - and impact on n_hr'                                        ,'-','tt'); nctname(cnt + idn_hr_rime_hs,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hghr_rime'  , 'riming growth for graupel with rain'                          ,'-','tt'); nctname(cnt + idq_hghr_rime ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_rime_hg' , ' - and impact on n_hr'                                        ,'-','tt'); nctname(cnt + idn_hr_rime_hg,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_rime_ri' , 'rain loss from riming of ice+rain->gr'                        ,'-','tt'); nctname(cnt + idq_hr_rime_ri,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_rime_ri' , 'ice loss from riming of ice+rain->gr'                         ,'-','tt'); nctname(cnt + idq_ci_rime_ri,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_rime_ri' , 'ice number loss from riming of ice+rain->gr'                  ,'-','tt'); nctname(cnt + idn_ci_rime_ri,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_rime_ri' , 'rain number loss from riming of ice+rain->gr'                 ,'-','tt'); nctname(cnt + idn_hr_rime_ri,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_col_rs'  , 'rain loss from riming of ice+snow->gr'                        ,'-','tt'); nctname(cnt + idq_hr_col_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_col_rs'  , 'rain number loss from riming of ice+snow->gr'                 ,'-','tt'); nctname(cnt + idq_hs_col_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_col_rs'  , 'snow loss from riming of ice+snow->gr'                        ,'-','tt'); nctname(cnt + idn_hr_col_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_col_rs'  , 'snow number loss from riming of ice+snow->gr'                 ,'-','tt'); nctname(cnt + idn_hs_col_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_col_ri'  , 'rain loss from riming of ice+rain->gr'                        ,'-','tt'); nctname(cnt + idq_hr_col_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_col_ri'  , 'ice loss from riming of ice+rain->gr'                         ,'-','tt'); nctname(cnt + idq_ci_col_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_col_ri'  , 'ice number loss from riming of ice+rain->gr'                  ,'-','tt'); nctname(cnt + idn_ci_col_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_col_ri'  , 'rain number loss from riming of ice+rain->gr'                 ,'-','tt'); nctname(cnt + idn_hr_col_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_cl_het'     , 'heterogeneou freezing of cloud water'                         ,'-','tt'); nctname(cnt + idq_cl_het    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_het'     , 'heterogeneou freezing of cloud water'                         ,'-','tt'); nctname(cnt + idn_cl_het    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_het'     , 'heterogeneou freezing of raindrops'                           ,'-','tt'); nctname(cnt + idq_hr_het    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_het'     , 'heterogeneou freezing of raindrops'                           ,'-','tt'); nctname(cnt + idn_hr_het    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_cl_hom'     , 'homogeneous freezing of cloud water'                          ,'-','tt'); nctname(cnt + idq_cl_hom    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_hom'     , 'homogeneous freezing of cloud water'                          ,'-','tt'); nctname(cnt + idn_cl_hom    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_col_iis' , 'self-collection of cloud ice'                                 ,'-','tt'); nctname(cnt + idq_ci_col_iis,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_col_iis' , 'self-collection of cloud ice'                                 ,'-','tt'); nctname(cnt + idn_ci_col_iis,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_col_sss' , 'self-collection of snow'                                      ,'-','tt'); nctname(cnt + idn_hs_col_sss,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hsci_col'   , 'collection s+i - trend in q_hs'                               ,'-','tt'); nctname(cnt + idq_hsci_col  ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_col_hs'  , 'collection s+i - trend in n_ci'                               ,'-','tt'); nctname(cnt + idn_ci_col_hs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hghs_col'   , 'collection g+s - trend in q_hg'                               ,'-','tt'); nctname(cnt + idq_hghs_col  ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_col_hg'  , 'collection g+s - trend in n_hs'                               ,'-','tt'); nctname(cnt + idn_hs_col_hg ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_cv'      , 'partial conversion ice -> graupel'                            ,'-','tt'); nctname(cnt + idq_ci_cv     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_cv'      , 'partial conversion ice -> graupel'                            ,'-','tt'); nctname(cnt + idn_ci_cv     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_cv'      , 'partial conversion snow-> graupel'                            ,'-','tt'); nctname(cnt + idq_hs_cv     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_cv'      , 'partial conversion snow-> graupel'                            ,'-','tt'); nctname(cnt + idn_hs_cv     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_sc'      , 'cloud self-collection'                                        ,'-','tt'); nctname(cnt + idn_cl_sc     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_mul'     , 'ice multiplication'                                           ,'-','tt'); nctname(cnt + idn_ci_mul    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_mul'     , 'ice multiplication'                                           ,'-','tt'); nctname(cnt + idq_ci_mul    ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_me'      , 'number tendency melting of cloud ice'                         ,'-','tt'); nctname(cnt + idn_ci_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_me'      , 'mass tendency melting of cloud ice'                           ,'-','tt'); nctname(cnt + idq_ci_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_me'      , 'number tendency melting of snow'                              ,'-','tt'); nctname(cnt + idn_hs_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_me'      , 'mass tendency melting of snow'                                ,'-','tt'); nctname(cnt + idq_hs_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hg_me'      , 'number tendency melting of graupel'                           ,'-','tt'); nctname(cnt + idn_hg_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_me'      , 'mass tendency melting of graupel'                             ,'-','tt'); nctname(cnt + idq_hg_me     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_ev'      , 'number tendency evaporation of cloud ice'                     ,'-','tt'); nctname(cnt + idn_ci_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_ev'      , 'mass tendency evaporation of cloud ice'                       ,'-','tt'); nctname(cnt + idq_ci_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_ev'      , 'number tendency evaporation of snow'                          ,'-','tt'); nctname(cnt + idn_hs_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_ev'      , 'mass tendency evaporation of snow'                            ,'-','tt'); nctname(cnt + idq_hs_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hg_ev'      , 'number tendency evaporation of graupel'                       ,'-','tt'); nctname(cnt + idn_hg_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_ev'      , 'mass tendency evaporation of graupel'                         ,'-','tt'); nctname(cnt + idq_hg_ev     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_eme_ic'  , 'number tendency enhanced melting of cloud ice by cloud water' ,'-','tt'); nctname(cnt + idn_ci_eme_ic ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_eme_ic'  , 'mass tendency enhanced melting of cloud ice by cloud water'   ,'-','tt'); nctname(cnt + idq_ci_eme_ic ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_eme_ri'  , 'number tendency enhanced melting of cloud ice by rain'        ,'-','tt'); nctname(cnt + idn_ci_eme_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_eme_ri'  , 'mass tendency enhanced melting of cloud ice  by rain'         ,'-','tt'); nctname(cnt + idq_ci_eme_ri ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_eme_sc'  , 'number tendency enhanced melting of snow by cloud water'      ,'-','tt'); nctname(cnt + idn_hs_eme_sc ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_eme_sc'  , 'mass tendency enhanced melting of snow by cloud water'        ,'-','tt'); nctname(cnt + idq_hs_eme_sc ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_eme_rs'  , 'number tendency enhanced melting of snow by rain'             ,'-','tt'); nctname(cnt + idn_hs_eme_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_eme_rs'  , 'mass tendency enhanced melting of snow by rain'               ,'-','tt'); nctname(cnt + idq_hs_eme_rs ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hg_eme_gc'  , 'number tendency enhanced melting of graupel by liquid clouds' ,'-','tt'); nctname(cnt + idn_hg_eme_gc ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_eme_gc'  , 'mass tendency enhanced melting of graupel by liquid clouds'   ,'-','tt'); nctname(cnt + idq_hg_eme_gc ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hg_eme_gr'  , 'number tendency enhanced melting of graupel by rain'          ,'-','tt'); nctname(cnt + idn_hg_eme_gr ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_eme_gr'  , 'mass tendency enhanced melting of graupel by rain'            ,'-','tt'); nctname(cnt + idq_hg_eme_gr ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_se'      , 'sedimentation for clouds water - number'                      ,'-','tt'); nctname(cnt + idn_cl_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_cl_se'      , '     -||-- mixing ration'                                     ,'-','tt'); nctname(cnt + idq_cl_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_ci_se'      , 'sedimentation for cloud ice - number'                         ,'-','tt'); nctname(cnt + idn_ci_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_ci_se'      , '      -||-- mixing ration'                                    ,'-','tt'); nctname(cnt + idq_ci_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hr_se'      , 'sedimentation for rain - number'                              ,'-','tt'); nctname(cnt + idn_hr_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hr_se'      , '      -||-- mixing ration'                                    ,'-','tt'); nctname(cnt + idq_hr_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hs_se'      , 'sedimentation for snow - number'                              ,'-','tt'); nctname(cnt + idn_hs_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hs_se'      , '      -||-- mixing ration'                                    ,'-','tt'); nctname(cnt + idq_hs_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_hg_se'      , 'sedimentation for graupel - number'                           ,'-','tt'); nctname(cnt + idn_hg_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_hg_se'      , '      -||-- mixing ration'                                    ,'-','tt'); nctname(cnt + idq_hg_se     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dq_cl_sa'      , 'saturation adjustment'                                        ,'-','tt'); nctname(cnt + idq_cl_sa     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'dn_cl_sa'      , 'change in n_cl due to saturation adjustment'                  ,'-','tt'); nctname(cnt + idn_cl_sa     ,:) = tmp4_t
          call ncinfo(tmp4_t, 'ret_cc'        , 'recovery of ccn'                                              ,'-','tt'); nctname(cnt + iret_cc       ,:) = tmp4_t
          cnt = cnt + 93

          ! -- finish
          if (cnt /= ntvar) then
            call finish(routine,  'wrong number of output tendencies')
          endif

          ! adding dimensions
          call open_nc(fname_tends,ncid_tends,nrec_tends,n3=kmax)
          if (nrec_tends == 0) then
            call define_nc(ncid_tends, 1, tnctname)
            call writestat_dims_nc(ncid_tends)
          end if
          call define_nc(ncid_tends, ntvar, nctname)
        endif ! l_tendencies
      end if
   end if
  end subroutine initbulkmicrostat3

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat3(l_write, l_sample)
    use modglobal,    only  : rk3step, timee, dt_lim
    implicit none

    ! Sampling the output has moved from dobulkmicrostat3 to bulkmicrostat3,
    ! because the data is more easily available there.
    ! To keep all the time keeping in one place, we do 2 calls per step to microstat3
    ! with the following parameters:
    logical, intent(out) :: l_sample ! If bulkmicro3 should sample the tendencies
    logical, intent(in)  :: l_write  ! If bulkmicro3 has sampled the tendencies

    ! by default, don't sample anything
    l_sample = .false.

    if (.not. lmicrostat) return
    if (rk3step /= 3) return
    if (timee == 0) return

    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if

    if (timee >= tnext) then
      tnext = tnext + idtav
      ! NOTE: this used to be "call dobulkmicrostat3" to average over (i,j) and sums over subdomains
      ! 1. the sum over (i,j) now happens in modbulkmicro3 if l_sample = .true.
      ! 2. the MPI call to sum over subdomains moved to writebulkmicrostat3
      l_sample = .true.
    endif

    ! only write if bulkmicro3 has sampled (ie. l_write = .true.)
    if (timee >= tnextwrite .and. l_write) then
      tnextwrite = tnextwrite + itimeav
      call writebulkmicrostat3
    end if

  end subroutine bulkmicrostat3


!------------------------------------------------------------------------------!
!> Write the stats to file
!>  * write a text file
!>  * write to a netcdf file using writestat_nc
  subroutine writebulkmicrostat3
    use modmpi,     only : myid, MPI_SUM, comm3d, mpierr, D_MPI_ALLREDUCE
    use modglobal,  only : rtimee, kmax, ijtot
    use modstat_nc, only : lnetcdf, writestat_nc
!    use modgenstat, only : nrec_prof=>nrec
    use modmicrodata3

    implicit none
    character(len=*), parameter :: routine = modname//'/writebulkmicrostat3'
    integer  :: nsecs, nhrs, nminut
    integer  :: cnt

    nsecs  = nint(rtimee)
    nhrs   = int (nsecs/3600)
    nminut = int (nsecs/60)-nhrs*60
    nsecs  = mod (nsecs,60)

    ! gather all columns
    if (l_statistics) then
      call D_MPI_ALLREDUCE(statistic_mphys    ,size(statistic_mphys    ),MPI_SUM,comm3d,mpierr)
      call D_MPI_ALLREDUCE(statistic_sv0_count,size(statistic_sv0_count),MPI_SUM,comm3d,mpierr)
      call D_MPI_ALLREDUCE(statistic_sv0_fsum ,size(statistic_sv0_fsum ),MPI_SUM,comm3d,mpierr)
      call D_MPI_ALLREDUCE(statistic_sv0_csum ,size(statistic_sv0_csum ),MPI_SUM,comm3d,mpierr)
      call D_MPI_ALLREDUCE(statistic_svp_fsum ,size(statistic_svp_fsum ),MPI_SUM,comm3d,mpierr)
      call D_MPI_ALLREDUCE(statistic_svp_csum ,size(statistic_svp_csum ),MPI_SUM,comm3d,mpierr)

      ! normalize
      if (nsamples <= 0) then
        call finish(routine, 'nsamples <= 0 before normalization in writebulkmicrostat3')
      end if
      if (ijtot <= 0.) then
        call finish(routine, 'ijtot <= 0 before normalization in writebulkmicrostat3')
      end if
      statistic_sv0_fsum = statistic_sv0_fsum / ijtot / nsamples
      statistic_svp_fsum = statistic_svp_fsum / ijtot / nsamples
      where(statistic_sv0_count .gt. 0)
        statistic_sv0_csum = statistic_sv0_csum / statistic_sv0_count / nsamples
        statistic_svp_csum = statistic_svp_csum / statistic_sv0_count / nsamples
      elsewhere
        statistic_sv0_csum = 0.
        statistic_svp_csum = 0.
      endwhere
      statistic_sv0_count = statistic_sv0_count / ijtot / nsamples

      call sanitize_nc_2d(statistic_mphys,     'statistic_mphys')
      call sanitize_nc_2d(statistic_sv0_count, 'statistic_sv0_count')
      call sanitize_nc_2d(statistic_sv0_fsum,  'statistic_sv0_fsum')
      call sanitize_nc_2d(statistic_sv0_csum,  'statistic_sv0_csum')
      call sanitize_nc_2d(statistic_svp_fsum,  'statistic_svp_fsum')
      call sanitize_nc_2d(statistic_svp_csum,  'statistic_svp_csum')
    endif

    if (l_tendencies) then
      call D_MPI_ALLREDUCE(tend_fsum, size(tend_fsum), MPI_SUM, comm3d, mpierr)

      ! normalize
      tend_fsum = tend_fsum / ijtot / nsamples
      call sanitize_nc_2d(tend_fsum, 'tend_fsum')
    endif

    if (myid == 0) then

      ! TODO: text output

      if (lnetcdf) then
        ! TODO units and scaling factors for output

        if (l_statistics) then
          call writestat_nc(ncid_mphys,1,tncmname,(/rtimee/),nrec_mphys,.true.)
          cnt = 0
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_freeze:imphys_freeze,:),statistic_mphys(imphys_freeze:imphys_freeze,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_melt  :imphys_melt  ,:),statistic_mphys(imphys_melt  :imphys_melt  ,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_cond  :imphys_cond  ,:),statistic_mphys(imphys_cond  :imphys_cond  ,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_ev    :imphys_ev    ,:),statistic_mphys(imphys_ev    :imphys_ev    ,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_dep   :imphys_dep   ,:),statistic_mphys(imphys_dep   :imphys_dep   ,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + imphys_sub   :imphys_sub   ,:),statistic_mphys(imphys_sub   :imphys_sub   ,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 6

          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hr :  cnt + in_hr,:),statistic_sv0_fsum(in_hr :  in_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hr :  cnt + iq_hr,:),statistic_sv0_fsum(iq_hr :  iq_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cl :  cnt + in_cl,:),statistic_sv0_fsum(in_cl :  in_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_cl :  cnt + iq_cl,:),statistic_sv0_fsum(iq_cl :  iq_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_ci :  cnt + in_ci,:),statistic_sv0_fsum(in_ci :  in_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_ci :  cnt + iq_ci,:),statistic_sv0_fsum(iq_ci :  iq_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hs :  cnt + in_hs,:),statistic_sv0_fsum(in_hs :  in_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hs :  cnt + iq_hs,:),statistic_sv0_fsum(iq_hs :  iq_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hg :  cnt + in_hg,:),statistic_sv0_fsum(in_hg :  in_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hg :  cnt + iq_hg,:),statistic_sv0_fsum(iq_hg :  iq_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cc :  cnt + in_cc,:),statistic_sv0_fsum(in_cc :  in_cc,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 11

          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hr : cnt + in_hr,:),statistic_svp_fsum(in_hr : in_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hr : cnt + iq_hr,:),statistic_svp_fsum(iq_hr : iq_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cl : cnt + in_cl,:),statistic_svp_fsum(in_cl : in_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_cl : cnt + iq_cl,:),statistic_svp_fsum(iq_cl : iq_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_ci : cnt + in_ci,:),statistic_svp_fsum(in_ci : in_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_ci : cnt + iq_ci,:),statistic_svp_fsum(iq_ci : iq_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hs : cnt + in_hs,:),statistic_svp_fsum(in_hs : in_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hs : cnt + iq_hs,:),statistic_svp_fsum(iq_hs : iq_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hg : cnt + in_hg,:),statistic_svp_fsum(in_hg : in_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hg : cnt + iq_hg,:),statistic_svp_fsum(iq_hg : iq_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cc : cnt + in_cc,:),statistic_svp_fsum(in_cc : in_cc,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 11

          call writestat_nc(ncid_mphys,1,ncmname(cnt + 1 : cnt + 1,:),statistic_sv0_count(iq_hr : iq_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + 2 : cnt + 2,:),statistic_sv0_count(iq_cl : iq_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + 3 : cnt + 3,:),statistic_sv0_count(iq_ci : iq_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + 4 : cnt + 4,:),statistic_sv0_count(iq_hs : iq_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + 5 : cnt + 5,:),statistic_sv0_count(iq_hg : iq_hg,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 5

          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hr : cnt + in_hr,:),statistic_sv0_csum(in_hr : in_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hr : cnt + iq_hr,:),statistic_sv0_csum(iq_hr : iq_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cl : cnt + in_cl,:),statistic_sv0_csum(in_cl : in_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_cl : cnt + iq_cl,:),statistic_sv0_csum(iq_cl : iq_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_ci : cnt + in_ci,:),statistic_sv0_csum(in_ci : in_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_ci : cnt + iq_ci,:),statistic_sv0_csum(iq_ci : iq_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hs : cnt + in_hs,:),statistic_sv0_csum(in_hs : in_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hs : cnt + iq_hs,:),statistic_sv0_csum(iq_hs : iq_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hg : cnt + in_hg,:),statistic_sv0_csum(in_hg : in_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hg : cnt + iq_hg,:),statistic_sv0_csum(iq_hg : iq_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cc : cnt + in_cc,:),statistic_sv0_csum(in_cc : in_cc,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 11

          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hr : cnt + in_hr,:),statistic_svp_csum(in_hr : in_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hr : cnt + iq_hr,:),statistic_svp_csum(iq_hr : iq_hr,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cl : cnt + in_cl,:),statistic_svp_csum(in_cl : in_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_cl : cnt + iq_cl,:),statistic_svp_csum(iq_cl : iq_cl,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_ci : cnt + in_ci,:),statistic_svp_csum(in_ci : in_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_ci : cnt + iq_ci,:),statistic_svp_csum(iq_ci : iq_ci,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hs : cnt + in_hs,:),statistic_svp_csum(in_hs : in_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hs : cnt + iq_hs,:),statistic_svp_csum(iq_hs : iq_hs,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_hg : cnt + in_hg,:),statistic_svp_csum(in_hg : in_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + iq_hg : cnt + iq_hg,:),statistic_svp_csum(iq_hg : iq_hg,1:kmax),nrec_mphys,kmax)
          call writestat_nc(ncid_mphys,1,ncmname(cnt + in_cc : cnt + in_cc,:),statistic_svp_csum(in_cc : in_cc,1:kmax),nrec_mphys,kmax)
          cnt = cnt + 11

          if (cnt /= nmvar) then
            call finish(routine, 'wrong number of output variables')
          endif
        endif

        if (l_tendencies) then
          call writestat_nc(ncid_tends,1,tnctname,(/rtimee/),nrec_tends,.true.)
          cnt = 0
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_nu     :cnt + idn_cl_nu     ,:),tend_fsum(idn_cl_nu     :idn_cl_nu     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_inu    :cnt + idn_ci_inu    ,:),tend_fsum(idn_ci_inu    :idn_ci_inu    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_au     :cnt + idn_cl_au     ,:),tend_fsum(idn_cl_au     :idn_cl_au     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_au     :cnt + idq_hr_au     ,:),tend_fsum(idq_hr_au     :idq_hr_au     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_au     :cnt + idn_hr_au     ,:),tend_fsum(idn_hr_au     :idn_hr_au     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_ac     :cnt + idq_hr_ac     ,:),tend_fsum(idq_hr_ac     :idq_hr_ac     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_ac     :cnt + idn_cl_ac     ,:),tend_fsum(idn_cl_ac     :idn_cl_ac     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_br     :cnt + idn_hr_br     ,:),tend_fsum(idn_hr_br     :idn_hr_br     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_sc     :cnt + idn_hr_sc     ,:),tend_fsum(idn_hr_sc     :idn_hr_sc     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_ev     :cnt + idq_hr_ev     ,:),tend_fsum(idq_hr_ev     :idq_hr_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_ev     :cnt + idn_hr_ev     ,:),tend_fsum(idn_hr_ev     :idn_hr_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_dep    :cnt + idq_ci_dep    ,:),tend_fsum(idq_ci_dep    :idq_ci_dep    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_dep    :cnt + idq_hs_dep    ,:),tend_fsum(idq_hs_dep    :idq_hs_dep    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_dep    :cnt + idq_hg_dep    ,:),tend_fsum(idq_hg_dep    :idq_hg_dep    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_rime   :cnt + idq_ci_rime   ,:),tend_fsum(idq_ci_rime   :idq_ci_rime   ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_rime_ci:cnt + idn_cl_rime_ci,:),tend_fsum(idn_cl_rime_ci:idn_cl_rime_ci,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_rime   :cnt + idq_hs_rime   ,:),tend_fsum(idq_hs_rime   :idq_hs_rime   ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_rime_hs:cnt + idn_cl_rime_hs,:),tend_fsum(idn_cl_rime_hs:idn_cl_rime_hs,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_rime   :cnt + idq_hg_rime   ,:),tend_fsum(idq_hg_rime   :idq_hg_rime   ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_rime_hg:cnt + idn_cl_rime_hg,:),tend_fsum(idn_cl_rime_hg:idn_cl_rime_hg,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hshr_rime :cnt + idq_hshr_rime ,:),tend_fsum(idq_hshr_rime :idq_hshr_rime ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_rime_hs:cnt + idn_hr_rime_hs,:),tend_fsum(idn_hr_rime_hs:idn_hr_rime_hs,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hghr_rime :cnt + idq_hghr_rime ,:),tend_fsum(idq_hghr_rime :idq_hghr_rime ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_rime_hg:cnt + idn_hr_rime_hg,:),tend_fsum(idn_hr_rime_hg:idn_hr_rime_hg,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_rime_ri:cnt + idq_hr_rime_ri,:),tend_fsum(idq_hr_rime_ri:idq_hr_rime_ri,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_rime_ri:cnt + idq_ci_rime_ri,:),tend_fsum(idq_ci_rime_ri:idq_ci_rime_ri,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_rime_ri:cnt + idn_ci_rime_ri,:),tend_fsum(idn_ci_rime_ri:idn_ci_rime_ri,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_rime_ri:cnt + idn_hr_rime_ri,:),tend_fsum(idn_hr_rime_ri:idn_hr_rime_ri,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_col_rs :cnt + idq_hr_col_rs ,:),tend_fsum(idq_hr_col_rs :idq_hr_col_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_col_rs :cnt + idq_hs_col_rs ,:),tend_fsum(idq_hs_col_rs :idq_hs_col_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_col_rs :cnt + idn_hr_col_rs ,:),tend_fsum(idn_hr_col_rs :idn_hr_col_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_col_rs :cnt + idn_hs_col_rs ,:),tend_fsum(idn_hs_col_rs :idn_hs_col_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_col_ri :cnt + idq_hr_col_ri ,:),tend_fsum(idq_hr_col_ri :idq_hr_col_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_col_ri :cnt + idq_ci_col_ri ,:),tend_fsum(idq_ci_col_ri :idq_ci_col_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_col_ri :cnt + idn_ci_col_ri ,:),tend_fsum(idn_ci_col_ri :idn_ci_col_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_col_ri :cnt + idn_hr_col_ri ,:),tend_fsum(idn_hr_col_ri :idn_hr_col_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_cl_het    :cnt + idq_cl_het    ,:),tend_fsum(idq_cl_het    :idq_cl_het    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_het    :cnt + idn_cl_het    ,:),tend_fsum(idn_cl_het    :idn_cl_het    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_het    :cnt + idq_hr_het    ,:),tend_fsum(idq_hr_het    :idq_hr_het    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_het    :cnt + idn_hr_het    ,:),tend_fsum(idn_hr_het    :idn_hr_het    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_cl_hom    :cnt + idq_cl_hom    ,:),tend_fsum(idq_cl_hom    :idq_cl_hom    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_hom    :cnt + idn_cl_hom    ,:),tend_fsum(idn_cl_hom    :idn_cl_hom    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_col_iis:cnt + idq_ci_col_iis,:),tend_fsum(idq_ci_col_iis:idq_ci_col_iis,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_col_iis:cnt + idn_ci_col_iis,:),tend_fsum(idn_ci_col_iis:idn_ci_col_iis,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_col_sss:cnt + idn_hs_col_sss,:),tend_fsum(idn_hs_col_sss:idn_hs_col_sss,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hsci_col  :cnt + idq_hsci_col  ,:),tend_fsum(idq_hsci_col  :idq_hsci_col  ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_col_hs :cnt + idn_ci_col_hs ,:),tend_fsum(idn_ci_col_hs :idn_ci_col_hs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hghs_col  :cnt + idq_hghs_col  ,:),tend_fsum(idq_hghs_col  :idq_hghs_col  ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_col_hg :cnt + idn_hs_col_hg ,:),tend_fsum(idn_hs_col_hg :idn_hs_col_hg ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_cv     :cnt + idq_ci_cv     ,:),tend_fsum(idq_ci_cv     :idq_ci_cv     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_cv     :cnt + idn_ci_cv     ,:),tend_fsum(idn_ci_cv     :idn_ci_cv     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_cv     :cnt + idq_hs_cv     ,:),tend_fsum(idq_hs_cv     :idq_hs_cv     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_cv     :cnt + idn_hs_cv     ,:),tend_fsum(idn_hs_cv     :idn_hs_cv     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_sc     :cnt + idn_cl_sc     ,:),tend_fsum(idn_cl_sc     :idn_cl_sc     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_mul    :cnt + idn_ci_mul    ,:),tend_fsum(idn_ci_mul    :idn_ci_mul    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_mul    :cnt + idq_ci_mul    ,:),tend_fsum(idq_ci_mul    :idq_ci_mul    ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_me     :cnt + idn_ci_me     ,:),tend_fsum(idn_ci_me     :idn_ci_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_me     :cnt + idq_ci_me     ,:),tend_fsum(idq_ci_me     :idq_ci_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_me     :cnt + idn_hs_me     ,:),tend_fsum(idn_hs_me     :idn_hs_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_me     :cnt + idq_hs_me     ,:),tend_fsum(idq_hs_me     :idq_hs_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hg_me     :cnt + idn_hg_me     ,:),tend_fsum(idn_hg_me     :idn_hg_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_me     :cnt + idq_hg_me     ,:),tend_fsum(idq_hg_me     :idq_hg_me     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_ev     :cnt + idn_ci_ev     ,:),tend_fsum(idn_ci_ev     :idn_ci_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_ev     :cnt + idq_ci_ev     ,:),tend_fsum(idq_ci_ev     :idq_ci_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_ev     :cnt + idn_hs_ev     ,:),tend_fsum(idn_hs_ev     :idn_hs_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_ev     :cnt + idq_hs_ev     ,:),tend_fsum(idq_hs_ev     :idq_hs_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hg_ev     :cnt + idn_hg_ev     ,:),tend_fsum(idn_hg_ev     :idn_hg_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_ev     :cnt + idq_hg_ev     ,:),tend_fsum(idq_hg_ev     :idq_hg_ev     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_eme_ic :cnt + idn_ci_eme_ic ,:),tend_fsum(idn_ci_eme_ic :idn_ci_eme_ic ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_eme_ic :cnt + idq_ci_eme_ic ,:),tend_fsum(idq_ci_eme_ic :idq_ci_eme_ic ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_eme_ri :cnt + idn_ci_eme_ri ,:),tend_fsum(idn_ci_eme_ri :idn_ci_eme_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_eme_ri :cnt + idq_ci_eme_ri ,:),tend_fsum(idq_ci_eme_ri :idq_ci_eme_ri ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_eme_sc :cnt + idn_hs_eme_sc ,:),tend_fsum(idn_hs_eme_sc :idn_hs_eme_sc ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_eme_sc :cnt + idq_hs_eme_sc ,:),tend_fsum(idq_hs_eme_sc :idq_hs_eme_sc ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_eme_rs :cnt + idn_hs_eme_rs ,:),tend_fsum(idn_hs_eme_rs :idn_hs_eme_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_eme_rs :cnt + idq_hs_eme_rs ,:),tend_fsum(idq_hs_eme_rs :idq_hs_eme_rs ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hg_eme_gc :cnt + idn_hg_eme_gc ,:),tend_fsum(idn_hg_eme_gc :idn_hg_eme_gc ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_eme_gc :cnt + idq_hg_eme_gc ,:),tend_fsum(idq_hg_eme_gc :idq_hg_eme_gc ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hg_eme_gr :cnt + idn_hg_eme_gr ,:),tend_fsum(idn_hg_eme_gr :idn_hg_eme_gr ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_eme_gr :cnt + idq_hg_eme_gr ,:),tend_fsum(idq_hg_eme_gr :idq_hg_eme_gr ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_se     :cnt + idn_cl_se     ,:),tend_fsum(idn_cl_se     :idn_cl_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_cl_se     :cnt + idq_cl_se     ,:),tend_fsum(idq_cl_se     :idq_cl_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_ci_se     :cnt + idn_ci_se     ,:),tend_fsum(idn_ci_se     :idn_ci_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_ci_se     :cnt + idq_ci_se     ,:),tend_fsum(idq_ci_se     :idq_ci_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hr_se     :cnt + idn_hr_se     ,:),tend_fsum(idn_hr_se     :idn_hr_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hr_se     :cnt + idq_hr_se     ,:),tend_fsum(idq_hr_se     :idq_hr_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hs_se     :cnt + idn_hs_se     ,:),tend_fsum(idn_hs_se     :idn_hs_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hs_se     :cnt + idq_hs_se     ,:),tend_fsum(idq_hs_se     :idq_hs_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_hg_se     :cnt + idn_hg_se     ,:),tend_fsum(idn_hg_se     :idn_hg_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_hg_se     :cnt + idq_hg_se     ,:),tend_fsum(idq_hg_se     :idq_hg_se     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idq_cl_sa     :cnt + idq_cl_sa     ,:),tend_fsum(idq_cl_sa     :idq_cl_sa     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + idn_cl_sa     :cnt + idn_cl_sa     ,:),tend_fsum(idn_cl_sa     :idn_cl_sa     ,1:kmax),nrec_tends,kmax)
          call writestat_nc(ncid_tends,1,nctname(cnt + iret_cc       :cnt + iret_cc       ,:),tend_fsum(iret_cc       :iret_cc       ,1:kmax),nrec_tends,kmax)
          cnt = cnt + 93

          if (cnt /= ntvar) then
            call finish(routine, 'wrong number of tendencies')
          endif
        endif
      end if
    end if

    ! Reset all accumulated fields
    if (l_statistics) then
      statistic_mphys = 0.
      statistic_sv0_fsum = 0.
      statistic_sv0_count = 0.
      statistic_sv0_csum = 0.
      statistic_svp_fsum = 0.
      statistic_svp_csum = 0.
    endif
    if (l_tendencies) then
      tend_fsum = 0.
    endif
  end subroutine writebulkmicrostat3

  ! Replace non-finite values and clamp to float32 range before netCDF writes.
  subroutine sanitize_nc_2d(arr, name)
    use modmpi, only : myid
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : real32
    implicit none
    real, intent(inout) :: arr(:,:)
    character(len=*), intent(in) :: name
    real, parameter :: nc_fill = -999.
    real, parameter :: nc_abs_max = real(huge(1.0_real32), kind=kind(arr))
    integer :: n_nonfinite, n_clamped

    n_nonfinite = count(.not. ieee_is_finite(arr))
    if (n_nonfinite > 0) then
      where (.not. ieee_is_finite(arr))
        arr = nc_fill
      end where
    end if

    n_clamped = count(abs(arr) > nc_abs_max)
    if (n_clamped > 0) then
      where (abs(arr) > nc_abs_max)
        arr = sign(nc_abs_max, arr)
      end where
    end if

    if (myid == 0 .and. (n_nonfinite > 0 .or. n_clamped > 0)) then
      write(6,'(A,1X,A,1X,A,I0,1X,A,I0)') 'WARNING modbulkmicrostat3:', trim(name), &
        'nonfinite=', n_nonfinite, 'clamped=', n_clamped
    end if
  end subroutine sanitize_nc_2d

!------------------------------------------------------------------------------!
  subroutine exitbulkmicrostat3
    use modstat_nc,    only : exitstat_nc
    use modmicrodata3, only : l_statistics, l_tendencies
    use modmpi,        only : myid
    implicit none

    if (.not. lmicrostat)  return

    if (myid==0) then
       if (l_statistics) then
          call exitstat_nc(ncid_mphys)
       endif
       if (l_tendencies) then
          call exitstat_nc(ncid_tends)
       endif
    endif
  end subroutine exitbulkmicrostat3

!------------------------------------------------------------------------------!

end module modbulkmicrostat3
