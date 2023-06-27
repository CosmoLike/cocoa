!2011/05/26: changed templates to Dls

module foregrounds
  !use cmbtypes
  use settings
  use FileUtils
  implicit none

  private 

  public :: foreground_params, utindex, InitForegroundData, dl_foreground, &
       GetForegroundParamsFromArray, GetAlphaPrior, GetCirrusFactor, &
       GetRadioAmpUncPrior, GetRadioAmpPrior, InitRadioAmpPrior, &
       GetRadioClAmpUncPrior, GetRadioClAmpPrior, InitRadioClAmpPrior, &
       dl_cib_foreground,cosmo_scale_ksz,cosmo_scale_tsz,index_tsz,index_ksz, pkSZ, &
       index_dg_po,index_dg_cl,index_cirrus,index_rg_po, &
       HaveForegroundsBeenInitialized,setForegroundsUninitialized,Bnu,&
       read_dl_template,read_cl_template,nForegroundParams,&
       getForegroundPriorLnL,InitFGModel,printForegrounds,&
       ReportFGLmax, &
       OpenReadBinaryFile,OpenWriteBinaryFile,OpenReadBinaryStreamFile,calFactorsToCIB,MaxNFreq
  
!set in settings.f90  
  integer, parameter :: lmax=13500
  integer, parameter :: MaxNFreq = 6
  integer, parameter :: NDecorrel = (MaxNFreq-1)*(MaxNFreq)/2
  integer, parameter :: nForegroundParams=39+NDecorrel
  integer, dimension(NDecorrel,NDecorrel) :: dmatrix_index
  logical :: cosmological_scaling_ksz
  logical :: cosmological_scaling_tsz,ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220
  logical :: single_clustered_freq_scaling
  logical :: only1HaloTszCib, CIB_decorrelation_matrix
  logical :: tSZ_CIB_logFreq
  logical :: ShangModelCorrelationShape,applyCIBCalToCirrus,calFactorsToCIB
  integer, parameter :: ntsz=1,nksz=2,ncirrus=1
  integer, parameter:: index_tsz=1, index_ksz = index_tsz+ntsz, &
       index_dg_po=index_ksz+nksz, &
       index_dg_cl=index_dg_po+2,&
       index_cirrus=index_dg_cl+4,&
       index_rg_po = index_cirrus + ncirrus

  type foreground_params
     !all czero's are D_3000 for that signal
     double precision czero_tsz
     double precision czero_ksz
     double precision czero_ksz2
     double precision czero_dg_po
     double precision czero_dg_cl
     double precision czero_dg_cl2
     double precision czero_cirrus
     double precision czero_rg_po
     double precision czero_rg_cl

     double precision T_dg_po
     double precision beta_dg_po
     double precision sigmasq_dg_po

     double precision T_dg_cl
     double precision beta_dg_cl
     double precision sigmasq_dg_cl

     double precision T_dg_cl2
     double precision beta_dg_cl2
     double precision sigmasq_dg_cl2

     double precision alpha_rg
     double precision sigmasq_rg

     double precision T_cirrus
     double precision beta_cirrus
     
     double precision tsz_dg_cor
     double precision tsz_rg_cor

     double precision dg_cl_ell_power

     double precision tsz_cib_slope

  end type foreground_params

  double precision, parameter :: d3000 = 3000*3001/(2*pi)

  ! Kinetic and thermal SZ effect templates
  double precision, dimension(2:lmax) :: ksz_templ, tsz_templ, ksz2_templ
  double precision, dimension(2:lmax) :: clust_dg_templ,clust2_dg_templ
  double precision, dimension(2:lmax) :: cirrus_templ
  double precision, dimension(2:lmax) :: clust_rg_templ
  double precision, dimension(2:lmax) :: l_divide_3000

!  double precision,dimension(2:lmax) :: cls_diff
  
  logical :: SuccessfulInitialization 
  
  double precision DelAlphaPrior,CirrusFactorPrior,AddAlphaPoisson
  double precision radio_amp,radio_unc,radio_dl_amp,radio_dl_unc
  double precision HFIPoissonScale
contains
  subroutine printForegrounds(fgs)
    type(foreground_params),intent(in):: fgs
    
     print*,'czero_tsz',fgs%czero_tsz
     print*,'czero_ksz',fgs%czero_ksz

     print*,'czero_ksz2',fgs%czero_ksz2
     print*,'czero_dg_po',fgs%czero_dg_po
     print*,'czero_dg_cl',fgs%czero_dg_cl
     print*,'czero_dg_cl2',fgs%czero_dg_cl2
     print*,'czero_cirrus',fgs%czero_cirrus
     print*,'czero_rg_po',fgs%czero_rg_po
     print*,'czero_rg_cl',fgs%czero_rg_cl

     print*,'T_dg_po',fgs%T_dg_po
     print*,'beta_dg_po',fgs%beta_dg_po
     print*,'sigmasq_dg_po',fgs%sigmasq_dg_po
     print*,'T_dg_cl',fgs%T_dg_cl
     print*,'beta_dg_cl',fgs%beta_dg_cl
     print*,'sigmasq_dg_cl',fgs%sigmasq_dg_cl
     print*,'T_dg_cl2',fgs%T_dg_cl2
     print*,'beta_dg_cl2',fgs%beta_dg_cl2
     print*,'sigmasq_dg_cl2',fgs%sigmasq_dg_cl2
     print*,'alpha_rg',fgs%alpha_rg
     print*,'sigmasq_rg',fgs%sigmasq_rg
     print*,'T_cirrus',fgs%T_cirrus
     print*,'beta_cirrus',fgs%beta_cirrus

     print*,'tsz_dg_cor_const',fgs%tsz_dg_cor
     print*,'tsz_rg_cor',fgs%tsz_rg_cor
     print*,'dg_cl_ell_power',fgs%dg_cl_ell_power
     print*,'tsz cib slope', fgs%tsz_cib_slope

  end subroutine printForegrounds

  subroutine setForegroundsUninitialized()
    SuccessfulInitialization=.false.
  end subroutine setForegroundsUninitialized
  function HaveForegroundsBeenInitialized()
    logical HaveForegroundsBeenInitialized
    HaveForegroundsBeenInitialized = SuccessfulInitialization
  end function HaveForegroundsBeenInitialized


  !returns Dl = l*(l+1)/2pi*Cl for foregrounds
  function dl_foreground(params,ifr,jfr,nfr,eff_fr,norm_fr,req_lmin,req_lmax,component_spectra)
    type(foreground_params) :: params
    integer :: i,ifr,jfr,nfr
    integer, intent(in) :: req_lmin,req_lmax
    double precision, dimension(5,nfr) :: eff_fr
    double precision, dimension(5) :: norm_fr
    double precision,dimension(req_lmin:req_lmax) :: dl_foreground
    double precision :: fri,frj,norm,fr0,frqdep
    double precision,dimension(2:lmax) :: dl_dg_po, dl_dg_cl, dl_rg_po, dl_k_sz, &
         dl_t_sz, dl_cirrus
    double precision,dimension(2:lmax) :: dl_tsz_rg_cor, dl_tsz_dgcl_cor
    double precision,dimension(2:lmax) :: dl_dg, dl_rg
    double precision,dimension(2:lmax) :: dltmpi,dltmpj
    double precision,dimension(req_lmax,7), optional, intent(out) :: component_spectra
    
    if (req_lmax .gt. lmax) then 
       print*,'asked for too large ell in dl_foreground'
       call mpistop()
    endif
    if (req_lmin .lt. 2) then 
       print*,'asked for too small ell in dl_foreground'
       call mpistop()
    endif
    
    
    dl_dg_cl = dl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,jfr),norm_fr(1),ifr,jfr)
    dl_dg_po = dl_dusty_poisson(params,eff_fr(2,ifr),eff_fr(2,jfr),norm_fr(2),ifr,jfr)
    dl_rg_po = dl_radio(params,eff_fr(3,ifr),eff_fr(3,jfr),norm_fr(3))
    dl_k_sz = dl_ksz(params)
    dl_t_sz = dl_tsz(params,eff_fr(5,ifr),eff_fr(5,jfr),norm_fr(5))
    !TSZ-Dusty correlation

    if (only1HaloTszCib) then
       dltmpi = dl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr,.true.) 
       dltmpj = dl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr,.true.) 
    else
       dltmpi = dl_dusty_clustered(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1),ifr,jfr) + &
            dl_dusty_poisson(params,eff_fr(2,ifr),eff_fr(2,ifr),norm_fr(2),ifr,jfr) 
       dltmpj = dl_dusty_clustered(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1),ifr,jfr) + &
            dl_dusty_poisson(params,eff_fr(2,jfr),eff_fr(2,jfr),norm_fr(2),ifr,jfr) 
    endif
    dl_tsz_dgcl_cor = -1 * tsz_dgcl_cor(params) * &
         ( tsz_cib(params,eff_fr(2,jfr)) * sqrt( dl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * dltmpj) + &
         tsz_cib(params,eff_fr(2,ifr)) * sqrt( dl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * dltmpi))
    

    !TSZ-Radio correlation
    if (params%tsz_rg_cor .eq. 0) then
       dl_tsz_rg_cor = 0.0
    else
       dl_tsz_rg_cor = -1 * params%tsz_rg_cor * tsz_rg_cor() * ( &
            sqrt(dl_tsz(params,eff_fr(5,ifr),eff_fr(5,ifr),norm_fr(5)) * &
            dl_radio(params,eff_fr(1,jfr),eff_fr(1,jfr),norm_fr(1))) &
            + sqrt(dl_tsz(params,eff_fr(5,jfr),eff_fr(5,jfr),norm_fr(5)) * &
            dl_radio(params,eff_fr(1,ifr),eff_fr(1,ifr),norm_fr(1))) )
    endif
    
    dl_cirrus = dl_galcirrus(params,eff_fr(1,ifr),eff_fr(1,jfr))
    
    dl_dg = dl_dg_po + dl_dg_cl
    dl_rg = dl_rg_po
    
    if (present(component_spectra)) then 

       component_spectra(2:req_lmax,1) = dl_dg_po(2:req_lmax)
       component_spectra(2:req_lmax,2) = dl_dg_cl(2:req_lmax)
       component_spectra(2:req_lmax,3) = dl_k_sz(2:req_lmax)
       component_spectra(2:req_lmax,4) = dl_t_sz(2:req_lmax)
       component_spectra(2:req_lmax,5) = dl_rg(2:req_lmax)
       component_spectra(2:req_lmax,6) = dl_tsz_dgcl_cor(2:req_lmax)
       component_spectra(2:req_lmax,7) = dl_cirrus(2:req_lmax) + dl_tsz_rg_cor(2:req_lmax)
    endif
    dltmpi = (dl_dg + dl_rg + dl_t_sz + dl_k_sz + dl_cirrus + dl_tsz_dgcl_cor + dl_tsz_rg_cor)
    dl_foreground = dltmpi(req_lmin:req_lmax)
  end function dl_foreground

  function tsz_cib(params,freq)
    double precision tsz_cib
    type(foreground_params) :: params
    double precision :: freq
    tsz_cib =  params%tsz_dg_cor

  end function tsz_cib



  function dl_cib_foreground(params,eff_fr,norm_fr,ifr,jfr)
    type(foreground_params) :: params
    integer,intent(in) :: ifr,jfr
    double precision :: eff_fr
    double precision :: norm_fr
    double precision,dimension(2:lmax) :: dl_cib_foreground
    double precision,dimension(2:lmax) :: dl_dg_po, dl_dg_cl
    
    dl_dg_cl = dl_dusty_clustered(params,eff_fr,eff_fr,norm_fr,ifr,jfr)
    dl_dg_po = dl_dusty_poisson(params,eff_fr,eff_fr,norm_fr,ifr,jfr)
    
    dl_cib_foreground = dl_dg_po + dl_dg_cl 

  end function dl_cib_foreground
  
  !updated to Dl & template being normalized
  function dl_radio(params,fri,frj,fr0)
    double precision,dimension(2:lmax)  :: dl_radio
    type(foreground_params) :: params
    double precision :: fri,frj,fr0
   
    dl_radio = (params%czero_rg_po/dBdT(fri,fr0)/dBdT(frj,fr0)*&
         (fri/fr0*frj/fr0)**(params%alpha_rg + &
         log(fri/fr0*frj/fr0)/2 * params%sigmasq_rg**2)) * &
         l_divide_3000**2

    if (params%czero_rg_cl .gt. 0) then 
       dl_radio(2:lmax) = dl_radio(2:lmax) + &
            ( params%czero_rg_cl/dBdT(fri,fr0) / &
            dBdT(frj,fr0)*&
            (fri/fr0*frj/fr0)**(params%alpha_rg) ) * &
            clust_rg_templ
    endif
  end function dl_radio

  function cirrus_power3000(params,fri,frj)
    double precision  :: cirrus_power3000
    type(foreground_params) :: params
    double precision :: fri, frj, fr0
    double precision :: frqdep
    
    fr0=220.0
    frqdep = ((fri*frj)/(fr0*fr0))**(params%beta_cirrus)
    frqdep = frqdep *&
         Bnu(fri,fr0,params%T_cirrus)*Bnu(frj,fr0,params%T_cirrus)
    
    frqdep = frqdep /( dBdT(fri,fr0) * dBdT(frj,fr0) )
    cirrus_power3000 = params%czero_cirrus*frqdep
  end function cirrus_power3000

  !updated to Dl & template being normalized
  function dl_galcirrus(params,fri,frj)
    double precision,dimension(2:lmax)  :: dl_galcirrus
    type(foreground_params) :: params
    double precision :: fri, frj, fr0
    double precision :: power

    power = cirrus_power3000(params,fri,frj)
    
    dl_galcirrus=power* cirrus_templ

  end function dl_galcirrus

  !turned to Dl and normalized templates
  function dl_dusty_poisson(params,fri,frj,fr0,ifr,jfr)
    double precision,dimension(2:lmax)  :: dl_dusty_poisson
    type(foreground_params) :: params
    integer, intent(in) :: ifr,jfr
    double precision :: fri, frj, fr0
    double precision :: frqdep
    double precision :: ff,ff1,ff2,decor
    double precision :: effalpha,effsigmasq,effalpha_2

    !basic powerlaw
    frqdep = ((fri*frj)/(fr0*fr0))**( params%beta_dg_po)

    !not <1 like decor in the other case...this boosts the autospectra
    decor =  ((fri*frj)/(fr0*fr0))**( &
         log(fri/fr0*frj/fr0)/2 * params%sigmasq_dg_po )

    frqdep = frqdep*decor
       
    frqdep = frqdep *&
         Bnu(fri,fr0,params%T_dg_po)*Bnu(frj,fr0,params%T_dg_po)
    
    dl_dusty_poisson = (params%czero_dg_po/dBdT(fri,fr0)/dBdT(frj,fr0)*frqdep) * &
         l_divide_3000**2


  end function dl_dusty_poisson

  !turned to Dl and normalized templates
  function dl_dusty_clustered(params,fri,frj,fr0,ifr,jfr,only1halo)
    integer, intent(in) :: ifr,jfr
    double precision,dimension(2:lmax)  :: dl_dusty_clustered
    type(foreground_params) :: params
    double precision :: fri,frj,fr0,frqdep,frqdep0,effalpha_cl, effalpha_cl_2,frqdeplin,effsigmasq_cl
    double precision :: effalpha,effsigmasq,effalpha_2,decor,ff,ff1,ff2
    logical, intent(in), optional :: only1halo
    logical :: Want2Halo

    if (present(only1halo)) then
        Want2Halo = .not. only1halo
    else
        Want2Halo = .true.
    end if

    effalpha_cl = params%T_dg_cl + AddAlphaPoisson*params%T_dg_po
    effsigmasq_cl = params%sigmasq_dg_cl + AddAlphaPoisson*params%sigmasq_dg_po
    effalpha_cl_2 = params%beta_dg_cl + AddAlphaPoisson*params%beta_dg_po
    frqdep = ((fri*frj)/(fr0*fr0))**(effalpha_cl_2)

    frqdep = frqdep *&
         Bnu(fri,fr0,effalpha_cl)*Bnu(frj,fr0,effalpha_cl)
    
    frqdep = frqdep / dBdT(fri,fr0) / dBdT(frj,fr0)
    

    !not <1 like decor in the other case...this boosts the autospectra
    decor =  ((fri*frj)/(fr0*fr0))**( &
         log(fri/fr0*frj/fr0)/2 * effsigmasq_cl )

    frqdep = frqdep*decor


    dl_dusty_clustered = clust_dg_templ * (params%czero_dg_cl*frqdep)
    if (params%dg_cl_ell_power /= 0) &
         dl_dusty_clustered =  dl_dusty_clustered * (l_divide_3000)**params%dg_cl_ell_power
    
    if (Want2Halo .and. params%czero_dg_cl2 /= 0) then
       if (single_clustered_freq_scaling ) then 
          dl_dusty_clustered = dl_dusty_clustered+ clust2_dg_templ*params%czero_dg_cl2*frqdep 
       else
          effalpha_cl = params%T_dg_cl2 + AddAlphaPoisson*params%T_dg_po
          effsigmasq_cl = params%sigmasq_dg_cl2 + AddAlphaPoisson*params%sigmasq_dg_po
          effalpha_cl_2 = params%beta_dg_cl2 + AddAlphaPoisson*params%beta_dg_po
          frqdep = ((fri*frj)/(fr0*fr0))**(effalpha_cl_2)
          frqdep = frqdep *&
               Bnu(fri,fr0,effalpha_cl)*Bnu(frj,fr0,effalpha_cl)
          frqdep = frqdep / dBdT(fri,fr0) / dBdT(frj,fr0)
          decor =  ((fri*frj)/(fr0*fr0))**( &
               log(fri/fr0*frj/fr0)/2 * effsigmasq_cl )
          frqdep = frqdep*decor
          
          dl_dusty_clustered = dl_dusty_clustered+ clust2_dg_templ*(params%czero_dg_cl2*frqdep )
          
       endif
       
    endif
  end function dl_dusty_clustered

!updated to Dl & template being normalized
  function dl_tsz(params,fri,frj,fr0)
    double precision,dimension(2:lmax) :: dl_tsz
    type(foreground_params) :: params
    double precision :: fri,frj,fr0

    dl_tsz = (params%czero_tsz * tszFreqDep(fri,fr0) * tszFreqDep(frj,fr0) ) * tsz_templ

  end function dl_tsz

!function to do a quick EoR version
function pkSZ(dZ)
  double precision :: pkSZ
  double precision :: dZ

  pkSZ = 1.445 * (dZ/1.05)** 0.51

end function pkSZ

!based on Laurie's email and templates!
!apply cosmological scaling
  function cosmo_scale_ksz(H0,sigma8,omegab,omegam,ns,tau)
    double precision :: H0,sigma8,omegab,omegam,ns,tau
    double precision :: cosmo_scale_ksz
    if (cosmological_scaling_ksz) then
       cosmo_scale_ksz = ((H0/71.0)**1.7 ) &
            * ( (sigma8/.8)**4.7 ) &
            * ( (omegab/.044)**2.1 ) &
            * ( (omegam/.264)**(-0.44) ) &
            * ( (ns/.96)**(-0.19) ) !&

    else
       cosmo_scale_ksz = 1.0
    endif
  end function cosmo_scale_ksz

!based on Laurie's email and templates!
!apply cosmological scaling
  function cosmo_scale_tsz(H0,sigma8,omegab)
    double precision :: H0,sigma8,omegab
    double precision :: cosmo_scale_tsz
    if (cosmological_scaling_tsz) then
       cosmo_scale_tsz = ((H0/71.0)**1.73 ) &
            * ( (sigma8/.8)**8.34 ) &
            * ( (omegab/.044)**2.81 ) 
    else
       cosmo_scale_tsz = 1.0
    endif
  end function cosmo_scale_tsz

  !updated to Dl & template being normalized
  function dl_ksz(params)
    double precision,dimension(2:lmax) :: dl_ksz
    type(foreground_params) :: params

    dl_ksz = params%czero_ksz *ksz_templ
    if (params%czero_ksz2 /= 0) &
         dl_ksz = dl_ksz + params%czero_ksz2 *ksz2_templ

  end function dl_ksz

  function flat_tsz_cor()
    logical flat_tsz_cor
    flat_tsz_cor = (.not. ShangModelCorrelationShape)
  end function flat_tsz_cor
  
  ! The correlation fraction between tSZ and dusty sources, normalized to 1 at high ell
  function tsz_dgcl_cor(params)
    double precision,dimension(2:lmax) :: tsz_dgcl_cor
    type(foreground_params) :: params
    
    !motivated by shaw analysis of Sehgal sims
!    tsz_dgcl_cor = max(0., (.3 - .2 * exp(-(l-500.)/1000.))/.3)
    !motivated by simplicity
    
    
    if (ShangModelCorrelationShape) then
       tsz_dgcl_cor = -0.0703 * (l_divide_3000*l_divide_3000) + &
            0.612 * l_divide_3000 + &
            0.458

    else
       tsz_dgcl_cor = 1.0
       if (params%tsz_cib_slope /= 0) then
          tsz_dgcl_cor = tsz_dgcl_cor + (l_divide_3000-1d0)*params%tsz_cib_slope
       endif
    endif
    
  end function tsz_dgcl_cor

  ! The correlation fraction between tSZ and dusty sources, normalized to 1 at high ell
  function tsz_rg_cor()

    double precision :: tsz_rg_cor

    tsz_rg_cor = 1.0

  end function tsz_rg_cor

  !
  ! nu,nu0 in GHz
  !
  ! dBdT is proportional to derivative of planck function
  ! but is normalized so its equal to 1 at nu0
  !
  function dBdT(nu,nu0)

    double precision x, x0, dBdT, dBdT0, nu, nu0

    x0 = nu0/56.78
    dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2

    x = nu/56.78
    dBdT = x**4 * exp(x) / (exp(x)-1)**2 / dbdT0

  end function dBdT
  
  !proportional to the Planck function normalized to 1 at nu0
  function Bnu(nu,nu0,T)
    double precision Bnu, nu,nu0,T
    !h/k
    !4.799237 × 10-11 s K
    !expect GHz
    ! so 4.799237e-2 K/GHz
    double precision, parameter :: hk = 4.799237e-2
    
    Bnu = (nu/nu0)**3
    Bnu = Bnu * (exp( hk*nu0/T)-1d0) / (exp( hk*nu/T)-1d0) 
    
  end function Bnu

  

  
  !
  ! nu, nu0 in GHz
  ! Gives the tsz frequency dependence normalized so its 1 at nu0
  !
  function tszFreqDep(nu,nu0)

    double precision :: tszFreqDep, tszFreqDep0, nu, nu0, x, x0

    x = nu / 56.78
    x0 = nu0 / 56.78

    tszFreqDep0 = x0*(exp(x0)+1)/(exp(x0)-1) - 4
    tszFreqDep = x*(exp(x)+1)/(exp(x)-1) - 4
    tszFreqDep = tszFreqDep/tszFreqDep0

  end function tszFreqDep
  
!read template in Dl form from file
!since templates are in Dl's internally, do nothing to it
  function read_dl_template(filename)
    character*(*), intent(in) :: filename
    double precision, dimension(2:lmax) :: read_dl_template
    double precision :: realtmp
    integer :: ll 
    Type(TTextFile) :: F
    logical wexist
    read_dl_template(2:lmax)=0.0

    if (filename/='') then
       inquire(FILE=trim(filename),EXIST=wexist)
       if (.not. wexist) then
          print*,'SPT hiell 2020, missing template file:', trim(filename)
          call mpistop()
       endif

       call F%Open(filename)
       do
          read(F%unit,*,end=2) ll, realtmp
          if (ll>=2 .and. ll<=lmax) &
               read_dl_template(ll) = realtmp
       end do
2      call F%Close()
    end if
  end function read_dl_template

!read template in Cl form from file, convert to Dl
  function read_cl_template(filename)
    character*(*), intent(in) :: filename
    double precision, dimension(2:lmax) :: read_cl_template
    double precision :: realtmp
    integer :: ll
    Type(TTextFile) :: F
    logical wexist
    read_cl_template(2:lmax)=0.0
    if (filename/='') then
       inquire(FILE=trim(filename),EXIST=wexist)
       if (.not. wexist) then
          print*,'SPT hiell 2020, missing template file:', trim(filename)
          call mpistop()
       endif
       call F%Open(filename)
       do
          read(F%unit,*,end=2) ll, realtmp
          if (ll>=2 .and. ll<=lmax) &
               read_cl_template(ll) = realtmp * (ll*(ll+1)/2d0/pi)
       end do
2      call F%Close()
    end if
  end function read_cl_template


  !
  ! Loads the clustered, ksz, and tsz templates which are expected to be in
  ! foreground_folder with names cluster_*.dat, ksz.dat, and tsz.dat
  !
  subroutine InitForegroundData(fClustered,fClustered2,fKSZ,fKSZ2,fTSZ,&
       del_alpha,relative_alpha_cluster)
    integer :: l,dum,i
    character(len=120) :: file
    character*(*) :: fClustered,fClustered2, fKSZ, fKSZ2, fTSZ
    double precision :: del_alpha
    logical :: relative_alpha_cluster
    SuccessfulInitialization=.true.
    
    ! clustered template
    clust_dg_templ = read_dl_template(fClustered)
    clust2_dg_templ = read_dl_template(fClustered2)
    clust_dg_templ  =  clust_dg_templ  / clust_dg_templ(3000)
    clust2_dg_templ =   clust2_dg_templ / clust2_dg_templ(3000)   

    ! KSZ template
    ksz_templ  = read_dl_template(fKSZ)
    ksz_templ  = ksz_templ / ksz_templ(3000)
    
    ! 2nd KSZ template
    ksz2_templ  = read_dl_template(fKSZ2)
    ksz2_templ  = ksz2_templ / ksz2_templ(3000)
    
    ! TSZ template
    tsz_templ  = read_dl_template(fTSZ)
    tsz_templ  = tsz_templ / tsz_templ(3000)
    
    if (MpiRank == 0) then
       print*,'1halo (or all clustered) DG template:',trim(fClustered)
       print*,'2halo (2nd) DG template:',trim(fClustered2)
       print*,'kSZ template:',trim(fKSZ)
       print*,'ksz template2:',trim(fKSZ2)
       print*,'tSZ template:',trim(fTSZ)
    endif
    !prior on separation between clsutered and Poisson terms
    ! ignored if <= 0
    DelAlphaPrior = del_alpha
    AddAlphaPoisson = 0
    if (relative_alpha_cluster) then
       AddAlphaPoisson = 1
    end if

    do i=2,lmax
       l_divide_3000(i) = real(i)/3000d0
    enddo

    cirrus_templ(200:lmax) = l_divide_3000(200:lmax)**(-1.2)
    cirrus_templ(2:199)=0.0
    
    clust_rg_templ = l_divide_3000**(0.6)
    clust_rg_templ(2:49)=0.0

  end subroutine InitForegroundData
  
  subroutine InitRadioAmpPrior(iradio_amp,iradio_unc)
    double precision::iradio_amp, iradio_unc
    radio_amp=iradio_amp
    radio_unc=iradio_unc
  end subroutine InitRadioAmpPrior
  
  function GetRadioAmpPrior()
    double precision :: GetRadioAmpPrior
    GetRadioAmpPrior=radio_amp
  end function GetRadioAmpPrior
  function GetRadioAmpUncPrior()
    double precision :: GetRadioAmpUncPrior
    GetRadioAmpUncPrior=radio_unc
  end function GetRadioAmpUncPrior

  subroutine InitRadioClAmpPrior(iradio_amp,iradio_unc)
    double precision::iradio_amp, iradio_unc
    radio_dl_amp=iradio_amp
    radio_dl_unc=iradio_unc
  end subroutine InitRadioClAmpPrior
  
  function GetRadioClAmpPrior()
    double precision :: GetRadioClAmpPrior
    GetRadioClAmpPrior=radio_dl_amp
  end function GetRadioClAmpPrior
  function GetRadioClAmpUncPrior()
    double precision :: GetRadioClAmpUncPrior
    GetRadioClAmpUncPrior=radio_dl_unc
  end function GetRadioClAmpUncPrior
  
  !
  ! returns alpha prior on Poisson-cluster index
  !
  function GetAlphaPrior()
    double precision :: GetAlphaPrior
    GetAlphaPrior = DelAlphaPrior
  end function GetAlphaPrior

  !
  ! returns cirrus prior pre-factor
  ! <0 is no prior.
  function GetCirrusFactor()
    double precision :: GetCirrusFactor
    GetCirrusFactor = CirrusFactorPrior
  end function GetCirrusFactor


  !
  ! Index of the (i,j) entries of the upper triangular part of an n-by-n matrix
  !
  recursive function utindex(i,j,n)

    integer :: i, j, n, utindex

    if (i <= j) then
       utindex = (i-1)*n+j - i*(i-1)/2
    else
       utindex = utindex(j,i,n)
    end if

  end function utindex



  ! Read parameters from an array of reals ordered in the same order as the type
  ! Basically just cast the array as the type and the parameters line up
  function GetForegroundParamsFromArray(array)

    double precision, dimension(:) :: array
    type(foreground_params) :: GetForegroundParamsFromArray

    GetForegroundParamsFromArray = transfer(array,GetForegroundParamsFromArray)

  end function GetForegroundParamsFromArray

  !calculate external foreground prior LnL
  function  getForegroundPriorLnL(foregrounds)
    double precision ::  getForegroundPriorLnL
    type(foreground_params) :: foregrounds
    double precision :: del_alpha
    integer i,j
    double precision :: radio_amp,radio_unc,fradio,radio_cl_amp,radio_cl_unc
    double precision :: cirrus90, cirrus150,cirrus220
    double precision :: prior90, prior150, prior220
    double precision :: cirrus_factor
    double precision :: freq
    integer anyhigh
    getForegroundPriorLnL=0
    
    del_alpha = GetAlphaPrior()
    ! dusty index prior
    if ( del_alpha > 0) then
       getForegroundPriorLnL = getForegroundPriorLnL + (foregrounds%T_dg_po-foregrounds%T_dg_cl)**2/(2*(del_alpha)**2)
!       print*,'dusty index prior:',getForegroundPriorLnL
    end if
    
    !cirrus prior
    !updated 2011/05/28 to RKs 08+09 values
!9/22 updated with adhoc 23h value
!also decided to leave out 90 GHz prior since effectively 0 and not well matched to SED shape.
    !9/23 updated to make all three optional. 150/220 will still be default if keyword's aren't in ini file
    cirrus_factor=1 !had been 1/3, but can't tell why
    if (feedback > 3)    print*,'do I do cirrus?',ApplyCirrusPrior90, ApplyCirrusPrior150, ApplyCirrusPrior220
    if (foregrounds%czero_cirrus .ne. 0) then
       if (ApplyCirrusPrior90) then 
          freq=97.9
          cirrus90=cirrus_power3000(foregrounds,freq,freq)
          
          prior90=0.16*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (cirrus90-prior90)**2 / (2* (.06)**2)
          !               (log(cirrus90 / prior90))**2 / (2*(0.43)**2)
              if (feedback > 3)  print*,'post cirrus90 prior:',getForegroundPriorLnL       
       endif
       if (ApplyCirrusPrior150) then 
          freq=153.44
          cirrus150=cirrus_power3000(foregrounds,freq,freq)
          
          prior150=0.21*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (cirrus150-prior150)**2 / (2*(.06)**2)
           !    (log(cirrus150/prior150))**2 / (2*(0.564)**2)
          if (feedback > 3)           print*,'post cirrus150 prior:',getForegroundPriorLnL,cirrus150,prior150
       endif
       if (ApplyCirrusPrior220) then 
          freq=219.67
          cirrus220=cirrus_power3000(foregrounds,freq,freq)
          
          prior220=2.19*cirrus_factor
          getForegroundPriorLnL = getForegroundPriorLnL + &
               (cirrus220-prior220)**2 / (2*(.7)**2)
!               (log(cirrus220/prior220))**2 / (2*(0.33)**2)
          if (feedback > 3) print*,'post cirrus220 prior:',getForegroundPriorLnL,cirrus220,prior220
       endif
    endif
       
    radio_amp = GetRadioAmpPrior()
    radio_unc = GetRadioAmpUncPrior()
    if (radio_amp >= 0 .and. radio_unc > 0) then
       fradio =(foregrounds%czero_rg_po - radio_amp) / radio_unc
       getForegroundPriorLnL = getForegroundPriorLnL + fradio**2 / 2
       if (foregrounds%czero_rg_po < 0) then
          getForegroundPriorLnL = getForegroundPriorLnL + 1e6
       endif
       if (feedback > 3) print*,'post rg poisson prior:',getForegroundPriorLnL       
    endif
    
    radio_cl_amp = GetRadioClAmpPrior()
    radio_cl_unc = GetRadioClAmpUncPrior()
    if (radio_dl_amp >= 0 .and. radio_dl_unc > 0) then
       fradio =(foregrounds%czero_rg_cl - radio_cl_amp) / radio_cl_unc
       getForegroundPriorLnL = getForegroundPriorLnL + fradio**2 / 2
       if (foregrounds%czero_rg_cl < 0) then
          getForegroundPriorLnL = getForegroundPriorLnL + 1e6
       endif
       if (feedback > 3) print*,'post rg clus prior:',getForegroundPriorLnL       
    endif
    
  end function getForegroundPriorLnL
  
  
  !subroutine InitLensingTemplate(LensingTemplate)
  !  character*(*),intent(in) :: LensingTemplate
  !  
  !  cls_diff = read_dl_template(LensingTemplate)
  !end subroutine InitLensingTemplate


  subroutine InitFGModel(Ini)
    use IniFile
    use IniObjects
    use settings
    implicit none

    class(TSettingIni) :: Ini
    character(LEN=Ini_max_string_len) SPTtSZTemplate, SPTkSZTemplate, &
         SPTkSZ2Template, SPTClus2Template, &
         SPTClusTemplate!, LensingTemplate
    logical relative_alpha_cluster
    integer i,j,k

!         ShangModelCorrelationShape
!    logical ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220
!    logical cosmological_scaling_tsz, cosmological_scaling_ksz
    double precision   radio_amp, radio_unc, del_alpha

    if (SuccessfulInitialization) then 
       return
       !skip - already called successfully
    endif

    !LensingTemplate = Ini_Read_String('add_lensing_template', .false.)
    !call InitLensingTemplate(LensingTemplate)
        
    SPTtSZTemplate = Ini%Read_String_Default('spt_dataset_tSZ', '')
    SPTkSZTemplate = Ini%Read_String_Default('spt_dataset_kSZ', '')
    SPTkSZ2Template = Ini%Read_String_Default('spt_dataset_kSZ2', '')
    SPTClusTemplate = Ini%Read_String_Default('spt_dataset_clustered', '')
    SPTClus2Template = Ini%Read_String_Default('spt_dataset_clustered2', '')
        
    del_alpha = Ini%Read_Real('spt_prior_clusterpoisson',-1.0)
    radio_amp = Ini%Read_Real('radio_ampl_mean',-1.0)
    radio_unc = Ini%Read_Real('radio_ampl_unc',-1.0)
    call InitRadioAmpPrior(radio_amp,radio_unc)
    if (MpiRank == 0) &
         print*,'priors, dAlpha,radio amp, radio sigma',&
         del_alpha,radio_amp,radio_unc
    
     only1HaloTszCib = Ini%Read_Logical('only_1halo_tsz_cib',.false.)

    ShangModelCorrelationShape = Ini%Read_Logical(&
         'shang_model_correlation_shape',.false.)
   

    if (MPIRank == 0) &
         print*,'tSZ-CIB assumptions:',only1HaloTszCib, &
         ShangModelCorrelationShape
    
    single_clustered_freq_scaling = Ini%Read_Logical(&
         'single_clustered_freq_scaling',.true.)
    relative_alpha_cluster = Ini%Read_Logical('relative_alpha_cluster',.true.)

    if (MPIRank == 0) &
         print*,'CIB freq dep. assumptions:',&
         single_clustered_freq_scaling ,relative_alpha_cluster
    
    cosmological_scaling_ksz = Ini%Read_Logical(&
         'cosmological_scaling_ksz',.false.)
    cosmological_scaling_tsz = Ini%Read_Logical(&
         'cosmological_scaling_tsz',.false.)
    if (MPIRank == 0) &
         print*,'Cosmo scaling ksz, tsz:',&
         cosmological_scaling_ksz,cosmological_scaling_tsz
    
    ApplyCirrusPrior90 = Ini%Read_Logical('apply_prior_cirrus_90ghz',.false.)
    ApplyCirrusPrior150 = Ini%Read_Logical('apply_prior_cirrus_150ghz',.true.)
    ApplyCirrusPrior220 = Ini%Read_Logical('apply_prior_cirrus_220ghz',.true.)
    if (MPIRank == 0) &
         print*,'Cirrus priors used:',&
         ApplyCirrusPrior90,ApplyCirrusPrior150,ApplyCirrusPrior220    
    
    if (SPTtSZTemplate/='' .AND. SPTkSZTemplate/='' .AND. SPTkSZ2Template/='' &
         .AND. SPTClusTemplate/='') then
       call InitForegroundData(SPTClusTemplate, SPTClus2Template, &
            SPTkSZTemplate,SPTkSZ2Template,SPTtSZTemplate,del_alpha, &
            relative_alpha_cluster)
    else
       print*,'called initFGModel w/o all args'
       print*,SPTtSZTemplate,SPTkSZTemplate,SPTkSZ2Template,SPTClusTemplate
       call mpistop()
       call   setForegroundsUninitialized()
    end if
  end subroutine InitFGModel

  function ReportFGLmax()
    integer :: ReportFGLmax
    ReportFGLmax = lmax
  end function ReportFGLmax

  subroutine OpenReadBinaryFile(aname,aunit,record_length)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    integer*8,intent(in) :: record_length
    open(unit=aunit,file=aname,form='unformatted',access='direct',recl=record_length,  err=500)
    return
    
500 call MpiStop('File not found: '//trim(aname))
  end subroutine openReadBinaryFile

  subroutine OpenReadBinaryStreamFile(aname,aunit)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    open(unit=aunit,file=aname,form='unformatted',access='stream', err=500)
    return
    
500 call MpiStop('File not found: '//trim(aname))
  end subroutine OpenReadBinaryStreamFile

  subroutine OpenWriteBinaryFile(aname,aunit,record_length)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: aunit
    integer*8,intent(in) :: record_length
    open(unit=aunit,file=aname,form='unformatted',status='replace',access='direct',recl=record_length, err=500)
    return
    
500 call MpiStop('File not able to be written to: '//trim(aname))
  end subroutine OpenWriteBinaryFile



end module
