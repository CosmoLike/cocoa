module egfs
  use keysvalues
  implicit none
  private

  public :: initialize, get_cib_clustering, get_cib_poisson, get_radio_poisson, get_tsz, get_ksz, get_egfs, get_egfs_component, NUM_FG_PARAMS, FG_NAMES,get_nfr,initialize_frommem,set_cib_decor_clust,set_cib_decor_poisson,get_next_instance,free_instance

  ! The total number of foreground parameters in
  !  1 - cib-clustered
  !  2 - cib-poisson
  !  3 - radio poisson
  !  4 - tsz
  !  5 - ksz
  !  6 - total
  !
  ! The dR matrix is nfr*nfr*NUM_FG_PARAMS at each ell
  !
  integer, parameter, dimension(6) :: NUM_FG_PARAMS = (/3,3,4,3,3,16/)

  character(*), dimension(6), parameter :: FG_NAMES = (/'cib-clustering','cib-poisson   ','radio-poison  ','tsz           ','ksz           ','total         '/)

  real(8), parameter :: pi = 3.1415926579
  real(8), parameter :: d3000 = 3000*3001/(2*pi)

  type ValIdx

    real v
    integer i

  end type

  ! Holds the templates and default parameters for each instance
  type foreground_params

    type(ValIdx) alpha_dg_cl
    type(ValIdx) tilt_dg_cl
    type(ValIdx) alpha_dg_po
    type(ValIdx) sigma_dg_po
    type(ValIdx) norm_dg_cl
    type(ValIdx) norm_dg_po
    type(ValIdx) fpol_dg_po
    type(ValIdx) alpha_rg
    type(ValIdx) sigma_rg
    type(ValIdx) norm_rg
    type(ValIdx) gamma_rg
    type(ValIdx) fpol_rg
    type(ValIdx) tsz_pca1
    type(ValIdx) tsz_pca2
    type(ValIdx) norm_ov
    type(ValIdx) norm_patchy
    type(ValIdx) shift_patchy
    type(ValIdx) tsz_dgcl_cor
    type(ValIdx) tsz_mean_scale

    ! Clustered template
    real(8), dimension(2:10000) :: clust_dg_templ

    ! Patchy template
    real(8), dimension(2:10000) :: ksz_patchy_templ

    ! OV template, {C_ell, n_s, omega_b, omega_dm, sigma_8, tau}
    real(8), dimension(2:10000,6) :: ksz_ov_templ

    ! TSZ template {FID, PCA_1, PCA_2, n_s, omega_b, omega_M, sigma_8}
    real(8), dimension(2:10000,7) :: tsz_templ

    ! The number of frequencies
    integer nfr

    ! The effective band centers
    real(8), allocatable, dimension(:,:) :: eff_fr

    ! Normalization frequencies
    real(8), dimension(4) :: norm_fr

    ! Radio flux-cut and normaliztion flux cut
    real(8) :: rg_flux_cut, norm_rg_flux_cut
		
		! cib decorrelation matrix
		real(8),allocatable, dimension(:,:) :: cib_decor_clust,cib_decor_poisson
  		 
  end type


  ! All the different instances
  type(foreground_params), dimension(100), target :: foreground_instances
  integer,dimension(100) :: avail_instance
  integer::m_instance=0
		
contains

	subroutine get_next_instance(instance_id)
		integer,intent(out)::instance_id
		integer::i
		if (m_instance==0) then
			avail_instance=0
			m_instance=1
		endif
		do i=1,100
			if (avail_instance(i)==0) then
				avail_instance(i) = 1
				instance_id=i
				return
			endif
		enddo
		instance_id = 0
	end subroutine get_next_instance
	
	subroutine free_instance(instance_id)
		integer,intent(in)::instance_id
		type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
		if (allocated(this%eff_fr)) deallocate(this%eff_fr)
		if (allocated(this%cib_decor_clust)) deallocate(this%cib_decor_clust)
		if (allocated(this%cib_decor_poisson)) deallocate(this%cib_decor_poisson)
		avail_instance(instance_id)=0
	end subroutine free_instance
  !
  ! Load foreground templates for a given instance
  ! and intialize the default parameters
  !

	subroutine initialize_internal(instance_id,keys,values,nkv,error)
		integer,intent(inout) :: instance_id, error
    integer,intent(in) :: nkv
    character(*), dimension(nkv),intent(in) :: keys
    character(*), dimension(nkv),intent(in) :: values
    type(foreground_params), pointer :: this
    real(8)::ns_cib
    integer::ell

    this => foreground_instances(instance_id)
  
    if (hasParam(keys,values,"ns_cib_clustering")) then
      ns_cib = getRealParam(keys,values,"ns_cib_clustering")
      do ell=2,10000
        this%clust_dg_templ(ell) = (1.*ell)**ns_cib  
      end do
    endif
    call NormalizeTemplate(this%clust_dg_templ)
      
    this%nfr=getIntParam(keys,values,"nfr")
    
    if (allocated(this%eff_fr)) deallocate(this%eff_fr)
    allocate(this%eff_fr(4,this%nfr))
    
    this%eff_fr(1,:) = getRealArrParam(keys,values,"eff_fr_cib_clustering",this%nfr)
    this%eff_fr(2,:) = getRealArrParam(keys,values,"eff_fr_cib_poisson",this%nfr)
    this%eff_fr(3,:) = getRealArrParam(keys,values,"eff_fr_radio_poisson",this%nfr)
    this%eff_fr(4,:) = getRealArrParam(keys,values,"eff_fr_tsz",this%nfr)
    this%norm_fr(1) = getRealParam(keys,values,"norm_fr_cib_clustering")
    this%norm_fr(2) = getRealParam(keys,values,"norm_fr_cib_poisson")
    this%norm_fr(3) = getRealParam(keys,values,"norm_fr_radio_poisson")
    this%norm_fr(4) = getRealParam(keys,values,"norm_fr_tsz")

    ! Radio galaxy flux cuts
    this%rg_flux_cut = getRealParam(keys,values,"rg_flux_cut")
    this%norm_rg_flux_cut = getRealParam(keys,values,"norm_rg_flux_cut")
	
		allocate(this%cib_decor_clust(this%nfr,this%nfr))
		allocate(this%cib_decor_poisson(this%nfr,this%nfr))
		this%cib_decor_clust = 1
		this%cib_decor_poisson = 1
		
  end subroutine 	initialize_internal

	subroutine set_cib_decor_clust(instance_id,mat_clust,nfr)
		integer,intent(in)::instance_id,nfr
		real(8),dimension(nfr,nfr),intent(in)::mat_clust
		integer::i,j
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
		
		
		do i=1,this%nfr
			do j=1,this%nfr
				this%cib_decor_clust(i,j) = mat_clust(i,j)
			enddo
		enddo	
		
	end subroutine set_cib_decor_clust
	
  subroutine set_cib_decor_poisson(instance_id,mat_poisson,nfr)
		integer,intent(in)::instance_id,nfr
		real(8),dimension(nfr,nfr),intent(in)::mat_poisson
		integer::i,j
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
		
		
		do i=1,this%nfr
			do j=1,this%nfr
				this%cib_decor_poisson(i,j) = mat_poisson(i,j)
			enddo
		enddo	
		
	end subroutine set_cib_decor_poisson
  
  subroutine initialize_frommem(instance_id,keys,values,nkv,cib_clustering,patchy_ksz,homogenous_ksz,tsz,error)
    integer,intent(inout) :: instance_id, error
    integer,intent(in) :: nkv
    character(*), dimension(nkv),intent(in) :: keys
    character(*), dimension(nkv),intent(in) :: values
		real(8),dimension(1:10000),intent(in)::cib_clustering,patchy_ksz
		real(8),dimension(1:10000*6),intent(in)::homogenous_ksz
		real(8),dimension(1:10000*7),intent(in)::tsz
    integer :: l,i
    real :: dum
    type(foreground_params), pointer :: this
    
    this => foreground_instances(instance_id)

    error=0

		do l=2,10000
			this%clust_dg_templ(l) = cib_clustering(l)
			this%ksz_patchy_templ(l) = patchy_ksz(l)
			do i=1,6
				this%ksz_ov_templ(l,i) = homogenous_ksz((l-1)*6+i)
			enddo
			do i=1,7
				this%tsz_templ(l,i) = tsz((l-1)*7+i)
			enddo
		enddo 
		
    ! Load effective band centers and normalization frequencies
		
		call initialize_internal(instance_id,keys,values,nkv,error)
		
  end subroutine initialize_frommem

  subroutine initialize(instance_id,keys,values,nkv,error)
    integer,intent(inout) :: instance_id, error
    integer,intent(in) :: nkv
    character(*), dimension(nkv),intent(in) :: keys
    character(*), dimension(nkv),intent(in) :: values
    integer :: l,i
    real :: dum
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)

    error=0

    ! Clustered template
    Open(10,file=getStrParam(keys,values,"template_cib_clustering"),err=20)
    Do l = 2,10000
      Read (10,*) dum, this%clust_dg_templ(l)
    End Do
    Close(10)

    ! Patchy
    Open(10,file=getStrParam(keys,values,"template_patchy_ksz"),err=20)
    Do l = 2,10000
      Read (10,*) dum, this%ksz_patchy_templ(l)
    End Do
    Close(10)

    ! OV
    Open(10,file=getStrParam(keys,values,"template_homogenous_ksz"),err=20)
    Do l = 2,10000
      Read (10,*) dum, this%ksz_ov_templ(l,:)
    End Do
    Close(10)


    ! TSZ template
    Open(10,file=getStrParam(keys,values,"template_tsz"),err=20)
    Do l = 2,10000
      Read (10,*) dum, this%tsz_templ(l,:)
    End Do
    Close(10)

    call initialize_internal(instance_id,keys,values,nkv,error)
    return

20  error=1
  end subroutine

  !
  ! Gets the number of parameters in the model,
  ! such that dR returned in other functions has last dimension np
  !
  ! The which integer specifies for which component:
  !  1 - cib-clustered
  !  2 - cib-poisson
  !  3 - radio poisson
  !  4 - tsz
  !  5 - ksz
  !  6 - total
  !
  subroutine get_np(instance_id,which,np)
    integer :: instance_id
    integer :: np
    integer :: which

    np = NUM_FG_PARAMS(which)
  end subroutine

  !
  ! Gets the number of frequencies
  !
  subroutine get_nfr(instance_id,nfr)
    integer :: instance_id, nfr
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)

    nfr = this%nfr

  end subroutine



  !
  ! Get the CIB-clustering contribution
  !
  subroutine get_cib_clustering(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR
    real(8), dimension(nfr,nfr) :: frqdep, logfac, normfac
    real(8), dimension(nfr) :: fr
    real(8) :: fr0
    real(8) :: alpha,tilt,norm
    integer :: i,j,l
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
    call update_foreground_params(this,keys,values)

    if (.not. checkSizes(this,nfr,1,np,error)) return

    
    fr = this%eff_fr(1,:)
    fr0 = this%norm_fr(1)
    alpha = this%alpha_dg_cl%v
    tilt = this%tilt_dg_cl%v
    norm = this%norm_dg_cl%v


    do i=1,this%nfr
      do j=1,this%nfr
        frqdep(i,j) = (fr(i)*fr(j)/fr0**2)**alpha/dBdT(fr(i),fr0)/dBdT(fr(j),fr0) * this%cib_decor_clust(i,j)
        logfac(i,j) = log(fr(i)*fr(j)/fr0**2)
      end do
    end do


    do l=lmin,lmax
      ! Make SURE things are properly normalized at l=3000, whether
       ! template is previously normalized or not...
      normfac = (max(1500,l)/3000.)**tilt * this%clust_dg_templ(l)&
           &/(d3000*this%clust_dg_templ(3000)) * frqdep
      R(:,:,l) = norm * normfac 

      dR(:,:,l,this%alpha_dg_cl%i) = norm * normfac * logfac
      dR(:,:,l,this%norm_dg_cl%i) = normfac
      if (l<1500) then
        dR(:,:,l,this%tilt_dg_cl%i) = norm * normfac * log(1500./3000.)
      else
        dR(:,:,l,this%tilt_dg_cl%i) = norm * normfac * log(l/3000.)
      end if
    end do
  end subroutine

  !
  ! Get the CIB-poisson contribution
  !
  subroutine get_cib_poisson(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR
    real(8), dimension(nfr,nfr) :: frqdep, logfac, normfac
    real(8), dimension(nfr) :: fr
    real(8) :: fr0
    real(8) :: alpha, sigma, norm
    integer :: i,j,l
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
    call update_foreground_params(this,keys,values)

    if (.not. checkSizes(this,nfr,2,np,error)) return

    fr = this%eff_fr(2,:)
    fr0 = this%norm_fr(2)
    alpha = this%alpha_dg_po%v
    sigma = this%sigma_dg_po%v
    norm = this%norm_dg_po%v

    do i=1,this%nfr
      do j=1,this%nfr
        frqdep(i,j) = (fr(i)*fr(j)/fr0**2)**(alpha + sigma**2*log(fr(i)*fr(j)/fr0**2)/2)/dBdT(fr(i),fr0)/dBdT(fr(j),fr0) * this%cib_decor_poisson(i,j)
        logfac(i,j) = log(fr(i)*fr(j)/fr0**2)
      end do
    end do

    do l=lmin,lmax
      normfac = 1./d3000 * frqdep
      R(:,:,l) = norm * normfac
      dR(:,:,l,this%alpha_dg_po%i) = norm * normfac * logfac
      dR(:,:,l,this%sigma_dg_po%i) = norm * normfac * sigma * logfac**2
      dR(:,:,l,this%norm_dg_po%i) = normfac
    end do
  end subroutine


  !
  ! Get the radio-poisson contribution
  !
  subroutine get_radio_poisson(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np,ff
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR
    real(8), dimension(nfr,nfr) :: frqdep, logfac, normfac
    real(8), dimension(nfr) :: fr
    real(8) :: fr0
    real(8) :: alpha, sigma, gam, norm
    integer :: i,j,l
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
    call update_foreground_params(this,keys,values)

    if (.not. checkSizes(this,nfr,3,np,error)) return

    fr = this%eff_fr(3,:)
    fr0 = this%norm_fr(3)
    alpha = this%alpha_rg%v
    sigma = this%sigma_rg%v
    norm = this%norm_rg%v
    gam = this%gamma_rg%v

    do i=1,this%nfr
      do j=1,this%nfr
        frqdep(i,j) = (fr(i)*fr(j)/fr0**2) ** (alpha + sigma**2 * log(fr(i)*fr(j)/fr0**2)/2 )/dBdT(fr(i),fr0)/dBdT(fr(j),fr0)
        !write(*,*) i,j,frqdep(i,j),fr(i),fr(j),dBdT(fr(i),fr0),dBdT(fr(j),fr0)
        logfac(i,j) = log(fr(i)*fr(j)/fr0**2)
      end do
    end do

    !do ff=1,nfr
    !  write(*,*) "in1",1001,1,ff,R(1,ff,1001)
    !  write(*,*) frqdep(1,ff),frqdep(ff,1),fr(ff)
    !enddo
    do l=lmin,lmax
      normfac = 1./d3000 * (this%rg_flux_cut/this%norm_rg_flux_cut)**(gam+2) * frqdep
      R(:,:,l) = norm * normfac
      dR(:,:,l,this%norm_rg%i) = normfac
      dR(:,:,l,this%alpha_rg%i) = norm * normfac * logfac
      dR(:,:,l,this%sigma_rg%i) = norm * normfac * sigma * logfac**2
      dR(:,:,l,this%gamma_rg%i) = norm * normfac * log(this%rg_flux_cut/this%norm_rg_flux_cut) * (this%rg_flux_cut/this%norm_rg_flux_cut)**2
    end do
    !do ff=1,nfr
    !  write(*,*) "in",1001,1,ff,R(1,ff,1001)
    !enddo
  end subroutine


  !
  ! Get the thermal SZ contribution
  !
  subroutine get_tsz(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR
    real(8), dimension(nfr,nfr) :: frqdep
    real(8), dimension(nfr) :: fr
    real(8) :: fr0
    integer :: i,j,l
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
    call update_foreground_params(this,keys,values)

    if (.not. checkSizes(this,nfr,4,np,error)) return

    fr = this%eff_fr(3,:)
    fr0 = this%norm_fr(3)

    do i=1,this%nfr
      do j=1,this%nfr
        frqdep(i,j) = tszFreqDep(fr(i),fr0)*tszFreqDep(fr(j),fr0)
      end do
    end do

    do l=lmin,lmax
      R(:,:,l) = (this%tsz_mean_scale%v * this%tsz_templ(l,1) + this%tsz_pca1%v * this%tsz_templ(l,2) + this%tsz_pca2%v * this%tsz_templ(l,3))
      dR(:,:,l,this%tsz_mean_scale%i) = this%tsz_templ(l,1)
      dR(:,:,l,this%tsz_pca1%i) = this%tsz_templ(l,2)
      dR(:,:,l,this%tsz_pca2%i) = this%tsz_templ(l,3)
    end do

  end subroutine


  !
  ! Get the kinetic SZ contribution
  !
  subroutine get_ksz(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR
    real(8), dimension(nfr,nfr) :: frqdep
    integer :: i,j,l,lshift
    type(foreground_params), pointer :: this
    this => foreground_instances(instance_id)
    call update_foreground_params(this,keys,values)

    if (.not. checkSizes(this,nfr,5,np,error)) return
    frqdep = 1
    lshift = max(2,min(10000,int(this%shift_patchy%v * l)))
    do l=lmin,lmax
      R(:,:,l) = this%norm_ov%v * this%ksz_ov_templ(l,1) + this%norm_patchy%v * this%ksz_patchy_templ(lshift)
      dR(:,:,l,this%norm_ov%i) = this%ksz_ov_templ(l,1)
      dR(:,:,l,this%norm_patchy%i) = this%ksz_patchy_templ(lshift)
      dR(:,:,l,this%shift_patchy%i) = 0
    end do
  end subroutine

  !
  ! Get the total extra-galactic foregroud contribution
  !
  subroutine get_egfs(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR

    real(8), dimension(:,:,:), allocatable :: Rdc, Rdp, Rrp, Rtsz, Rksz
    real(8), dimension(:,:,:,:), allocatable :: dRdc, dRdp, dRrp, dRtsz, dRksz
    allocate(Rdc(nfr,nfr,lmin:lmax),&
              Rdp(nfr,nfr,lmin:lmax),&
              Rrp(nfr,nfr,lmin:lmax),&
              Rtsz(nfr,nfr,lmin:lmax),&
              Rksz(nfr,nfr,lmin:lmax))
    allocate(dRdc(nfr,nfr,lmin:lmax,NUM_FG_PARAMS(1)),&
              dRdp(nfr,nfr,lmin:lmax,NUM_FG_PARAMS(2)),&
              dRrp(nfr,nfr,lmin:lmax,NUM_FG_PARAMS(3)),&
              dRtsz(nfr,nfr,lmin:lmax,NUM_FG_PARAMS(4)),&
              dRksz(nfr,nfr,lmin:lmax,NUM_FG_PARAMS(5)))

    call get_cib_clustering(instance_id,keys,values,nkv,Rdc,dRdc,nfr,lmin,lmax,NUM_FG_PARAMS(1),error)
    call get_cib_poisson(instance_id,keys,values,nkv,Rdp,dRdp,nfr,lmin,lmax,NUM_FG_PARAMS(2),error)
    call get_radio_poisson(instance_id,keys,values,nkv,Rrp,dRrp,nfr,lmin,lmax,NUM_FG_PARAMS(3),error)
    call get_tsz(instance_id,keys,values,nkv,Rtsz,dRtsz,nfr,lmin,lmax,NUM_FG_PARAMS(4),error)
    call get_ksz(instance_id,keys,values,nkv,Rksz,dRksz,nfr,lmin,lmax,NUM_FG_PARAMS(5),error)

    R = Rdc + Rdp + Rrp + Rtsz + Rksz

    dR(:,:,:,1:NUM_FG_PARAMS(1)) = dRdc
    dR(:,:,:,NUM_FG_PARAMS(1)+1:NUM_FG_PARAMS(1)+NUM_FG_PARAMS(2)) = dRdp
    dR(:,:,:,sum(NUM_FG_PARAMS(:2))+1:sum(NUM_FG_PARAMS(:2))+NUM_FG_PARAMS(3)) = dRrp
    dR(:,:,:,sum(NUM_FG_PARAMS(:3))+1:sum(NUM_FG_PARAMS(:3))+NUM_FG_PARAMS(4)) = dRtsz
    dR(:,:,:,sum(NUM_FG_PARAMS(:4))+1:sum(NUM_FG_PARAMS(:4))+NUM_FG_PARAMS(5)) = dRksz

    deallocate(Rdc, Rdp, Rrp, Rtsz, Rksz)
    deallocate(dRdc, dRdp, dRrp, dRtsz, dRksz)
  end subroutine


  subroutine get_egfs_component(instance_id,comp,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    integer :: instance_id, error
    integer :: nkv, nfr, lmin, lmax, np, comp
    character(*), dimension(nkv) :: keys
    character(*), dimension(nkv) :: values
    real(8), dimension(nfr,nfr,lmin:lmax) :: R
    real(8), dimension(nfr,nfr,lmin:lmax,np) :: dR

    if (comp==1) then
      call get_cib_clustering(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    elseif (comp==2) then
      call get_cib_poisson(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    elseif (comp==3) then
      call get_radio_poisson(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    elseif (comp==4) then
      call get_tsz(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    elseif (comp==5) then
      call get_ksz(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    elseif (comp==6) then
      call get_egfs(instance_id,keys,values,nkv,R,dR,nfr,lmin,lmax,np,error)
    end if
  end subroutine



  subroutine update_foreground_params(this,keys,values)
    integer :: instance_id
    character(*), dimension(:) :: keys
    character(*), dimension(:) :: values
    type(foreground_params), pointer :: this

    !StringJoin["if (hasParam(keys,values,\"" <> # <> "\")) this%" <> # <> "%v = getRealParam(keys,values,\"" <> # <> "\")\n" & /@ ps]
    !StringJoin["if (hasParam(keys,values,\"idx_"<>#<>"\")) this%"<>#<>"%i = getIntParam(keys,values,\"idx_"<>#<>"\")\n"&/@ps]

    if (hasParam(keys,values,"alpha_dg_cl")) this%alpha_dg_cl%v = getRealParam(keys,values,"alpha_dg_cl")
    if (hasParam(keys,values,"alpha_dg_po")) this%alpha_dg_po%v = getRealParam(keys,values,"alpha_dg_po")
    if (hasParam(keys,values,"alpha_rg")) this%alpha_rg%v = getRealParam(keys,values,"alpha_rg")
    if (hasParam(keys,values,"fpol_dg_po")) this%fpol_dg_po%v = getRealParam(keys,values,"fpol_dg_po")
    if (hasParam(keys,values,"fpol_rg")) this%fpol_rg%v = getRealParam(keys,values,"fpol_rg")
    if (hasParam(keys,values,"gamma_rg")) this%gamma_rg%v = getRealParam(keys,values,"gamma_rg")
    if (hasParam(keys,values,"norm_dg_cl")) this%norm_dg_cl%v = getRealParam(keys,values,"norm_dg_cl")
    if (hasParam(keys,values,"norm_dg_po")) this%norm_dg_po%v = getRealParam(keys,values,"norm_dg_po")
    if (hasParam(keys,values,"norm_ov")) this%norm_ov%v = getRealParam(keys,values,"norm_ov")
    if (hasParam(keys,values,"norm_patchy")) this%norm_patchy%v = getRealParam(keys,values,"norm_patchy")
    if (hasParam(keys,values,"norm_rg")) this%norm_rg%v = getRealParam(keys,values,"norm_rg")
    if (hasParam(keys,values,"shift_patchy")) this%shift_patchy%v = getRealParam(keys,values,"shift_patchy")
    if (hasParam(keys,values,"sigma_dg_po")) this%sigma_dg_po%v = getRealParam(keys,values,"sigma_dg_po")
    if (hasParam(keys,values,"sigma_rg")) this%sigma_rg%v = getRealParam(keys,values,"sigma_rg")
    if (hasParam(keys,values,"tilt_dg_cl")) this%tilt_dg_cl%v = getRealParam(keys,values,"tilt_dg_cl")
    if (hasParam(keys,values,"tsz_dgcl_cor")) this%tsz_dgcl_cor%v = getRealParam(keys,values,"tsz_dgcl_cor")
    if (hasParam(keys,values,"tsz_mean_scale")) this%tsz_mean_scale%v = getRealParam(keys,values,"tsz_mean_scale")
    if (hasParam(keys,values,"tsz_pca1")) this%tsz_pca1%v = getRealParam(keys,values,"tsz_pca1")
    if (hasParam(keys,values,"tsz_pca2")) this%tsz_pca2%v = getRealParam(keys,values,"tsz_pca2")

    if (hasParam(keys,values,"idx_alpha_dg_cl")) this%alpha_dg_cl%i = getIntParam(keys,values,"idx_alpha_dg_cl")
    if (hasParam(keys,values,"idx_tilt_dg_cl")) this%tilt_dg_cl%i = getIntParam(keys,values,"idx_tilt_dg_cl")
    if (hasParam(keys,values,"idx_alpha_dg_po")) this%alpha_dg_po%i = getIntParam(keys,values,"idx_alpha_dg_po")
    if (hasParam(keys,values,"idx_sigma_dg_po")) this%sigma_dg_po%i = getIntParam(keys,values,"idx_sigma_dg_po")
    if (hasParam(keys,values,"idx_norm_dg_cl")) this%norm_dg_cl%i = getIntParam(keys,values,"idx_norm_dg_cl")
    if (hasParam(keys,values,"idx_norm_dg_po")) this%norm_dg_po%i = getIntParam(keys,values,"idx_norm_dg_po")
    if (hasParam(keys,values,"idx_fpol_dg_po")) this%fpol_dg_po%i = getIntParam(keys,values,"idx_fpol_dg_po")
    if (hasParam(keys,values,"idx_alpha_rg")) this%alpha_rg%i = getIntParam(keys,values,"idx_alpha_rg")
    if (hasParam(keys,values,"idx_sigma_rg")) this%sigma_rg%i = getIntParam(keys,values,"idx_sigma_rg")
    if (hasParam(keys,values,"idx_norm_rg")) this%norm_rg%i = getIntParam(keys,values,"idx_norm_rg")
    if (hasParam(keys,values,"idx_gamma_rg")) this%gamma_rg%i = getIntParam(keys,values,"idx_gamma_rg")
    if (hasParam(keys,values,"idx_fpol_rg")) this%fpol_rg%i = getIntParam(keys,values,"idx_fpol_rg")
    if (hasParam(keys,values,"idx_tsz_pca1")) this%tsz_pca1%i = getIntParam(keys,values,"idx_tsz_pca1")
    if (hasParam(keys,values,"idx_tsz_pca2")) this%tsz_pca2%i = getIntParam(keys,values,"idx_tsz_pca2")
    if (hasParam(keys,values,"idx_norm_ov")) this%norm_ov%i = getIntParam(keys,values,"idx_norm_ov")
    if (hasParam(keys,values,"idx_norm_patchy")) this%norm_patchy%i = getIntParam(keys,values,"idx_norm_patchy")
    if (hasParam(keys,values,"idx_shift_patchy")) this%shift_patchy%i = getIntParam(keys,values,"idx_shift_patchy")
    if (hasParam(keys,values,"idx_tsz_dgcl_cor")) this%tsz_dgcl_cor%i = getIntParam(keys,values,"idx_tsz_dgcl_cor")
    if (hasParam(keys,values,"idx_tsz_mean_scale")) this%tsz_mean_scale%i = getIntParam(keys,values,"idx_tsz_mean_scale")

  end subroutine


  !
  ! nu,nu0 in GHz
  !
  ! dBdT is proportional to derivative of planck function
  ! but is normalized so its equal to 1 at nu0
  !
  function dBdT(nu,nu0)
    real(8) x, x0, dBdT, dBdT0, nu, nu0
    x0 = nu0/56.78
    dBdT0 = x0**4 * exp(x0) / (exp(x0)-1)**2
    x = nu/56.78
    dBdT = x**4 * exp(x) / (exp(x)-1)**2 / dbdT0
  end function

  !
  ! nu, nu0 in GHz
  ! Gives the tsz frequency dependence normalized so its 1 at nu0
  !
  function tszFreqDep(nu,nu0)
    real(8) :: tszFreqDep, tszFreqDep0, nu, nu0, x, x0
    x = nu / 56.78
    x0 = nu0 / 56.78
    tszFreqDep0 = x0*(exp(x0)+1)/(exp(x0)-1) - 4
    tszFreqDep = x*(exp(x)+1)/(exp(x)-1) - 4
    tszFreqDep = tszFreqDep/tszFreqDep0
  end function


  !
  ! Check that R and dR matrices are the right size
  ! Sets error = 1 and returns false if they are not
  !
  function checkSizes(this,nfr,which,np,error)
    integer nfr, np, which, error
    type(foreground_params), pointer :: this
    logical :: checkSizes

    checkSizes=.true.

    if (this%nfr /= nfr) then
      print *, "R matrix size",nfr," not the same as nfr parameter", this%nfr
      error=1
      checkSizes=.false.
    end if

    if (np /= NUM_FG_PARAMS(which)) then
      print *, "dR matrix size",np," for component",which," is not the same as the number of parameters", NUM_FG_PARAMS(which)
      error=1
      checkSizes=.false.
    end if
  end function

  ! Normalizes a template so D_3000 = 1
  subroutine NormalizeTemplate(templ)
      real(8), dimension(:) :: templ
      templ = templ/(d3000*templ(3000))
  end subroutine

end module


