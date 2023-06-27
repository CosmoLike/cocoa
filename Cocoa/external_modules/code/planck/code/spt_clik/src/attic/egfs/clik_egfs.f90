module clik_egfs_extra
	use egfs
	integer::cefgs=1
	character(*),dimension(3),parameter:: cib_clustering_paramname = (/ 'alpha_dg_cl', 'tilt_dg_cl ','norm_dg_cl '/)
	character(*),dimension(4),parameter:: cib_poisson_paramname = (/ 'alpha_dg_po', 'sigma_dg_po','norm_dg_po ','fpol_dg_po '/)
	character(*),dimension(5),parameter:: radio_poisson_paramname = (/ 'alpha_rg', 'sigma_rg','norm_rg ','gamma_rg','fpol_rg '/)
	character(*),dimension(3),parameter:: tsz_paramname = (/ 'tsz_pca1      ', 'tsz_pca2      ','tsz_mean_scale'/)
	character(*),dimension(3),parameter:: ksz_paramname = (/ 'norm_ov     ', 'norm_patchy ','shift_patchy'/)
contains
	subroutine fill_kv_from_c(c_keys,c_values,nkv,keys,values)
		character(len=*),intent(in)::c_keys
		character(len=*),intent(in)::c_values
	  integer,intent(inout)::nkv
		character(len=265),dimension(:),pointer,intent(out)::keys
		character(len=265),dimension(:),pointer,intent(out)::values
		character(len=256)::dum
		integer::ikv,okv
		integer::tkv

		tkv = nkv + size(cib_clustering_paramname) +size(cib_poisson_paramname) +size(radio_poisson_paramname) +size(tsz_paramname) +size(ksz_paramname)

		allocate(keys(tkv))
		allocate(values(tkv))

		do ikv=0,nkv-1
			dum = c_keys(ikv*256+1:ikv*256+256)
			dum(len_trim(dum):len_trim(dum))=' '
			keys(ikv+1)=TRIM(ADJUSTL(dum))
			dum = c_values(ikv*256+1:ikv*256+256)
			dum(len_trim(dum):len_trim(dum))=' '
			values(ikv+1)=TRIM(ADJUSTL(dum))
		enddo
		okv = nkv
		do ikv=1,size(cib_clustering_paramname)
			keys(okv+ikv) = 'idx_'//trim(cib_clustering_paramname(ikv))
			write(dum,'(I2)') ikv 
			values(okv+ikv) = TRIM(ADJUSTL(dum))
		enddo
		okv = okv + size(cib_clustering_paramname)
		do ikv=1,size(cib_poisson_paramname)
			keys(okv+ikv) = 'idx_'//trim(cib_poisson_paramname(ikv))
			write(dum,'(I2)') ikv 
			values(okv+ikv) = TRIM(ADJUSTL(dum))
		enddo
		okv = okv + size(cib_poisson_paramname)
		do ikv=1,size(radio_poisson_paramname)
			keys(okv+ikv) = 'idx_'//trim(radio_poisson_paramname(ikv))
			write(dum,'(I2)') ikv 
			values(okv+ikv) = TRIM(ADJUSTL(dum))
		enddo
		okv = okv + size(radio_poisson_paramname)
		do ikv=1,size(tsz_paramname)
			keys(okv+ikv) = 'idx_'//trim(tsz_paramname(ikv))
			write(dum,'(I2)') ikv 
			values(okv+ikv) = TRIM(ADJUSTL(dum))
		enddo
		okv = okv + size(tsz_paramname)
		do ikv=1,size(ksz_paramname)
			keys(okv+ikv) = 'idx_'//trim(ksz_paramname(ikv))
			write(dum,'(I2)') ikv 
			values(okv+ikv) = TRIM(ADJUSTL(dum))
		enddo
		!do ikv=1,tkv
		!	WRITE(*,*) trim(keys(ikv)),'=',trim(values(ikv))
		!enddo
		nkv = tkv
	end subroutine fill_kv_from_c
	subroutine clik_egfs_inhere(key,li,id)
		integer,intent(out)::id
		character(len=*),intent(in)::key
		character(len=*),dimension(:),intent(in)::li
		integer::i

		do i=1,size(li)
			if (trim(key)==trim(li(i))) THEN
				id = id+i
				return
			endif
		enddo	
	end subroutine
end module clik_egfs_extra

subroutine clik_egfs_init(instance_id,c_keys,c_values,nkv)
	use clik_egfs_extra
	integer,intent(out) :: instance_id
	integer,intent(in) :: nkv
	character(len=*),intent(in)::c_keys
	character(len=*),intent(in)::c_values
	integer::errid
	integer::tkv
	character(len=265),dimension(:),pointer::keys
	character(len=265),dimension(:),pointer::values
	
	call get_next_instance(instance_id)
	if (instance_id==0) THEN
		instance_id=0
		return
	endif
	tkv=nkv
	call fill_kv_from_c(c_keys,c_values,tkv,keys,values)
	call initialize(instance_id,keys,values,tkv,errid)
	deallocate(keys)
	deallocate(values)
	if (errid/=0) THEN
		instance_id=0
	endif

end subroutine

subroutine clik_egfs_init_frommem(instance_id,c_keys,c_values,nkv,cib_clustering,patchy_ksz,homogenous_ksz,tsz)
	use clik_egfs_extra
	integer,intent(out) :: instance_id
	integer,intent(in) :: nkv
	character(len=*),intent(in)::c_keys
	character(len=*),intent(in)::c_values
	real(8),dimension(1:10000),intent(in)::cib_clustering,patchy_ksz
	real(8),dimension(1:10000*6),intent(in)::homogenous_ksz
	real(8),dimension(1:10000*7),intent(in)::tsz

	
	integer::errid
	integer::tkv
	character(len=265),dimension(:),pointer::keys
	character(len=265),dimension(:),pointer::values
	
	call get_next_instance(instance_id)
	if (instance_id==0) THEN
		instance_id=0
		return
	endif
	tkv=nkv
	call fill_kv_from_c(c_keys,c_values,tkv,keys,values)
	call initialize_frommem(instance_id,keys,values,tkv,cib_clustering,patchy_ksz,homogenous_ksz,tsz,errid)
	deallocate(keys)
	deallocate(values)
	if (errid/=0) THEN
		instance_id=0
	endif

end subroutine


subroutine clik_get_egfs_component(instance_id,cmp,c_keys,c_values,nkv,R,dR,nfr,lmin,lmax,np,errid)
	use clik_egfs_extra
  integer,intent(in) :: instance_id,nkv, nfr, lmin, lmax, np
  character(*),intent(in) :: c_keys
  character(*),intent(in) :: c_values
  real(8), dimension(nfr,nfr,lmin:lmax),intent(in) :: R
  real(8), dimension(nfr,nfr,lmin:lmax,np),intent(in) :: dR
	integer,intent(out)::errid
	integer,intent(in)::cmp
	character(len=265),dimension(:),pointer::keys
	character(len=265),dimension(:),pointer::values
	integer::tkv,ll,ff
	
	tkv=nkv
	call fill_kv_from_c(c_keys,c_values,tkv,keys,values)
	call get_egfs_component(instance_id,cmp,keys,values,tkv,R,dR,nfr,lmin,lmax,np,errid)
	!do ff=1,nfr
	!	write(*,*) 1001,1,ff,R(1,ff,1001)
	!enddo
	!do ff=1,nfr
	!	do ll=lmin,lmax
	!		write(*,*) ff,ll,R(ff,ff,ll),dR(ff,ff,ll,1)
	!	enddo
	!enddo
	deallocate(keys)
	deallocate(values)
end subroutine

subroutine clik_get_nder(nm, nd)
	use clik_egfs_extra
	integer,intent(in)::nm
	integer,intent(out)::nd
	integer::i
	nd = NUM_FG_PARAMS(nm)
end subroutine



subroutine clik_egfs_order(nac,c_keys,cid)
	use clik_egfs_extra
	integer,intent(in)::nac
	integer,intent(out),dimension(nac)::cid
	character(*),intent(in)::c_keys
	character(256)::key
	integer::i,j,id
	
	do i=1,nac
		key = c_keys((i-1)*256+1:i*256)
		key(len_trim(key):len_trim(key))=' '
		id = 0
		call clik_egfs_inhere(key,cib_clustering_paramname,id)
		if (id/=0) then
			cid(i) = id
			cycle
		endif
		id = 0
		call clik_egfs_inhere(key,cib_poisson_paramname,id)
		if (id/=0) then
			cid(i) = id +10
			cycle
		endif
		id = 0
		call clik_egfs_inhere(key,radio_poisson_paramname,id)
		if (id/=0) then
			cid(i) = id +20
			cycle
		endif
		id = 0
		call clik_egfs_inhere(key,tsz_paramname,id)
		if (id/=0) then
			cid(i) = id +30
			cycle
		endif
		id = 0
		call clik_egfs_inhere(key,ksz_paramname,id)
		if (id/=0) then
			cid(i) = id +40
			cycle
		endif
		id = 0
		cid(i)=id
	enddo
end subroutine clik_egfs_order

subroutine clik_get_nfr(instance_id,nfr)
	use clik_egfs_extra
	integer,intent(in)::instance_id
	integer,intent(out)::nfr
	call get_nfr(instance_id,nfr)
end subroutine

subroutine clik_free_instance(instance_id)
	use clik_egfs_extra
	integer,intent(in)::instance_id
	call free_instance(instance_id)
end subroutine clik_free_instance

subroutine clik_set_cib_decor_clust(instance_id,mat_clust,nfr)
	use clik_egfs_extra
	integer,intent(in)::instance_id,nfr
	real(8),dimension(nfr,nfr),intent(in)::mat_clust
	call set_cib_decor_clust(instance_id,mat_clust,nfr)
end subroutine clik_set_cib_decor_clust

subroutine clik_set_cib_decor_poisson(instance_id,mat_poisson,nfr)
	use clik_egfs_extra
	integer,intent(in)::instance_id,nfr
	real(8),dimension(nfr,nfr),intent(in)::mat_poisson
	call set_cib_decor_poisson(instance_id,mat_poisson,nfr)
end subroutine clik_set_cib_decor_poisson
!