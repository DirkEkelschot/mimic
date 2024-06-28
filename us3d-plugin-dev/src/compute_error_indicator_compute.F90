module computeErrorIndicator_compute
    use us3d_global
    use us3d_varkinds
    use us3d_component_module
    use us3d_signals_module
    use us3d_handler_class
    use us3d_register_slots_module
    use computeErrorIndicator_global
    implicit none
    integer                              :: nsvars 
    integer, dimension(:), allocatable   :: ismap
    real(rk),allocatable,dimension(:,:,:),target :: HessianT
    contains

        
    subroutine error_indicator_postupdate(ier)
        use mpivars, only : id
        use switches, only : nit 
        use sizing, only   : net,nel,neg,nbf,nff,nif,nmap_grad_invf,bu_size
        use simvars, only  : bigu,grad,rhop
        use geometry, only : xce,xcf
        use connect, only : ief,ife,iee
        
        integer i,ii,j,jj,epg,n
        !logical, intent(OUT) :: new
        integer(KIND=US3D_GINT) m
        real(rk) wei(1)
        real(rk) rhs(3), xcc(3,1), dT2dxi2(3), dT2dyi2(3), dT2dzi2(3)
        real(rk) dTdxi(3), dTdxip1(3),dTdxi_diff(3)
        !real(rk),allocatable,dimension(:,:,:),target :: HessianT
        !real(rk) :: bigu1(bu_size),bigu2(bu_size)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        !implicit none
        integer :: ier 
        ier = 0
        if(id.eq.0)print'("----compute error indicator----")' 
        !new= .false.
        epg= nel + neg
        if(id.eq.0)write(*,*)"busize",bu_size
        allocate(HessianT(3,3,net))
        HessianT = 0._rk
        do m = 1, nff
            i  = ife(m,1)
            ii = ife(m,2)

            !! This is the vector from cell i to cell ii
            xcc(:,1) = xce(:,ii) - xce(:,i) ! local coordinates

            !! This is the weight (1/distance squared)
            wei(1) = 1._rk/(dot_product(xcc(:,1),xcc(:,1)))

            !! Set-up the BigU variables for cell i and ii
            dTdxi       = grad(:,5,i)
            dTdxip1     = grad(:,5,ii)
            dTdxi_diff  = dTdxip1(:) - dTdxi(:)
            !if(id.eq.0)write(*,*)"busize",dTdxi_diff
            ! For the ife bucketing, grad(:,:,i) serves the same purpose
            ! as big_rhs.
            !< When adding to the i side of the face,  xcc points i to ii, so this
            !  the direction that the vector is pointing. When adding to the ii
            !  side of the face, both xcc and bigud should be degated (opposite
            !  direction).  Since they are only used when multiplied together,
            !  this cancels and it can just be straight added.
            !>
            do n=1,3
                rhs(:)= xcc(:,1)*dTdxi_diff(n)*wei(1)
                HessianT(:,n,i)  = HessianT(:,n,i)  + rhs(:)
                HessianT(:,n,ii) = HessianT(:,n,ii) + rhs(:)
            enddo
         
        enddo
        do i = 1, nel
            do n= 1,3
               rhs(:)= HessianT(:,n,i)
               HessianT(1,n,i)= dot_product(rhop(1,:,i),rhs(:))
               HessianT(2,n,i)= dot_product(rhop(2,:,i),rhs(:))
               HessianT(3,n,i)= dot_product(rhop(3,:,i),rhs(:))
            enddo
        enddo

        !do i = 1, nel 
        !    do n= 1,3 
        !       write(*,*)"Hessian",HessianT(1,n,i),HessianT(2,n,i),HessianT(3,n,i)
        !    enddo
        !enddo

       ! if(mod(nit,error_indicator_write_frequency)==0 .and. do_ErrorIndicator) call write_error_data()
        return
    end subroutine error_indicator_postupdate


end module computeErrorIndicator_compute
