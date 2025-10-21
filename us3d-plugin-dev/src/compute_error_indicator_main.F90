module computeErrorIndicator_main_m
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

     subroutine error_indicator_readin(ier)
            use mpivars, only : icomw,id
            use switches
            use US3D_EXTRAS, only : us3d_scanto_sect, us3d_check_int
            use TXT_FUNCS
            use IO_FUNCS
            implicit none
            integer             :: ier
            integer             :: lun,iline,ieof,ios
            integer             :: i,j,k,m,tmp(100)
            logical             :: iexist
            character(len=200)  :: text, key, val, cval

            do while(id.eq.0)
                  call io_find_unused_lun(lun)
                  open(unit=lun,file=inpf,status='old')
                  iline= 0

                  write(6,*) '-- Scanning for COMPUTE ERROR ESTIMATION section'


            enddo
            return
    end subroutine error_indicator_readin

    subroutine error_indicator_main(ier)
        use simvars,  only : t
        integer :: ier

        
        return
    end subroutine error_indicator_main


    function compute_hessian() result(HessianT)
        use computeErrorIndicator_io, only : write_error_data
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
        real(rk),allocatable,dimension(:,:,:) :: HessianT
        !real(rk) :: bigu1(bu_size),bigu2(bu_size)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        !implicit none
        integer :: ier 
        ier = 0
        !if(id.eq.0)print'("----compute error indicator----")' 
        !new= .false.
        epg= nel + neg
        !if(id.eq.0)write(*,*)"busize huhuh ",bu_size
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
    end function compute_hessian   



    function test_compute_hessian() result(Temp)
        use computeErrorIndicator_io, only : write_error_data
        use mpivars, only : id
        use switches, only : nit 
        use sizing, only   : net,nel,neg,nbf,nff,nif,nmap_grad_invf,bu_size
        use simvars, only  : bigu,grad,rhop,t
        use geometry, only : xce,xcf
        use connect
        use models, only : ns,bcs,jbcmap
        
        integer i,ii,j,jj,epg,n
        !logical, intent(OUT) :: new
        integer(KIND=US3D_GINT) m
        real(rk) wei(1)
        real(rk) rhs(3), xcc(3,1), dT2dxi2(3), dT2dyi2(3), dT2dzi2(3)
        real(rk) dTdxi(3), dTdxip1(3),dTdxi_diff(3)
        real(rk) dT(1), dTp1(1),dT_diff(1)
        real(rk),allocatable,dimension(:,:,:) :: Temp
        !real(rk) :: bigu1(bu_size),bigu2(bu_size)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        !implicit none
        integer :: ier 
        ier = 0
        !if(id.eq.0)print'("----compute error indicator----")' 
        !new= .false.
        epg= nel + neg
        !if(id.eq.0)write(*,*)"busize",bu_size,ugrid%glob_nfmax
        allocate(Temp(3,1,net))
        Temp = 0._rk
        !write(*,*)"bu_size ",bu_size,ugrid%glob_nfmax
        do m = 1, nff
            i  = ife(m,1)
            ii = ife(m,2)

            !! This is the vector from cell i to cell ii
            xcc(:,1) = xce(:,ii) - xce(:,i) ! local coordinates

            !! This is the weight (1/distance squared)
            wei(1) = 1._rk/(dot_product(xcc(:,1),xcc(:,1)))

            !! Set-up the BigU variables for cell i and ii
            dT       = t(i)
            dTp1     = t(ii)
            
            dT_diff  = dTp1(:) - dT(:)
            !if(id.eq.0)write(*,*)"busize",dTdxi_diff
            ! For the ife bucketing, grad(:,:,i) serves the same purpose
            ! as big_rhs.
            !< When adding to the i side of the face,  xcc points i to ii, so this
            !  the direction that the vector is pointing. When adding to the ii
            !  side of the face, both xcc and bigud should be degated (opposite
            !  direction).  Since they are only used when multiplied together,
            !  this cancels and it can just be straight added.
            !>

            !do n=1,3
            rhs(:)= xcc(:,1)*dT_diff(1)*wei(1)
            Temp(:,1,i)  = Temp(:,1,i)  + rhs(:)
            Temp(:,1,ii) = Temp(:,1,ii) + rhs(:)
            !enddo
         
        enddo
        do i = 1, nel
            !do n= 1,3
            rhs(:)= Temp(:,1,i)
            Temp(1,1,i)= dot_product(rhop(1,:,i),rhs(:))
            Temp(2,1,i)= dot_product(rhop(2,:,i),rhs(:))
            Temp(3,1,i)= dot_product(rhop(3,:,i),rhs(:))
            
            !enddo
        enddo

       !do i = 1, nel 
       !     do n= 1,3 
       !        write(*,*)"temps ", Temp(1,1,i), Temp(2,1,i), Temp(3,1,i), grad(1,5,i), grad(2,5,i), grad(3,5,i)
       !     enddo
       ! enddo

       ! if(mod(nit,error_indicator_write_frequency)==0 .and. do_ErrorIndicator) call write_error_data()

        
        
        
        
        
        return
    end function test_compute_hessian  




    function test_compute_hessian_v2() result(Temp)
        use computeErrorIndicator_io, only : write_error_data
        use us3d_global
        use mpivars
        use sizing
        use connect
        use US3D_GRID
        use geometry
        use switches
        use us3d_core_module
        !use mpivars, only : id
        !use switches, only : nit 
        !use sizing, only   : net,nel,neg,nbf,nff,nif,nmap_grad_invf,bu_size
        use simvars, only  : bigu,grad,rhop,t
        !use geometry, only : xce,xcf
        !use connect, only : ief,ife,iee,ugrid,ilst
        
        integer i,ii,j,jj,epg,n
        !logical, intent(OUT) :: new
        integer(KIND=US3D_GINT) m,ib,ie,nfaces,k
        real(rk) wei(1)
        real(rk) rhs(3), xcc(3,1), dT2dxi2(3), dT2dyi2(3), dT2dzi2(3)
        real(rk) dTdxi(3), dTdxip1(3),dTdxi_diff(3)
        real(rk) dT(1), dTp1(1),dT_diff(1),big_rhs(1,3)
        real(rk),allocatable,dimension(:,:,:) :: Temp
        !real(rk) :: bigu1(1),bigu2(1)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        integer(KIND=US3D_GINT) list(0:4)
        !integer(KIND=US3D_GINT) list(0:6)
        !implicit none
        integer :: ier 
        ier = 0
        !if(id.eq.0)print'("----compute error indicator----")' 
        !new= .false.
        epg= nel + neg
        !if(id.eq.0)write(*,*)"busize",bu_size,ugrid%glob_nfmax
        allocate(Temp(3,1,net))
        Temp = 0._rk
        
        !write(*,*)"ugrid%glob_nfmax", ugrid%glob_nfmax

        do m = 1,2
            if(m.eq.1) then                        ! bounds for interior cells
            ib = 1
            ie = nel
            else                                   ! bounds for shared cells
            ib = nel+neg+1
            ie = net
            endif

            Do i = ib,ie                       ! loop over selected cells
            list = ilst(:,i)                   ! pick up the list

            if(list(0) .ge. 4) then            ! not a singular operator
                do k = 1,list(0)                ! loop over cells in list
                   ii = list(k)

                   xcc(:,k) = xce(:,ii) - xce(:,i) ! local coordinates
                   wei(k) = 1._rk/(dot_product(xcc(:,k),xcc(:,k))) 
                enddo

                big_rhs= 0._rk
                dT       = t(i)
                dTp1     = t(ii)
                

                rhs = 0._rk                     ! initialize right hand side
                do k = 1,list(0)                ! loop over cells in list
                    ii = list(k)                 ! get cell number

                    ! bigud= bigu2(:) - bigu1(:)

                    dT_diff  = dTp1(:) - dT(:)

                    big_rhs(:,1)= big_rhs(:,1) + xcc(1,k)*dT_diff(:)*wei(k)
                    big_rhs(:,2)= big_rhs(:,2) + xcc(2,k)*dT_diff(:)*wei(k)
                    big_rhs(:,3)= big_rhs(:,3) + xcc(3,k)*dT_diff(:)*wei(k)
                enddo


                do n= 1,1
                    rhs(:)= big_rhs(imap_grad_invf(n),:)
                    Temp(1,n,i)= dot_product(rhop(1,:,i),rhs(:))
                    Temp(2,n,i)= dot_product(rhop(2,:,i),rhs(:))
                    Temp(3,n,i)= dot_product(rhop(3,:,i),rhs(:))
                enddo                
            else
                Temp(:,:,i) = 0._rk
            endif

            Enddo
        enddo
        
        
        return
    end function test_compute_hessian_v2


    subroutine error_indicator_postupdate(ier)
        use computeErrorIndicator_io, only : write_error_data
        use mpivars, only : id
        use switches, only : nit 
        use sizing, only   : net,nel,neg,nbf,nff,nif,nmap_grad_invf,bu_size
        use simvars, only  : bigu,grad,rhop,t
        use geometry, only : xce,xcf
        use connect, only : ief,ife,iee,ugrid
        
        integer i,ii,j,jj,epg,n
        !logical, intent(OUT) :: new
        integer(KIND=US3D_GINT) m
        real(rk) wei(1)
        real(rk) rhs(3), xcc(3,1)
        real(rk) dTdxi(1), dTdxip1(1),dTdxi_diff(1)
        real(rk),allocatable,dimension(:,:,:),target :: Temp
        !real(rk) :: bigu1(bu_size),bigu2(bu_size)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        !implicit none
        integer :: ier 
        ier = 0
        if(id.eq.0)print'("postupdate----compute error indicator----")'
        !new= .false.
        epg= nel + neg
        !if(id.eq.0)write(*,*)"busize",bu_size
        allocate(Temp(3,1,net))
        Temp = 0._rk
        do m = 1, nff
            i  = ife(m,1)
            ii = ife(m,2)
            !write(*,*)ugrid%glob_nfmax
            !! This is the vector from cell i to cell ii
            xcc(:,1) = xce(:,ii) - xce(:,i) ! local coordinates

            !! This is the weight (1/distance squared)
            wei(1) = 1._rk/(dot_product(xcc(:,1),xcc(:,1)))

            !! Set-up the BigU variables for cell i and ii
            dTdxi       = t(i)
            dTdxip1     = t(ii)
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
            !do n=1,3
            rhs(:)= xcc(:,1)*dTdxi_diff(1)*wei(1)
            Temp(:,1,i)  = Temp(:,1,i)  + rhs(:)
            Temp(:,1,ii) = Temp(:,1,ii) + rhs(:)
            !enddo
         
        enddo
        do i = 1, nel
            !do n= 1,3
               rhs(:)= Temp(:,1,i)
               Temp(1,1,i)= dot_product(rhop(1,:,i),rhs(:))
               Temp(2,1,i)= dot_product(rhop(2,:,i),rhs(:))
               Temp(3,1,i)= dot_product(rhop(3,:,i),rhs(:))
            !enddo
        enddo

        !do i = 1, nel 
        !    do n= 1,3 
        !       write(*,*)"Hessian",HessianT(1,n,i),HessianT(2,n,i),HessianT(3,n,i)
        !    enddo
        !enddo

       ! if(mod(nit,error_indicator_write_frequency)==0 .and. do_ErrorIndicator) call write_error_data()
        return
    end subroutine error_indicator_postupdate



    subroutine error_indicator_postupdate_v2(ier)
        use computeErrorIndicator_io, only : write_error_data
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
        real(rk),allocatable,dimension(:,:,:),target :: HessianT
        !real(rk) :: bigu1(bu_size),bigu2(bu_size)
        real(8), dimension(:,:), allocatable :: my_svars,my_bvars
        !implicit none
        integer :: ier 
        ier = 0
        if(id.eq.0)print'("postupdate_v2----compute error indicator----")' 
        !new= .false.
        epg= nel + neg
        if(id.eq.0)write(*,*)"busizei error_indicator_postupdate_v2 ",bu_size
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
    end subroutine error_indicator_postupdate_v2


    !  ***********************************************************************************
    !  Listener routine for adding particle-related variables to US3D as output.
    !  (only for version 1.1.7 and not used for RC22)
    !  Currently, hardcoded for bx, by, bz, ex, ey, ez, psi, phi
    !  ***********************************************************************************
    subroutine dTdxi_listener(signal,package)

      Use mpivars
      Use switches, only : us3d_debug

    
      implicit none

      integer                           :: ier, ib
      character(LEN=*),      intent(IN) :: signal
      integer, dimension(:), intent(IN) :: package
      character(20), allocatable        :: VarNames(:)
      if (id==0.and.us3d_debug) then
         write(6,*) ' MaxListener was triggered:'
         write(6,*) '   signal: '//trim(signal)
      endif

      ! Set variabels required for storing extra Maxwell variables
      nsvars = 9

      ALLOCATE(VarNames(nsvars))
      ALLOCATE(   ismap(nsvars))

      VarNames( 1) ='dT2dxx'
      VarNames( 2) ='dT2dxy'
      VarNames( 3) ='dT2dxz'
      VarNames( 4) ='dT2dyx'
      VarNames( 5) ='dT2dyy'
      VarNames( 6) ='dT2dyz'
      VarNames( 7) ='dT2dzx'
      VarNames( 8) ='dT2dzy'
      VarNames( 9) ='dT2dzz' 
      select case(signal)

      case('us3d_save_cb_request')

         !! At this signal, requesting additional storage for saving solution 
         !! variables (in slots_save_cells). Then in AddPartVars handler routine, 
         !! these variables are assigned
         !! Interior + ghost cell data (should be handled through the IO object in future)
         call us3d_register_slots('save_cells', 'dTdxi', VarNames, ier, locs = ismap)

         if (ier/=0) stop
         if (id==0.and.us3d_debug) write(6,*) '   - Got ismap= ', ismap
          
      end select

    end subroutine dTdxi_listener




    subroutine my_addvars_sub(this)

      use sizing, only   : net
      Use mpivars
      Use switches, only : us3d_debug
      use simvars, only: grad
!      use computeErrorIndicator_compute only: HessianT
      !real(rk),allocatable,dimension(:,:,:) :: HessianT
      !real(rk) HessianT(3,3,net)
      implicit none
      class(us3d_handler) :: this
      integer :: i,j,k,n,ii,jj
      real(rk) HessianT(3,3,net)
      real(rk) Temp(3,1,net)
      !! Note that on reading, we should really compare the string this%text to
      !! see the position of old svars or bvars, in case something has changed
      !! and the variables are in different locations in the restart file than
      !! where we expect them for writing.  This could happen if a new plugin is
      !! activated on restart which also adds variables which weren't present
      !! when we wrote the solution.  To be safe, we should check.  However, as
      !! a stand-alone example, this works for now.

      this%ier= 0
      !HessianT  = compute_hessian()
      Temp      = test_compute_hessian()
      !Temp     = test_compute_hessian_v2()
      select case(this%iloc)
      
      ! -- Save user solution variables in interior
!      case (US3D_ILOC_DATAIO_SVARS)
!     index = row*3+col

!         do ii= 1,this%i1
!            i = this%ivars(ii)
!              do j= 1,3
!                do k=1,3
!                   n = (j-1)*3+k
!                   this%rv2(ismap(n),ii)=HessianT(j,k,i)
!                enddo
!              enddo
!         enddo
!      case default
         ! Currently this plugin will display a warning here because the data IO
         ! routines will call the handler to add variables to slices too, but
         ! we have not yet implemented that in this plugin code.
         !if (id==0) write(6,*) 'WWW unhandled '//us3dh_iloc_to_text_f(this%iloc)
!      end select


      case (US3D_ILOC_DATAIO_SVARS)

      do ii= 1,this%i1
         i = this%ivars(ii)
           do j= 1,3
             !write(*,*)"ismap(j)",ismap(j),ismap(j)+3,ismap(j)+6
             this%rv2(ismap(j),ii)=Temp(j,1,i)
             this%rv2(ismap(j)+3,ii)=grad(j,5,i)
             this%rv2(ismap(j)+6,ii)=Temp(j,1,i)-grad(j,5,i)
             !write(*,*) Temp(j,1,i), grad(j,5,i)
          enddo
     enddo
    case default
         ! Currently this plugin will display a warning here because the data IO
         ! routines will call the handler to add variables to slices too, but
         ! we have not yet implemented that in this plugin code.
         !if (id==0) write(6,*) 'WWW unhandled '//us3dh_iloc_to_text_f(this%iloc)
    end select


    end subroutine my_addvars_sub






    subroutine error_indicator_dataio_post(irnum,fname,cio,ier)
        integer             :: irnum,ier
        character(1)        :: cio
        character(len=200)  :: fname

    end subroutine error_indicator_dataio_post


end module computeErrorIndicator_main_m
