#include "us3d_header.h"
module computeErrorIndicator_io
    use us3d_global
    use us3d_varkinds
    use us3d_component_module
    use us3d_signals_module
    use us3d_handler_class
    use us3d_register_slots_module
    implicit none
#include "interface_us3d_core.f"
#include "interface_us3d_handler_sub.f"
    save
    contains
    subroutine write_error_data()
        use computeErrorIndicator_global
        use mpivars, only : id
        use switches, only : nit 
        use sizing, only   : net,nel,neg,nbf,nif
        use simvars, only  : bigu
        use geometry, only : xce,xcf
        integer i,j,jj,k,epg,ier
        integer, parameter :: nuvars=2  ! Number of BigU variables
        integer, parameter :: nsvars=2  ! Number of interior (save/load) variables
        integer, parameter :: nbvars=4  ! Number of boundary (save/load) variables
        integer :: iumap(nuvars),ismap(nsvars),ibmap(nbvars)
        real(rk), dimension(:,:), allocatable :: my_svars,my_bvars
!        class(us3d_handler) :: this
        ier = 0 
!        if (.NOT. do_ErrorIndicator) return
!        ! Write error indicator data to file
!        if(mod(nit,error_indicator_write_frequency)==0 .and. do_ErrorIndicator) call write_error_data()
    
!        integer i,j,jj,k,epg,ier
    
         epg= nel + neg 

         !! Here we're just storing some dummy values for testing

         do i = 1,net
            do k= 1,nuvars
               bigu(iumap(k),i)= dble(k)
            enddo
         enddo

         do i= 1,epg
            do k= 1,nsvars
               my_svars(k,i)= xce(k,i)
               write(*,*) "xce! ",xce(k,i)
            enddo
         enddo

         do jj= 1,nbf
            j= jj + nif 
            do k= 1,nsvars
               my_bvars(k,jj)= xcf(k,j)
            enddo
            do k= nsvars+1,nbvars
               my_bvars(k,jj)= dble(2*k)
            enddo
         enddo
    
    
        return
        end subroutine write_error_data

end module computeErrorIndicator_io
