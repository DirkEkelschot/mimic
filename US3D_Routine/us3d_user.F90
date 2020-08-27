#include "us3d_header.h"
module userdata
      implicit none
      save
      contains
      !  ***********************************************************************************
      ! The following subroutine is intended to be called after everything is done by the
      ! solver within a timestep. Once the variables are updated, we can repalce the residual 
      ! variable with our Jacobian
      Subroutine my_postupdate(istage)
            use sizing,   only : nel ! number of elements
            use geometry, only : xcn,dw,voli ! node centroids (3,nnds)
            use connect,  only : ien ! element to node mapping
            Use mpivars,  only : id  ! processor ID
            use simvars,  only : res ! L2 norm of density residual that we will replace with Jacobian
            use, intrinsic :: ISO_C_Binding
            Implicit none
            ! We create an interface for the desired C++ function that we wish to call
            interface
                  function ComputeJAtCenter (P, np) result(J) bind(C, name="ComputeJAtCenter")
                  use, intrinsic :: ISO_C_Binding
                  implicit none
                  integer(c_int),value :: np
                  real(C_double), dimension(3*8) :: P
                  real(C_double) :: J
                  end function ComputeJAtCenter
            end interface
            integer :: i,j 
            integer, intent(IN) :: istage
            integer(C_INT) :: numpoints
            double precision :: Points(3*8)
            real(C_double) :: Jacobian
            type(C_ptr) :: Pptr,npPtr
            double precision :: maxJac,OldJac

            res = 0.0d0
            OldJac = 0.0d0
            do i = 1, nel
                  if(dw(i).gt.1.0D-3)then
                        numpoints = ien(0,i)
                        do j = 1, ien(0,i)
                              Points((j-1)*3+1) = xcn(1,ien(j,i))
                              Points((j-1)*3+2) = xcn(2,ien(j,i))
                              Points((j-1)*3+3) = xcn(3,ien(j,i))
                        enddo
                        Jacobian = 0.0d0
                        Jacobian = ComputeJAtCenter(Points,numpoints)
                        res(i) = Jacobian*voli(i)
                        maxJac = max(Jacobian,OldJac)
                        OldJac = Jacobian
                  endif
            enddo
            print*,maxJac 
            return
      end Subroutine my_postupdate
end module userdata
!  ***********************************************************************************
subroutine user_initialize(component,debug,ier)
      use us3d_user_module, only: user_update_post
      use userdata, only: my_postupdate
      implicit none
      character(LEN=*), intent(IN) :: component
      logical, intent(IN) :: debug
      integer, intent(OUT) :: ier
      ier= 0
      user_update_post => my_postupdate
      return
end subroutine user_initialize