#include "us3d_header.h"
module userdata
      implicit none
      double precision, parameter :: pi=4.0d0*atan(1.0d0)
      save
      contains
      !  ***********************************************************************************
      ! The following subroutine is intended to be called after everything is done by the
      ! solver within a timestep. Once the variables are updated, we can repalce the residual 
      ! variable with our Jacobian
      Subroutine my_postupdate(istage)
            use sizing,   only : nel ! number of elements
            use geometry, only : xcn ! node centroids (3,nnds)
            use connect,  only : ien ! element to node mapping
            Use mpivars,  only : id  ! processor ID
            use simvars, only : res  ! L2 norm of density residual that we will replace with Jacobian
            use, intrinsic :: ISO_C_Binding
            Implicit none
            ! We create an interface for the desired C++ function that we wish to call
            interface
                  function ComputeJ (P, ElType) result(J) bind(C,name="ComputeJ")
                  use, intrinsic :: ISO_C_Binding
                  import
                  integer(C_int) :: ElType
                  real(C_double) :: P(0:3*8-1)
                  real(C_double) :: J
                  end function ComputeJ
            end interface
            integer :: i,j 
            integer, intent(IN) :: istage
            integer(C_int)  :: ElType
            real(C_double) :: Points(0:3*8-1)
            real(C_double) :: Jacobian
            double precision :: maxJac,OldJac

            OldJac = 0.0d0
            do i = 1, nel
                  do j = 0, ien(0,i) - 1
                        Points(j*3+0) = xcn(1,ien(j+1,i))
                        Points(j*3+1) = xcn(2,ien(j+1,i))
                        Points(j*3+2) = xcn(3,ien(j+1,i))
                  enddo
                  Jacobian = 0.0d0
                  Jacobian = ComputeJ(Points,ElType)
                  res(i) = Jacobian
                  maxJac = max(Jacobian,OldJac)
                  OldJac = Jacobian
            enddo
            if(id.eq.0)print*,maxJac
            return
      end Subroutine my_postupdate
end module userdata
!  ***********************************************************************************
!     User initialization routine
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
