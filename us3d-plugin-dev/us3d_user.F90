function us3d_plugin_init(component,idebug,iactive,comm) result(ier) &
   bind(c,name="us3d_plugin_init")
   use, intrinsic :: iso_c_binding
   use computeErrorIndicator_setup, only : error_indicator_initialize
   use mpi
   use mpivars
   implicit none
   character(c_char), intent(IN) :: component(*)
   integer(c_int), intent(IN)    :: idebug
   integer(c_int), intent(INOUT) :: iactive
   integer(c_int), intent(IN)    :: comm
   integer(c_int)                :: ier
   integer                       :: i
   logical                       :: debug
   character(LEN=40)             :: mycomp
      ier = 0   
      mycomp(:)= ' '
      print'("Are we registered?")'
      do i = 1,len(mycomp)
         if (component(i)==c_null_char) exit
         mycomp(i:i)= component(i)
      enddo
      select case(mycomp)
      case('us3d-solver')
         debug= .false.
         if (idebug==1) debug= .true.
         call error_indicator_initialize(mycomp,debug,ier)
         if (ier==0)iactive= 1
      end select
   return
end function us3d_plugin_init
!%-------------------------------------------------------------------%
function us3d_plugin_dest(component) result(ier) &
   bind(c,name="us3d_plugin_dest")
   use, intrinsic :: iso_c_binding
   implicit none
   character(c_char), intent(IN) :: component(*)
   integer(c_int)                :: ier
   ier = 0
   return
end function us3d_plugin_dest
!%-------------------------------------------------------------------%
function us3d_plugin_query(component,qstring) result(ptr) &
   bind(c,name="us3d_plugin_query")
   use, intrinsic :: iso_c_binding
   implicit none
   character(c_char), intent(IN) :: component(*)
   character(c_char), intent(IN) :: qstring(*)
   type(c_ptr)                   :: ptr
   return
end function us3d_plugin_query
