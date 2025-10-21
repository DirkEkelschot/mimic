module computeErrorIndicator_setup

   use us3d_handler_class, only : us3d_handler   
   use, intrinsic :: iso_c_binding
   implicit none
   save
   private
   character(LEN=20),public :: plugin_name   = 'Error Indicator Routines'
   character(LEN=43),public :: plugin_author = 'AMA Dirk Ekelschot, NASA Ames'
   ! Extend the us3d_handler with  type
   type, extends(us3d_handler), public :: ff_handler_type
   contains
   private
   end type ff_handler_type
   ! Define initialization handler
   type(ff_handler_type), pointer :: ff_handle_init 
   procedure(us3d_handler_sub), pointer :: ei_readin    => NULL()
   procedure(us3d_handler_sub), pointer :: ei_startup   => NULL()
   procedure(us3d_handler_sub), pointer :: ei_main_pre   => NULL()
   procedure(us3d_handler_sub), pointer :: ei_main_post   => NULL()
   procedure(us3d_handler_sub), pointer :: ei_main_postupdate => NULL()
   procedure(us3d_handler_sub), pointer :: ei_dataio    => NULL()
   procedure(us3d_handler_sub), pointer :: ei_util      => NULL()
   procedure(us3d_handler_sub),pointer :: my_hsub
#include "interface_us3d_core.f"

   procedure(us3d_core_signal),pointer :: fptr => NULL()
   public :: error_indicator_initialize
   !%-------------------------------------------------------------------%
   abstract interface
   subroutine us3d_handler_sub(h)
   import :: us3d_handler
   class(us3d_handler) :: h
   end subroutine us3d_handler_sub

   end interface
   contains
!  ***********************************************************************************
!     User initialization routine
!  ***********************************************************************************
   subroutine error_indicator_initialize(component,debug,ier)
        use us3d_user_module
        use us3d_register_procs_module
        use us3d_signals_module
        use us3d_register_handlers_module
        !use us3d_register_handlers_module, only : us3d_add_handler
        use mpivars, only : id
        use computeErrorIndicator_driver_m, only : error_indicator_us3d_postupdate, &
                                                 error_indicator_us3d_dataio_post, &
                                                dTdxi_listener, &
                                                my_addvars_sub
        implicit none
        character(LEN=*), intent(IN) :: component
        logical, intent(IN) :: debug
        integer, intent(OUT) :: ier 
        !procedure(us3d_core_signal), pointer :: fptr
        


        ei_readin => ei_readin
        call us3d_add_handler('us3dh_readbc_post',istat=ier,priority=1,hsub=ei_readin)

        ei_main_post => error_indicator_us3d_postupdate
        call us3d_add_handler('us3dh_main_postupdate',istat=ier,priority=10,hsub=ei_main_post)

        ei_dataio    => error_indicator_us3d_dataio_post
        call us3d_add_handler('us3dh_dataio_post',istat=ier,priority=10,hsub=ei_dataio)

        fptr => dTdxi_listener
        ! Listen for signals to add particle varaibles
        call us3d_signal_listen('us3d_save_cb_request', 'dTdxi-data', fptr, cname = 'us3d-solver')

        
        my_hsub => my_addvars_sub

        call us3d_add_handler('us3dh_dataio_addvars', &
         hsub=my_hsub,desc='Addvars example',author='Heath Johnson',istat=ier)
        if (ier/=0) stop


   end subroutine error_indicator_initialize

end module computeErrorIndicator_setup


