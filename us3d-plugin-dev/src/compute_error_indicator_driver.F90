module computeErrorIndicator_driver_m
    use us3d_handler_class, only : us3d_handler
    use computeErrorIndicator_main_m
    use computeErrorIndicator_global
    implicit none
    contains


    subroutine error_indicator_us3d_readin(h)
        implicit none
        integer :: ier 
        class(us3d_handler) :: h
        ier = 0 
        call error_indicator_readin(ier)
        h%ier = ier 
        return
    end subroutine error_indicator_us3d_readin




    subroutine error_indicator_us3d_postupdate(h)
        implicit none
        integer :: ier 
        class(us3d_handler) :: h
        ier = 0 
        call error_indicator_postupdate(ier)
        h%ier = ier 
        return
    end subroutine error_indicator_us3d_postupdate



    subroutine error_indicator_us3d_dataio_post(h)
        implicit none
        class(us3d_handler) :: h
        integer             :: irnum,ier
        character(len=200)  :: fname
        character(1)        :: cio
        ier = 0
        if(.not.do_ErrorIndicator) return
        irnum = h%irnum
        fname = h%text
        Select Case(h%iop)
        Case(1)
            cio = 'w'
        Case(2)
            cio = 'r'
        End Select
        call error_indicator_dataio_post(irnum,fname,cio,ier)
        h%ier = ier
        return
    end subroutine error_indicator_us3d_dataio_post

end module computeErrorIndicator_driver_m
