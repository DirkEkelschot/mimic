        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 28 14:34:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE US3D_PLUGIN_INIT__genmod
          INTERFACE 
            FUNCTION US3D_PLUGIN_INIT(COMPONENT,IDEBUG,IACTIVE,COMM)    &
     & RESULT(IER) BIND(C, NAME = 'us3d_plugin_init')
              CHARACTER(LEN=1), INTENT(IN) :: COMPONENT(*)
              INTEGER(KIND=4), INTENT(IN) :: IDEBUG
              INTEGER(KIND=4), INTENT(INOUT) :: IACTIVE
              INTEGER(KIND=4), INTENT(IN) :: COMM
              INTEGER(KIND=4) :: IER
            END FUNCTION US3D_PLUGIN_INIT
          END INTERFACE 
        END MODULE US3D_PLUGIN_INIT__genmod
