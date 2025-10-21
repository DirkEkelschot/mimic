        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 28 14:34:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE US3D_PLUGIN_DEST__genmod
          INTERFACE 
            FUNCTION US3D_PLUGIN_DEST(COMPONENT) RESULT(IER)            &
     & BIND(C, NAME = 'us3d_plugin_dest')
              CHARACTER(LEN=1), INTENT(IN) :: COMPONENT(*)
              INTEGER(KIND=4) :: IER
            END FUNCTION US3D_PLUGIN_DEST
          END INTERFACE 
        END MODULE US3D_PLUGIN_DEST__genmod
