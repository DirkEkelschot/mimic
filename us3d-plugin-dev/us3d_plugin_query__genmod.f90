        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 28 14:34:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE US3D_PLUGIN_QUERY__genmod
          INTERFACE 
            FUNCTION US3D_PLUGIN_QUERY(COMPONENT,QSTRING) RESULT(PTR)   &
     & BIND(C, NAME = 'us3d_plugin_query')
              USE ISO_C_BINDING
              CHARACTER(LEN=1), INTENT(IN) :: COMPONENT(*)
              CHARACTER(LEN=1), INTENT(IN) :: QSTRING(*)
              TYPE (C_PTR) :: PTR
            END FUNCTION US3D_PLUGIN_QUERY
          END INTERFACE 
        END MODULE US3D_PLUGIN_QUERY__genmod
