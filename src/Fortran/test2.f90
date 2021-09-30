! Program to test the function getLogAddExp_RK

#include "Constants_mod.f90"
#include "log_mod.f90"

program test2

    use iso_fortran_env, only: output_unit
    use Constants_mod, only: IK, RK
    use log_mod, only: getLogAddExp_RK

    Implicit None

    real(RK)                    :: directResult
    real(RK)                    :: functionResult
    real(RK)                    :: term1
    real(RK)                    :: term2
    
    term1 = log(2._RK) + 3._RK
    term2 = log(4._RK) + 5._RK
    
    functionResult = getLogAddExp_RK(term2, term1)
    directResult = log(2._RK*exp(3._RK) + 4._RK*exp(5._RK))
    
    write(output_unit,"(*(g0))") "Function result: ", functionResult
    write(output_unit,"(*(g0))") "Direct result: ", directResult

end