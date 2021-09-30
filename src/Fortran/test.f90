! Program to test the functions getCholeskyFactor and getDeterminant

#include "Constants_mod.f90"
#include "Matrix_mod.f90"

program test

    use iso_fortran_env, only: output_unit
    use Constants_mod, only: IK, RK
    use Matrix_mod, only: getCholeskyFactor, getDeterminant

    Implicit None

    integer(IK) , parameter     :: NDIM = 3_IK

    real(RK)                    :: matrix(NDIM,NDIM)
    real(RK)                    :: CholeskyDiag(NDIM)
    real(RK)                    :: sqrtDetMatrix
    real(RK)                    :: dummy
    
    matrix(1,1) = 4_RK
    matrix(1,2) = 12_RK
    matrix(1,3) = -16_RK
    matrix(2,1) = 12_RK
    matrix(2,2) = 37_RK
    matrix(2,3) = -43_RK
    matrix(3,1) = -16_RK
    matrix(3,2) = -43_RK
    matrix(3,3) = 98_RK
    
    call getCholeskyFactor(nd = NDIM, PosDefMat = matrix, Diagonal = CholeskyDiag)
    sqrtDetMatrix = product(CholeskyDiag)
    
    write(output_unit,"(*(g0))") "sqrt of det of matrix via getCholeskyFactor is: ", sqrtDetMatrix
    
    sqrtDetMatrix = sqrt(getDeterminant(nd = NDIM, InputMat = matrix))
    
    write(output_unit,"(*(g0))") "sqrt of det of matrix via getDeterminant is: ", sqrtDetMatrix
    
end
