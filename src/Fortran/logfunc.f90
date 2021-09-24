!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%  Description:
!%      +   Return the natural logarithm of the NPAR-dimensional logLikelihood function.
!%  Input:
!%      +   npar:       The number of dimensions of the domain of the objective function.
!%      +   point:      The input 64-bit real-valued vector of length ndim,
!%                      at which the natural logarithm of objective function is computed.
!%  Output:
!%      +   logFunc:    A 64-bit real scalar number representing the natural logarithm of the objective function.
!%  Author:
!%      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
!%  Visit:
!%      +   https://www.cdslab.org/paramonte
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Constants_mod.f90"
#include "Matrix_mod.f90"

module LogFunc_mod

    use iso_fortran_env, only: output_unit
    use Constants_mod, only: IK, RK, SQRT2PI

    implicit none

    integer(IK) , parameter     :: NDIM = 2_IK ! number of observational attributes: Dur, Epk
    integer(IK) , parameter     :: NPAR = 12_IK ! number of parameters of the model
    integer(IK) , parameter     :: NDATA = 1966_IK ! Batse_orig data size
    real(RK)    , parameter     :: COEF0 = NDIM * log( 1._RK / SQRT2PI )

    real(RK)                    :: mv_meanVec(NDIM)
    real(RK)                    :: mv_stdDur
    real(RK)                    :: mv_stdEpk
    real(RK)                    :: mv_rho
    real(RK)                    :: mv_amp
    real(RK)                    :: mv_meanVec2(NDIM)
    real(RK)                    :: mv_stdDur2
    real(RK)                    :: mv_stdEpk2
    real(RK)                    :: mv_rho2
    real(RK)                    :: mv_amp2
    real(RK)                    :: mv_CovMatUpperCholeskyLower(NDIM,NDIM)
    real(RK)                    :: mv_CholeskyDiag(NDIM)
    real(RK)                    :: mv_CovMatUpperCholeskyLower2(NDIM,NDIM)
    real(RK)                    :: mv_CholeskyDiag2(NDIM)
    real(RK)                    :: mv_InvCovMat(NDIM,NDIM)
    real(RK)                    :: mv_InvCovMat2(NDIM,NDIM)
    
    ! data attributes

    type :: GRB_type
        real(RK) :: T50(NDATA)
        real(RK) :: T90(NDATA)
        real(RK) :: Epk(NDATA)
        real(RK) :: LogT50(NDATA)
        real(RK) :: LogT90(NDATA)
        real(RK) :: LogEpk(NDATA)
    end type GRB_type 

    type(GRB_type) :: GRB

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine readData() ! Read BATSE data

        implicit none

        character(:), allocatable   :: header
        integer(IK)                 :: inputFileUnit
        integer(IK)                 :: idata

        ! Read L19 data from input data

        open( newunit = inputFileUnit &
            , file = "../../in/Batse_orig.csv" &
            , status = "old" &
            )

        allocate(character(100) :: header)

        read(inputFileUnit,*) header; header = trim(adjustl(header))
        !write(*,*) header
        do idata = 1, NDATA
            read(inputFileUnit,*) GRB%T50(idata), GRB%T90(idata), GRB%Epk(idata)
            GRB%LogT50(idata) = log(GRB%T50(idata))
            GRB%LogT90(idata) = log(GRB%T90(idata))
            GRB%LogEpk(idata) = log(GRB%Epk(idata))
        end do

        close(inputFileUnit)

    end subroutine readData

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogFunc(npar,param) result(logFunc)
        ! Return the negative natural logarithm of MVN distribution evaluated at the input vector point.

        use Constants_mod, only: NEGINF_RK
        use Matrix_mod, only: getCholeskyFactor, getInvMatFromCholFac

        implicit none

        integer(IK) , intent(in)    :: npar
        real(RK)    , intent(in)    :: param(npar)
        real(RK)                    :: logFunc
        integer(IK)                 :: idata, i, j
        real(RK)                    :: sqrtDetCov
        real(RK)                    :: sqrtDetCov2
        real(RK)                    :: NormedPoint(NDIM)
        real(RK)                    :: NormedPoint2(NDIM)

        ! Order of variables in param:
        !
        !   param(1) = logAvg(x)
        !   param(2) = logAvg(y)
        !   param(3) = logStd(x)
        !   param(4) = logStd(y)
        !   param(5) = FirsherTrans(self._rho)
        !   param(6) = amplitude
        !   param(7) = logAvg2(x)
        !   param(8) = logAvg2(y)
        !   param(9) = logStd2(x)
        !   param(10) = logStd2(y)
        !   param(11) = FirsherTrans(self._rho2)
        !   param(12) = amplitude2

        mv_meanVec(1) = param(1)
        mv_meanVec(2) = param(2)
        mv_stdDur = exp(param(3))
        mv_stdEpk = exp(param(4))
        mv_rho = tanh(param(5))
        mv_amp = param(6)
        mv_meanVec2(1) = param(7)
        mv_meanVec2(2) = param(8)
        mv_stdDur2 = exp(param(9))
        mv_stdEpk2 = exp(param(10))
        mv_rho2 = tanh(param(11))
        mv_amp2 = param(12)
        
        mv_CovMatUpperCholeskyLower(1,1) = mv_stdDur**2
        mv_CovMatUpperCholeskyLower(2,2) = mv_stdEpk**2
        mv_CovMatUpperCholeskyLower(1,2) = mv_rho * mv_stdDur * mv_stdEpk
        mv_CovMatUpperCholeskyLower(2,1) = mv_CovMatUpperCholeskyLower(1,2)
        mv_CovMatUpperCholeskyLower2(1,1) = mv_stdDur2**2
        mv_CovMatUpperCholeskyLower2(2,2) = mv_stdEpk2**2
        mv_CovMatUpperCholeskyLower2(1,2) = mv_rho2 * mv_stdDur2 * mv_stdEpk2
        mv_CovMatUpperCholeskyLower2(2,1) = mv_CovMatUpperCholeskyLower2(1,2)
        
        ! compute the Cholescky factor of mv_CovMatUpperCholeskyLower
        call getCholeskyFactor(nd = NDIM, PosDefMat = mv_CovMatUpperCholeskyLower, Diagonal = mv_CholeskyDiag)
        if (mv_CholeskyDiag(1)<0._RK) then
            logFunc = NEGINF_RK
            return
        end if
        call getCholeskyFactor(nd = NDIM, PosDefMat = mv_CovMatUpperCholeskyLower2, Diagonal = mv_CholeskyDiag2)
        if (mv_CholeskyDiag2(1)<0._RK) then
            logFunc = NEGINF_RK
            return
        end if
        
        ! The following is the log of the coefficient used in the definition of the MVN.

        ! coef = self.coef0 + np.log( np.sqrt(np.linalg.det(invCovMat)) )
        sqrtDetCov = product(mv_CholeskyDiag) ! sqrt of determinant of the covariance matrix
        if (sqrtDetCov<=0) then
            write(output_unit,"(*(g0))") "sqrt of covariance determinant is <=0: ", sqrtDetCov
            write(output_unit,"(*(g0))") "Cholesky mv_DiagonalLogNormModel: "
            write(output_unit,"(*(g0))") mv_CholeskyDiag
            write(output_unit,"(*(g0))") "mv_CholeskyLowerLogNormModel/CovarianceMatrix: "
            write(output_unit,"(*(g20.13))") ((mv_CovMatUpperCholeskyLower(i,j),j=1,NDIM),new_line("A"),i=1,NDIM)
            stop
        end if
        sqrtDetCov2 = product(mv_CholeskyDiag2) ! sqrt of determinant of the covariance matrix
        if (sqrtDetCov2<=0) then
            write(output_unit,"(*(g0))") "sqrt of covariance determinant is <=0: ", sqrtDetCov2
            write(output_unit,"(*(g0))") "Cholesky mv_DiagonalLogNormModel: "
            write(output_unit,"(*(g0))") mv_CholeskyDiag2
            write(output_unit,"(*(g0))") "mv_CholeskyLowerLogNormModel/CovarianceMatrix: "
            write(output_unit,"(*(g20.13))") ((mv_CovMatUpperCholeskyLower2(i,j),j=1,NDIM),new_line("A"),i=1,NDIM)
            stop
        end if
        
        ! get the full Inverse covariance matricies
        mv_InvCovMat = getInvMatFromCholFac( nd = NDIM, CholeskyLower = mv_CovMatUpperCholeskyLower, Diagonal = mv_CholeskyDiag )
        mv_InvCovMat2 = getInvMatFromCholFac( nd = NDIM, CholeskyLower = mv_CovMatUpperCholeskyLower2, Diagonal = mv_CholeskyDiag2 )

        logFunc = COEF0 * NDATA
        do idata = 1, NDATA
            NormedPoint(1) = GRB%LogT90(idata) - mv_meanVec(1)
            NormedPoint(2) = GRB%LogEpk(idata) - mv_meanVec(2)
            NormedPoint2(1) = GRB%LogT90(idata) - mv_meanVec2(1)
            NormedPoint2(2) = GRB%LogEpk(idata) - mv_meanVec2(2)
            logFunc = logFunc + log(mv_amp * sqrtDetCov * exp(-0.5_RK * dot_product(NormedPoint,matmul(mv_invCovMat,NormedPoint))) &
                    + mv_amp2 * sqrtDetCov2 * exp(-0.5_RK * dot_product(NormedPoint2,matmul(mv_invCovMat2,NormedPoint2))) )
        end do

    end function getLogFunc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module LogFunc_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
