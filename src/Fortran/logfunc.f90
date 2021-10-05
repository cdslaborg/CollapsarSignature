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
#include "Math_mod.f90"

module LogFunc_mod

    use iso_fortran_env, only: output_unit
    use Constants_mod, only: IK, RK, SQRT2PI

    implicit none

    integer(IK) , parameter     :: NDIM = 2_IK ! number of observational attributes: Dur, Epk.
    integer(IK) , parameter     :: NPAR = 11_IK ! number of parameters of the model.
   !integer(IK) , parameter     :: NPAR = 10_IK ! number of parameters of the model.
    integer(IK) , parameter     :: NDATA = 1966_IK ! Batse_orig data size
   !real(RK)    , parameter     :: COEF0 = NDIM * log( 1._RK / SQRT2PI )

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
        integer(IK)                 :: outputFileUnit
        integer(IK)                 :: idata

        ! Read L19 data from input data

        open( newunit = inputFileUnit &
            , file = "../../in/Batse_orig.csv" &
            , status = "old" &
            )

        open( newunit = outputFileUnit &
            , file = "Batse_orig.out" &
            , status = "replace" &
            )

        allocate(character(100) :: header)

        read(inputFileUnit,*) header; header = trim(adjustl(header))
        do idata = 1, NDATA
            read(inputFileUnit,*) GRB%T50(idata), GRB%T90(idata), GRB%Epk(idata)
            write(outputFileUnit,*) GRB%T50(idata), GRB%T90(idata), GRB%Epk(idata)
            GRB%LogT50(idata) = log(GRB%T50(idata))
            GRB%LogT90(idata) = log(GRB%T90(idata))
            GRB%LogEpk(idata) = log(GRB%Epk(idata))
        end do

        close(outputFileUnit)
        close(inputFileUnit)

    end subroutine readData

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the negative natural logarithm of MVN distribution evaluated at the input vector point.
    function getLogFunc(npar,param) result(logFunc)

        use Constants_mod, only: NEGINF_RK
        use Matrix_mod, only: getCholeskyFactor, getInvMatFromCholFac
        use Math_mod, only: getLogAddExp

        implicit none

        ! Function interface arguments.

        integer(IK) , intent(in)    :: npar
        real(RK)    , intent(in)    :: param(npar)

        ! Other variables for MVN mixture modeling.

        real(RK)                    :: MeanVecSGRB(NDIM)
        real(RK)                    :: stdDurSGRB
        real(RK)                    :: stdEpkSGRB
        real(RK)                    :: rhoDurEpkSGRB
        real(RK)                    :: logAmplitudeSGRB
        real(RK)                    :: MeanVecLGRB(NDIM)
        real(RK)                    :: stdDurLGRB
        real(RK)                    :: stdEpkLGRB
        real(RK)                    :: rhoDurEpkLGRB
        real(RK)                    :: logAmplitudeLGRB
        real(RK)                    :: ChoLowCovUppSGRB(NDIM,NDIM)
        real(RK)                    :: ChoDiaSGRB(NDIM)
        real(RK)                    :: ChoLowCovUppLGRB(NDIM,NDIM)
        real(RK)                    :: ChoDiaLGRB(NDIM)
        real(RK)                    :: InvCovMatSGRB(NDIM,NDIM)
        real(RK)                    :: InvCovMatLGRB(NDIM,NDIM)


        real(RK)                    :: logFunc
        real(RK)                    :: logProbDenSGRB
        real(RK)                    :: logProbDenLGRB
        integer(IK)                 :: idata, i, j
        real(RK)                    :: logSqrtDetInvCovSGRB
        real(RK)                    :: logSqrtDetInvCovLGRB
       !real(RK)                    :: invSqrtDetCov1
       !real(RK)                    :: invSqrtDetCov2
        real(RK)                    :: NormedDurEpkSGRB(NDIM)
        real(RK)                    :: NormedDurEpkLGRB(NDIM)

        MeanVecSGRB(1) = param(1)       ! SGRB: mean(log(T90))
        MeanVecSGRB(2) = param(2)       ! SGRB: mean(log(Epk))
        stdDurSGRB = exp(param(3))      ! SGRB: log(std(log(T90)))
        stdEpkSGRB = exp(param(4))      ! SGRB: log(std(log(Epk)))
        rhoDurEpkSGRB = tanh(param(5))  ! SGRB: FisherTrans(rho( log(T90), log(Epk) ))
                                        
        MeanVecLGRB(1) = param(6)       ! LGRB: mean(log(T90))
        MeanVecLGRB(2) = param(7)       ! LGRB: mean(log(Epk))
        stdDurLGRB = exp(param(8))      ! LGRB: log(std(log(T90)))
        stdEpkLGRB = exp(param(9))      ! LGRB: log(std(log(Epk)))
        rhoDurEpkLGRB = tanh(param(10)) ! LGRB: FisherTrans(rho( log(T90), log(Epk) ))

        ! Here we assume the amplitude of SGRBS MVN distribution is 1.
        ! Therefore,

        !logAmplitudeSGRB = 0._RK
        !logAmplitudeLGRB = log(2._RK)
        !logAmplitudeLGRB = param(11) + logAmplitudeSGRB
        logAmplitudeSGRB = log(0.5_RK + 0.5_RK * tanh(param(11)))
        logAmplitudeLGRB = log(1._RK - logAmplitudeSGRB)

        ChoLowCovUppSGRB(1,1) = stdDurSGRB**2
        ChoLowCovUppSGRB(2,2) = stdEpkSGRB**2
        ChoLowCovUppSGRB(1,2) = rhoDurEpkSGRB * stdDurSGRB * stdEpkSGRB
        ChoLowCovUppSGRB(2,1) = ChoLowCovUppSGRB(1,2)

        ChoLowCovUppLGRB(1,1) = stdDurLGRB**2
        ChoLowCovUppLGRB(2,2) = stdEpkLGRB**2
        ChoLowCovUppLGRB(1,2) = rhoDurEpkLGRB * stdDurLGRB * stdEpkLGRB
        ChoLowCovUppLGRB(2,1) = ChoLowCovUppLGRB(1,2)

        ! Compute the Cholescky factor of the covariance matrices of the two MVNs.

        call getCholeskyFactor(nd = NDIM, PosDefMat = ChoLowCovUppSGRB, Diagonal = ChoDiaSGRB)
        if (ChoDiaSGRB(1)<0._RK) then
            logFunc = NEGINF_RK
            error stop
            return
        end if

        call getCholeskyFactor(nd = NDIM, PosDefMat = ChoLowCovUppLGRB, Diagonal = ChoDiaLGRB)
        if (ChoDiaLGRB(1)<0._RK) then
            logFunc = NEGINF_RK
            error stop
            return
        end if

        ! The following is the log of the coefficient used in the definition of the MVN.

        logSqrtDetInvCovSGRB = -log(product(ChoDiaSGRB)) ! sqrt of determinant of the covariance matrix.
        !if (logSqrtDetInvCovSGRB<=0) then
        !    write(output_unit,"(*(g0))") "sqrt of covariance determinant is <=0: ", logSqrtDetInvCovSGRB
        !    write(output_unit,"(*(g0))") "Cholesky DiagonalLogNormModel: "
        !    write(output_unit,"(*(g0))") ChoDiaSGRB
        !    write(output_unit,"(*(g0))") "CholeskyLowerLogNormModel/CovarianceMatrix: "
        !    write(output_unit,"(*(g20.13))") ((ChoLowCovUppSGRB(i,j),j=1,NDIM),new_line("A"),i=1,NDIM)
        !    error stop
        !end if

        logSqrtDetInvCovLGRB = -log(product(ChoDiaLGRB)) ! sqrt of determinant of the covariance matrix
        !if (logSqrtDetInvCovLGRB<=0) then
        !    write(output_unit,"(*(g0))") "sqrt of covariance determinant is <=0: ", logSqrtDetInvCovLGRB
        !    write(output_unit,"(*(g0))") "Cholesky DiagonalLogNormModel: "
        !    write(output_unit,"(*(g0))") ChoDiaLGRB
        !    write(output_unit,"(*(g0))") "CholeskyLowerLogNormModel/CovarianceMatrix: "
        !    write(output_unit,"(*(g20.13))") ((ChoLowCovUppLGRB(i,j),j=1,NDIM),new_line("A"),i=1,NDIM)
        !    error stop
        !end if

        ! Get the full Inverse covariance matrices.

        InvCovMatSGRB = getInvMatFromCholFac( nd = NDIM, CholeskyLower = ChoLowCovUppSGRB, Diagonal = ChoDiaSGRB )
        InvCovMatLGRB = getInvMatFromCholFac( nd = NDIM, CholeskyLower = ChoLowCovUppLGRB, Diagonal = ChoDiaLGRB )

        logFunc = 0._RK
       !logFunc = COEF0 * NDATA
        do idata = 1, NDATA
            NormedDurEpkSGRB(1) = GRB%LogT90(idata) - MeanVecSGRB(1)
            NormedDurEpkSGRB(2) = GRB%LogEpk(idata) - MeanVecSGRB(2)
            NormedDurEpkLGRB(1) = GRB%LogT90(idata) - MeanVecLGRB(1)
            NormedDurEpkLGRB(2) = GRB%LogEpk(idata) - MeanVecLGRB(2)
            logProbDenSGRB = logAmplitudeSGRB + logSqrtDetInvCovSGRB - 0.5_RK * dot_product(NormedDurEpkSGRB,matmul(InvCovMatSGRB,NormedDurEpkSGRB))
            logProbDenLGRB = logAmplitudeLGRB + logSqrtDetInvCovLGRB - 0.5_RK * dot_product(NormedDurEpkLGRB,matmul(InvCovMatLGRB,NormedDurEpkLGRB))
            if (logProbDenSGRB >= logProbDenLGRB) then
                logFunc = logFunc + getLogAddExp(logValueLarger = logProbDenSGRB, logValueSmaller = logProbDenLGRB)
            else
                logFunc = logFunc + getLogAddExp(logValueLarger = logProbDenLGRB, logValueSmaller = logProbDenSGRB)
            endif
            !logFunc = logFunc + log(logAmplitudeSGRB * invSqrtDetCov1 * exp(-0.5_RK * dot_product(NormedDurEpkSGRB,matmul(InvCovMatSGRB,NormedDurEpkSGRB))) &
            !                      + logAmplitudeLGRB * invSqrtDetCov2 * exp(-0.5_RK * dot_product(NormedDurEpkLGRB,matmul(InvCovMatLGRB,NormedDurEpkLGRB))) )
        end do

    end function getLogFunc

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module LogFunc_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
