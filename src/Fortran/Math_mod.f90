module Math_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Math_mod"

    interface getLogAddExp
        module procedure :: getLogAddExp_RK, getLogSumExp_RK
    end interface getLogAddExp

contains

    PURE function getLogAddExp_RK(logValueLarger, logValueSmaller) result(logSumExp)
        use Constants_mod, only: RK, LOGTINY_RK
        real(RK)    , intent(in)    :: logValueLarger
        real(RK)    , intent(in)    :: logValueSmaller
        real(RK)                    :: logSumExp
        real(RK)                    :: logRatio
        block
        use iso_fortran_env, only: output_unit
        if (logValueSmaller > logValueLarger) then
            write(output_unit,"(*(g0,:,' '))") "logValueSmaller > logValueLarger:", logValueSmaller, logValueLarger
            error stop
        end if
        end block
        logRatio = logValueSmaller - logValueLarger
        logSumExp = logValueLarger
        if (logRatio > LOGTINY_RK) logSumExp = logSumExp + log(1._RK + exp(logRatio))
    end function getLogAddExp_RK

    pure function getLogSumExp_RK(LogValue, maxLogValue) result(logSumExp)
#if INTEL_COMPILER_ENABLED && DLL_ENABLED && (WINDOWS_ENABLED || DARWIN_ENABLED)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExp_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        real(RK)    , intent(in), contiguous    :: LogValue(:)
        real(RK)    , intent(in)                :: maxLogValue
        real(RK)                                :: logSumExp
        real(RK)                                :: logRatio
        integer(IK)                             :: i
        logSumExp = 0._RK
        do i = 1, size(LogValue)
            logRatio = LogValue(i) - maxLogValue
            if (logRatio > LOGTINY_RK) then
                logSumExp = logSumExp + exp(logRatio)
            end if
        end do
        logSumExp = maxLogValue + log(logSumExp)
    end function getLogSumExp_RK

end module Math_mod