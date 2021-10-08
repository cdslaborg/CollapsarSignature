
!>  \brief
!>  This code fits a bivariate Gaussian mixture to the BATSE catalog data of logEpk and logT90.
program main

!use iso_fortran_env, only: IK => int32, RK => real64
use paramonte, only: runParaDRAM
use LogFunc_mod, only: getLogFunc, NPAR, readData
use Constants_mod, only: IK, RK

implicit none

! Read BATSE data and load it into the global variable `LogFunc_mod::GRB`

call readData()

! Run ParaDRAM MCMC simulation.

call runParaDRAM( ndim = NPAR &
                , getLogFunc = getLogFunc &
                , inputFile = "./paramonte.in" & ! this is optional argument
                )

end program main
