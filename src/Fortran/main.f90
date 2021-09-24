!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%  
!%  Description:
!%      +   Run the Monte Carlo sampler of the ParaMonte library given the input log-target density function `getLogFunc()`.
!%  Output:
!%      +   The simulation output files.
!%  Author:
!%      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
!%  Visit:
!%      +   https://www.cdslab.org/paramonte
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main

!use iso_fortran_env, only: IK => int32, RK => real64
use paramonte, only: runParaDRAM
use LogFunc_mod, only: getLogFunc, NPAR, readData
use Constants_mod, only: IK, RK

implicit none

! read L19 data

call readData()

! run ParaDRAM MCMC simulation

call runParaDRAM( ndim = NPAR &
                , getLogFunc = getLogFunc &
                , inputFile = "./paramonte.in" & ! this is optional argument
                )

end program main
