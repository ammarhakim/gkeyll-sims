2x2v GK ETG nonlinear cases with R/LT = 4, 6, 8, 10. kinetic electrons, adiabatic ions.

> pgkyl -f etg4d-rlt4_phi2_ -f etg4d-rlt6_phi2_ -f etg4d-rlt8_phi2_ -f etg4d-rlt10_phi2_ log plot -f0

the resulting plot is saved as etg4d_nonlinear_saturation.png

note: R/LT=10 case required higher resolution and larger box size. ran on 48 cores on edison in ~30min (batch file included).
