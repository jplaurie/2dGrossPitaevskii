reset

set contour base
set cntrparam levels discrete 0
unset surface
set table "real=0.dat"

splot "./Intensity.00200.dat" u 2:3:4 
unset table

set contour base
set cntrparam levels discrete 0
unset surface
set table "imag=0.dat"
splot "./Intensity.00200.dat" u 2:3:5
unset table
