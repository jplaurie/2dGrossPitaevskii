#!/bin/sh
gnuplot << EOF
set terminal png size 1000,1000  transparent truecolor
set output "$1.png"

set size square

set pm3d map

set palette defined ( -2 'blue', 0 'yellow', 2 'red' )
#set palette defined ( -2 'white', -1 'purple', 0 'black', 1 'orange', 2 'white' )
set cbrange[0:1]
pi = 3.14159265358979323846

set xlabel ''
set ylabel ''
unset colorbox
#set xtics ''
#set ytics ''
unset xtics
unset ytics


set xrange[0:2*pi]
set yrange[0:2*pi]
#unset colorbox
#splot "$1.dat" u 1:2:3 w pm3d notitle
plot "../data/$1" u 1:2:((\$3**2.0+\$4**2.0)**0.5) w image notitle

set output
EOF



