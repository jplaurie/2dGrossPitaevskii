#!/bin/sh
gnuplot << EOF
set terminal png transparent truecolor
set output "$1.png"

set size square

set pm3d map

set palette defined ( -2 'blue', 0 'yellow', 2 'red' )
#set palette defined ( -2 'white', -1 'purple', 0 'black', 1 'orange', 2 'white' )
set cbrange[-10:10]
pi = 3.14159265358979323846

set xlabel ''
set ylabel ''
unset colorbox
set xtics ''
set ytics ''

set xrange[0:32*pi]
set yrange[0:32*pi]
#unset colorbox
#splot "$1.dat" u 1:2:3 w pm3d notitle
plot "../Velocity/$1" u 1:2:3 w image notitle

set output
EOF



