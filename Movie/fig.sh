#!/bin/sh
gnuplot << EOF
set terminal epslatex standalone color colourtext

set output "$1.tex"

set size square

set pm3d map
#set pm3d interpolate 4,4

set palette defined ( 0 'blue', 1 'yellow', 2 'red' )
#set palette defined ( -10 'white', -4 'purple', 0 'black', 4 'orange', 10 'white' )
set cbrange[0:2]
pi = 3.14159265358979323846

set xlabel '\$ x\$'
set ylabel '\$ y\$'
set xrange[0:2*pi+0.01]
set yrange[0:2*pi+0.01]
#unset colorbox
splot "$1.dat" w pm3d notitle

set output
EOF



