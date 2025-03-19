#!/bin/bash

for (( i=100 ; i <= 200 ; i++ ))
do
./fig_png.sh "Intensity.00$i"

#mv "w.00$i-inc.eps" "w.00$i-inc_copy.eps"
convert -density 300x300 "w.00$i-inc.eps" EPS2:"w.00$i-inc"

latex "w.00$i.tex"

dvips "w.00$i.dvi"
ps2eps -r 300x300 -f "w.00$i.ps"
convert -density 300x300 "w.00$i.eps" "w.00$i.png"
rm -rf *.eps *inc *.ps *.tex *.aux *.log



#declare -i n
#n=${#i}
#echo "n = $n"
#if [ $n -eq 1 ]  ; then
#mv "./w@time=$i.png" "./w@time=000$i.png"
#elif [ $n -eq 2 ] ; then
#mv "./w@time=$i.png" "./w@time=00$i.png"
#elif [ $n -eq 3 ] ; then
#mv "./w@time=$i.png" "./w@time=0$i.png"

#fi



echo "$i"
done

./mencoder mf://*.png -mf w=50:h=50:fps=5:type=png -ovc lavc \-lavcopts vcodec=msmpeg4:mbd=2:trell -oac copy -o Movie.avi

rm *eps *inc *ps
