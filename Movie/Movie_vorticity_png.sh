#!/bin/bash

for (( i=700 ; i <750 ; i++ ))

do

filenumber=$(printf "%05d" $i)
./fig_vorticity_png.sh "W.$filenumber"

echo "$i"
done

#./mencoder mf://*.png -mf w=50:h=50:fps=25:type=png -ovc lavc \-lavcopts vcodec=msmpeg4:mbd=2:trell -oac copy -o Movie.avi

#./mencoder mf://*.png -mf fps=25 -o Movie.avi -ovc xvid -xvidencopts pass=2:bitrate=3000 -vf scale=640:400,eq2=1.5

./mencoder mf://*.png -mf w=50:h=50:fps=25:type=png -ovc xvid -xvidencopts bitrate=3000 -oac copy -o Vorticity.avi

#./mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc raw -oac copy -o output.avi
rm -rf *.png
