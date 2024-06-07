#!/bin/bash

for (( i=0 ; i < 300 ; i++ ))

do

filenumber=$(printf "%05d" $i)
./fig_png.sh "psi.$filenumber"

echo "$i"
done

#ffmpeg -f image2 -pattern_type glob -framerate 12 -i 'psi.*.png' -preset veryslow -crf 25 -vcodec libx264     Movie.mp4

ffmpeg -f image2 -pattern_type glob -framerate 12 -i 'psi.*.png' -preset veryslow -vcodec libx265 -x265-params crf=15:psy-rd=1 output.mp4




#rm -rf *.png
