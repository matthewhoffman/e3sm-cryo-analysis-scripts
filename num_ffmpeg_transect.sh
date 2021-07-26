#!/bin/bash
STARTYR=70
ENDYR=100
STARTMO=1
ENDMO=12
PREFIX=/global/homes/c/cbegeman/weddell_output/
SUFFIX=.png

for i in $(seq -f "%04g" ${STARTYR} ${ENDYR} )
do
   for j in $( seq -f "%02g" ${STARTMO} ${ENDMO} )
   do
      k=$(($k + 1))
      cp ${PREFIX}ISMF_Tsigma1_trough_crossshelf_${i}-${j}_zmax1000${SUFFIX} ${PREFIX}ISMF_Tsigma1_trough_crossshelf_${k}${SUFFIX}
   done
done
ffmpeg -framerate 5 -i ${PREFIX}ISMF_Tsigma1_trough_crossshelf_%d${SUFFIX} -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ${PREFIX}ISMF_Tsigma1_trough_crossshelf${STARTYR}-${ENDYR}_zmax1000.mov
