#!/bin/bash

STARTYR=70
ENDYR=100
STARTMO=1
ENDMO=12

#loc=( 76S70W-76S0W 77S70W-77S0W 78S48W-78S30W 79S50W-79S30W 82S38W-73S38W 82S32W-73S32W ) 
loc1=76S70W-76S0W
loc2=77S70W-77S0W
loc3=78S48W-78S30W
loc4=79S50W-79S30W
loc5=82S38W-73S38W
loc6=82S32W-73S32W
loc7=76S30W-60S30W
loc=${loc4} 
#loc=frisEAcoast_0m

var1=T
var2=S
var3=rho
var4=u
var5=v
#var=cmp_tau_quiver
var=( T S rho u v )
#var=( ${var2} ${var3} )
run1=ISMF
run2=ISMF-noEAIS
#run=( ISMF ) 
run=( ISMF-noEAIS )
#run=( ISMF_cmp ) #ISMF-noEAIS )

PREFIX=/global/homes/c/cbegeman/weddell_output/
#PREFIX=/global/homes/c/cbegeman/weddell_output/new_transects/
#SUFFIX=.png
SUFFIX=_limTrue.png

for r in "${run[@]}"
do
   for v in "${var[@]}"
   do
      for l in "${loc[@]}"
      do
         for i in $(seq -f "%04g" ${STARTYR} ${ENDYR} )
         do
            for j in $( seq -f "%02g" ${STARTMO} ${ENDMO} )
            do
               k=$(($k + 1))
               cp ${PREFIX}${r}_${v}_${l}_${i}-${j}${SUFFIX} ${PREFIX}temp/${r}_${v}_${l}_${k}${SUFFIX}
            done
         done
         ffmpeg -framerate 5 -i ../weddell_output/${r}_${v}_${l}_%d_limTrue.png -pix_fmt yuv420p ../weddell_output/${r}_${v}_${l}_70-100_limTrue.mov
         echo $k
         k=0
      done
   done
done

