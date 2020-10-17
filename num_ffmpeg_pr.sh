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
#loc=${loc3} 
loc=( 75S30W 74S30W )

var1=T
var2=S
var3=rho
var4=u
var5=v
var=( TS uv )
#var=( T S rho u v )

run1=ISMF
run2=ISMF-noEAIS
run=( ISMF ISMF-noEAIS )

PREFIX=/global/homes/c/cbegeman/weddell_output/
#PREFIX=/global/homes/c/cbegeman/weddell_output/new_transects/
SUFFIX=.png
#SUFFIX=_limTrue.png

for r in "${run[@]}"
do
   for v in "${var[@]}"
   do
      for l in "${loc[@]}"
      do
         for i in $(seq -f "%g" ${STARTYR} ${ENDYR} )
         do
           cp ${PREFIX}${r}_profiles_${v}_${l}_${i}-${i}${SUFFIX} ${PREFIX}${r}_profiles_${v}_${l}_${i}${SUFFIX}
         done
      done
   done
done
