#!/bin/bash
STARTYR=95
ENDYR=101
STARTMO=1
ENDMO=12
loc=( 76S70W-76S0W 77S70W-77S0W 78S48W-78S30W 79S50W-79S30W 82S38W-73S38W 82S32W-73S32W ) 
loc1=76S70W-76S0W
loc2=77S70W-77S0W
loc3=78S48W-78S30W
loc4=79S50W-79S30W
loc5=82S38W-73S38W
loc6=82S32W-73S32W
var=( T S rho u v )
var1=T
var2=S
var3=rho
var4=u
var5=v
run=( ISMF-noEAIS ) #ISMF-noEAIS )
run1=ISMF
run2=ISMF-noEAIS
PREFIX=/global/homes/c/cbegeman/weddell_output/
SUFFIX=_limTrue.png

for i in $(seq -f "%04g" ${STARTYR} ${ENDYR} )
do
   for j in $( seq -f "%02g" ${STARTMO} ${ENDMO} )
   do
      k=$(($k + 1))
      cp ${PREFIX}ISMF_cmp_tau_frisEAcoast_0m_${i}-${j}${SUFFIX} ${PREFIX}ISMF_cmp_tau_frisEAcoast_${k}${SUFFIX}
   done
done
