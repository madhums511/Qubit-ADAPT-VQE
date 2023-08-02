#!/bin/bash

#for file in *; do
for file in *v0*R1_00.erg; do

   #echo $file

   radius=${file/*R/}
   radius=${radius/\.*/}
   radius=`echo $radius | sed 's/_/\./g'`

   #echo $radius

   energy=`head -2 $file`
   #echo $energy

   echo "$radius $energy"

done

exit
