#!/bin/bash

#for file in *; do
for file in *vgo2*erg; do

   #echo $file

   radius=${file/*R/}
   radius=${radius/\.*/}
   radius=`echo $radius | sed 's/_/\./g'`

   #echo $radius

   energy=`head -1 $file | awk '{print $6}'`

   #echo $energy

   echo "$radius $energy"

done

exit
