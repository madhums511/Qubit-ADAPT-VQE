#!/bin/bash

#for file in *; do
for file in *vgo2*erg; do

   #echo $file

   radius=${file/*R/}
   radius=${radius/\.*/}
   radius=`echo $radius | sed 's/_/\./g'`

   #echo $radius
   energygr=`head -1 $file | awk '{print $6}'`
   energy=`head -50 $file|grep 'odd'|head -2  |tail -1| awk '{print $6}'`
   #echo $energy

   echo "$radius $energygr $energy"

done

exit
