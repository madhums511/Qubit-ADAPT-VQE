#!/bin/bash

#for file in *; do
for file in *v08f*dat; do

   #echo $file

   radius=${file/*R/}
   radius=${radius/\.*/}
   radius=`echo $radius | sed 's/_/\./g'`

   #echo $radius

   energy=`head -2 $file |tail -1`
   #echo $energy

   echo "$radius $energy"

done

exit
