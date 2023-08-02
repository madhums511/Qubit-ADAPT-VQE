#!/bin/bash

#for file in *; do
for file in *v09*erg; do

   #echo $file

   radius=${file/*R/}
   radius=${radius/\.*/}
   radius=`echo $radius | sed 's/_/\./g'`

   #echo $radius

   energy=`head -8 $file`

   number1=`head -1 $file|tail -1| awk '{print $NF}'`
   number1=${number1/\)*/}
   number1=${number1/*-/}
   number2=`head -2 $file|tail -1| awk '{print $NF}'`
   number2=${number2/\)*/}
   number2=${number2/*-/}
   number3=`head -3 $file|tail -1| awk '{print $NF}'`
   number3=${number3/\)*/}
   number3=${number3/*-/}
   number4=`head -4 $file|tail -1| awk '{print $NF}'`
   number4=${number4/\)*/}
   number4=${number4/*-/}
   number5=`head -5 $file|tail -1| awk '{print $NF}'`
   number5=${number5/\)*/}
   number5=${number5/*-/}
   number6=`head -6 $file|tail -1| awk '{print $NF}'`
   number6=${number6/\)*/}
   number6=${number6/*-/}
   number7=`head -7 $file|tail -1| awk '{print $NF}'`
   number7=${number7/\)*/}
   number7=${number7/*-/}
   number8=`head -8 $file|tail -1| awk '{print $NF}'`
   number8=${number8/\)*/}
   number8=${number8/*-/}
   
   schnitt=`echo "scale=2; 8 - $number1 - $number2 - $number3 - $number4 - $number5 - $number6 - $number7 - $number8 " | bc `
   #echo $energy

   echo "$radius $energy"
   echo "$schnitt"
done

exit
