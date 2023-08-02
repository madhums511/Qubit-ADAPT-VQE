#!/bin/tcsh
#
#     HHSP02.CSH
#     ==========
#
#     Perform an H\bar{H} "standard" calculation.
#      
#     ****   Non-symmetry adapted basis functions !!!   ****
#
#     ****   Perform a symmetry analysis.               ****
#
#     ****   Reduced output compared to HHSP01.CSH !!!  **** 
#
  set echo
#
#
#    Define the compiler options:
#    ----------------------------
#
LINUX:
  set FCMP="$FORTRAN_COMPILER $FORTRAN_TEST_OPT $FORTRAN_COMP_OPT"
  set FCMP77="$FORTRAN_COMPILER $FORTRAN_TEST_OPT $FORTRAN_COMP_OPT"
#
#
#    Define some paths:
#    ------------------
#
  set FORT=$AMOHOME/H-H/hh
  set INP=$AMOHOME/H-H/hh/input
  set OUT=$AMOHOME/H-H/hh/out
  set SCR=$AMOHOME/H-H/hh/scr


#
#
if ($#argv <> 5) then
   echo " Incorrect number of arguments ==> Please enter 5 arguments"
   echo "     1. Name of the molecule (hh, h2, het,...) "
   echo "     2. Number of basis functions (b501, b868, 400,...) "
   echo "     3. Basis set version (v01, v05a, v2f,...) "
   echo "     4. First part of internuclear distance (-1- for 1.40) "
   echo "     5. Second part of internuclear distance (-40- for 1.40) "
   exit
endif
#
set MOL=$1
set BAS1=$2
set BAS2=$3
set R1VAL=$4
set R2VAL=$5
#
set RVAL="R""$R1VAL""_""$R2VAL"
set BASIS="$BAS1""$BAS2"
set MOLBAS1="$MOL""$BAS1"
#
set HHBAS="$MOL""$BASIS""$RVAL"
#
#
#    Check, if the input files exist:
#    --------------------------------
#
  if -e $INP/$HHBAS.dat then
    goto endcheck
  else 
    echo " Input file $INP/$HHBAS.dat doesn't exist ! "
    exit
  endif
#
endcheck:
#
#
#
#        Cleanup from previous (interrupted) runs:
#        -----------------------------------------
#
  if -e $SCR/$HHBAS.f10  rm $SCR/$HHBAS.f10
  if -e $SCR/$HHBAS.evv  rm $SCR/$HHBAS.evv
  if -e $SCR/$HHBAS.hmt  rm $SCR/$HHBAS.hmt
  if -e $SCR/$HHBAS.smt  rm $SCR/$HHBAS.smt
  if -e $SCR/$HHBAS.vmt  rm $SCR/$HHBAS.vmt
  if -e $SCR/$HHBAS.qmt  rm $SCR/$HHBAS.qmt
  if -e $SCR/$HHBAS.q    rm $SCR/$HHBAS.q
#
  if -e $SCR/$HHBAS.vir    rm $SCR/$HHBAS.vir
  if -e $SCR/$HHBAS.h2err  rm $SCR/$HHBAS.h2err
#
#
#    Compile the H\bar{H} program
#    ----------------------------
#
#        Compiler-dependent definitions: 
#        -------------------------------
#
# AIX
#  f77 -o vorh2u_$HHBAS.x vorh2u.f
#
  $FCMP -o vorh2u_$HHBAS.x vorh2u.f
#
#
  setenv INPDAT $INP/$HHBAS.dat
  setenv OUTDAT h2upar.for
  setenv OUT2DAT h2upar2.for
#
  vorh2u_$HHBAS.x
#
#
#    Compiler-dependent definitions: 
#    -------------------------------
#
# AIX
#  f77 -o hhsp02_$HHBAS.x hhsp02.f
#  f77 -o hhsp01_q_$HHBAS.x hhsp01_q.f
#  f77 -o hhqsym_$HHBAS.x hhqsym.f
###
#
# Intel/Linux
  $FCMP -o hhsp02_$HHBAS.x hhsp02.f
  $FCMP -o hhsp01_q_$HHBAS.x hhsp01_q.f
  $FCMP -o hhqsym_$HHBAS.x hhqsym.f
#
#
  rm h2upar.for
  rm h2upar2.for
#
#
#      Start the H\bar{H} program:
#      ---------------------------
#
     setenv VIRIAL $SCR/$HHBAS.vir
     setenv H2ERR  $SCR/$HHBAS.h2err
     setenv EIGVEC $SCR/$HHBAS.evv
     setenv OVLAP  $SCR/$HHBAS.smt
     setenv QMAT   $SCR/$HHBAS.qmt
     setenv VMAT   $SCR/$HHBAS.vmt
     setenv QHQMAT $SCR/$HHBAS.hmt
#
     setenv F10DAT $SCR/$HHBAS.f10
     setenv F09DAT $SCR/$HHBAS.f09
     setenv F18DAT $SCR/$HHBAS.f18
     setenv F23DAT $SCR/$HHBAS.f23
#
     hhsp02_$HHBAS.x < $INP/"$HHBAS".dat
#
#
#      Cleanup:
#      --------
#
     if -e $SCR/$HHBAS.vir    rm $SCR/$HHBAS.vir
     if -e $SCR/$HHBAS.h2err  rm $SCR/$HHBAS.h2err
     if -e $SCR/$HHBAS.f10    rm $SCR/$HHBAS.f10
     if -e $SCR/$HHBAS.f09    rm $SCR/$HHBAS.f09
     if -e $SCR/$HHBAS.f18    rm $SCR/$HHBAS.f18
     if -e $SCR/$HHBAS.f23    rm $SCR/$HHBAS.f23
     #if -e $SCR/$HHBAS.hmt    rm $SCR/$HHBAS.hmt
     if -e $SCR/$HHBAS.smt    rm $SCR/$HHBAS.smt
     if -e $SCR/$HHBAS.vmt    rm $SCR/$HHBAS.vmt
     if -e $SCR/$HHBAS.qmt    rm $SCR/$HHBAS.qmt
#
#
#
#       Calculate the q-symmetry matrix elements:
#       -----------------------------------------
#
     setenv F10DAT $SCR/$HHBAS.f10
     setenv F09DAT $SCR/$HHBAS.f09
     setenv F18DAT $SCR/$HHBAS.f18
     setenv F23DAT $SCR/$HHBAS.f23
#
     setenv QSYM   $SCR/$HHBAS.q
#
     hhsp01_q_$HHBAS.x < $INP/$HHBAS.dat
#
#
#        Cleanup:
#        --------
#
     if -e $SCR/$HHBAS.f10  rm $SCR/$HHBAS.f10
     if -e $SCR/$HHBAS.f09  rm $SCR/$HHBAS.f09
     if -e $SCR/$HHBAS.f18  rm $SCR/$HHBAS.f18
     if -e $SCR/$HHBAS.f23  rm $SCR/$HHBAS.f23
#
#
#       Check the symmetries using program HHQSYM:
#       ------------------------------------------
#
  setenv XMTDAT $SCR/$HHBAS.evv
  setenv MMRDAT $SCR/$HHBAS.q
#
  setenv OUTDAT $OUT/$HHBAS.erg
#
  hhqsym_$HHBAS.x
#
#
#      Cleanup:
#      --------
#
  rm $SCR/$HHBAS.evv
  rm $SCR/$HHBAS.q
#
#
  rm hhsp02_$HHBAS.x
  rm vorh2u_$HHBAS.x
  rm hhqsym_$HHBAS.x
  rm hhsp01_q_$HHBAS.x
#
 unset echo
 date
 time
