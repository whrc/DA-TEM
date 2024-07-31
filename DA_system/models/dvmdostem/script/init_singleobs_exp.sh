#!/bin/sh
#########################################################
# The purpose of this script:
#
# Prepare neccessary ensemble files for testing DART
#
# created by Chu-Chun Chang    May 2024
########################################################

# dvmdostem-DART directory
DART=~/DART/DA-TEM/DA_system/models/dvmdostem/work  # where the DART /work is 

# main directory for ensemble runs
DATADIR=/home/cchang/dvmdostem-workflows/bonanzacreeklter_ens40 

# EXP directory
EXPDIR=$DART/exp02


# ======== Parameters: 
nens=40           # number of ensemble

# ============ create WORKDIR if needed ================
if [ ! -d ${EXPDIR} ]
then
   mkdir -p ${EXPDIR}
   echo " working directory does not exist, mkdir : "${EXPDIR}
   cd ${EXPDIR}/
   mkdir obs
   mkdir truth
   mkdir bkgd
   mkdir anal
   mkdir log
fi


cd ${EXPDIR}/
# ============ create ensemble directories =====================
iens=1
while [ $iens -le $nens ]
do
      
   if [ $iens -lt 10 ] 
   then
       cens=0$iens
   else
       cens=$iens
   fi

   currentdir=ens${cens}
   
   # copy the model ensembles to EXPDIR    
   cp $DATADIR/${currentdir}/output/restart-eq.nc ${DART}/restart-eq.nc

   # run ./tem_to_dart to split leaf_C and stem_C, and generate model_bkgd.nc
   cd $DART/
   
   # copy restart file as a template for tem_to_dart output
   if [ -f $DART/model_bkgd.nc ]
   then
      rm $DART/model_bkgd.nc	   
   fi
   cp $DART/restart-eq.nc $DART/model_bkgd.nc

   echo " Start running tem_to_dart for ensemble" $cens
   ./tem_to_dart > tem_to_dart_ens${cens}.log &
   pid=$!

   # Wait for the background process to complete
    wait $pid

   # move and store data in EXPDIR	   
   mv $DART/model_bkgd.nc ${EXPDIR}/bkgd/model_bkgd_${cens}.nc
   mv $DART/restart-eq.nc ${EXPDIR}/restart-eq-${cens}.nc
   mv $DART/tem_to_dart_ens${cens}.log ${EXPDIR}/log/tem_to_dart_ens${cens}.log      

   echo "ens $iens completed!! "
     
   iens=`expr $iens + 1`
   
done
