#/bin/sh

WORKDIR=~/DART/DA-TEM/DA_system/models/dvmdostem/work
EXPDIR=$WORKDIR/exp02
nens=40
var=biomass_root
ANALDIR=$EXPDIR/anal_${var}

echo "ANALDIR = "$ANALDIR

if [ ! -d $ANALDIR ]
then
   mv $EXPDIR/anal $ANALDIR	
   mkdir $ANALDIR/mean
   mkdir $ANALDIR/analysis
   mkdir $ANALDIR/preassim   
fi

# move mean state and single files
mv $WORKDIR/analysis_mean.nc $ANALDIR/mean/analysis_mean.nc
mv $WORKDIR/preassim_mean.nc $ANALDIR/mean/preassim_mean.nc
mv $WORKDIR/output_mean.nc $ANALDIR/mean/output_mean.nc
mv $WORKDIR/output_priorinf_mean.nc $ANALDIR/mean/output_priorinf_mean.nc
mv $WORKDIR/preassim_priorinf_mean.nc $ANALDIR/mean/preassim_priorinf_mean.nc
mv $WORKDIR/analysis_priorinf_mean.nc $ANALDIR/mean/analysis_priorinf_mean.nc
mv $WORKDIR/analysis_sd.nc $ANALDIR/mean/analysis_sd.nc
mv $WORKDIR/preassim_sd.nc $ANALDIR/mean/preassim_sd.nc
mv $WORKDIR/output_sd.nc $ANALDIR/mean/output_sd.nc

# move obs
mv $WORKDIR/obs_seq.out $EXPDIR/obs/obs_seq_${var}.out
mv $WORKDIR/obs_seq.final $EXPDIR/obs/obs_seq_${var}.final
#mv $WORKDIR/set_def.out $EXPDIR/obs/set_def_${var}.out
mv $WORKDIR/obs_seq.in $EXPDIR/obs/obs_seq_${var}.in

# Loop over ensembles
iens=1
while [ $iens -le $nens ]
do

   if [ $iens -lt 10 ]
   then
       cens=0$iens
   else
       cens=$iens
   fi

   mv $WORKDIR/preassim_member_00${cens}.nc $ANALDIR/preassim/preassim_member_00${cens}.nc
   mv $WORKDIR/analysis_member_00${cens}.nc $ANALDIR/analysis/analysis_member_00${cens}.nc

   echo " ============== completed coping data for ens =  "$cens

   iens=`expr $iens + 1`
done


