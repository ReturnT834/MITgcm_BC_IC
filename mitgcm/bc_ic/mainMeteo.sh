#! /bin/bash


# Step 0. Get Mask info
# Step 1. Space Interpolation

RUNDATE=ECMWF-ERAINT
#DATESTART=$( date -d "${RUNDATE}  -  7  days " +%Y%m%d-%H:%M:%S )
#DATE__END=$( date -d "${RUNDATE}  + 36 hours " +%Y%m%d-%H:%M:%S )
#DOWNSTART=$( date -d "${RUNDATE}  -  8  days " +%Y%m%d-%H:%M:%S )
DATESTART=20000101-00:30:00
#DATE__END=20000101-23:30:00
DATE__END=20101231-23:30:00
. ./profile.inc

RUNDIR=/g100_scratch/userexternal/plazzari/RETURN_P-o-C/${RUNDATE}
         MASK_ARSO=$RUNDIR/mask.arso.nc
         MASK_RCM_V0=$RUNDIR/mask.rcm_v0.nc
          MASK_128=$RUNDIR/mask128.nc          # Cadeau mask
          BATHYGZ=bathy.gz
            BATHY=$RUNDIR/bathy.nc
DOWNLOADED_METEO=$RUNDIR/DOWNLOAD
      ORIG_METEO=/g100_scratch/userexternal/plazzari/REGCM/ECMWF-ERAINT/
        BC_METEO=$RUNDIR/BC/
   FTP_BATCHFILE=$RUNDIR/batchfile
DOWNTIMELISTFILE=$RUNDIR/to_download_timelist.txt
    TIMELISTFILE=$RUNDIR/meteo_timelist.txt
export PATH=$PATH:/share/scratch/backup_root/usr/bin/ # to have ftp on fluxus
export PATH=$PATH:/marconi/home/usera07ogs/a07ogs00/OPA/V3C/HOST/marconi/bin/

mkdir -p $RUNDIR
mkdir -p $DOWNLOADED_METEO
mkdir -p $ORIG_METEO


### Step 0.  GET MASK INFO  #####################################

medmit_prex_or_die "gzip -cd static-data/masks/METEO/sftlf_ALP-3_MOHC-HadGEM2-ES_historical_r0i0p0_ICTP-RegCM4-7_fpsconv-x2yn2-v1_fx.nc.gz > $MASK_RCM_V0 "
medmit_prex_or_die "gzip -cd static-data/masks/CADEAU/${BATHYGZ} > $BATHY   "
medmit_prex_or_die "python static-data/masks/CADEAU/maskgen.py -b $BATHY -o $MASK_128  "

### Step 1. Space Interpolation   ###############################

medmit_prex_or_die " python TimeList_generator.py -s $DATESTART -e $DATE__END --hours 1 > $TIMELISTFILE "
medmit_prex_or_die " mpirun -np 7 python meteo_generator_RegCM_LOW_RAM.py -i $ORIG_METEO -o $BC_METEO -m $MASK_128 --nativemask $MASK_RCM_V0 -t $TIMELISTFILE "

##################################################################

medmit_prex_or_die " mv $BC_METEO/CHECK $RUNDIR "

