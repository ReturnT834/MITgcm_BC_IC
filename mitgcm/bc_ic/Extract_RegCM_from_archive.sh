SIMULATION=ECMWF-ERAINT


YEAR=$1


ORIG_METEO=/gss/gss_work/DRES_OGS_BiGe/plazzari/ICT23_ESP/${SIMULATION}/evaluation/r1i1p1/ICTP-RegCM4-7/fpsconv-x2yn2-v1/1hr/

SCRATCH_METEO=/g100_scratch/userexternal/plazzari/REGCM/${SIMULATION}

echo ${SCRATCH_METEO}
for var in tas huss uas vas rsds rlus pr; do 
   mkdir -p ${SCRATCH_METEO}/${var}
   rsync -PravzHS ${ORIG_METEO}/${var}/${var}*${YEAR}*nc ${SCRATCH_METEO}/${var}/
done
