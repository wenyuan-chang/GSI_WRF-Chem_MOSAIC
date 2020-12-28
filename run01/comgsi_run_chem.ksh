#!/bin/ksh
#####################################################
# machine set up (users should change this part)
#####################################################

set -x
#
# GSIPROC = processor number used for GSI analysis
#------------------------------------------------
  GSIPROC=1
  ARCH='LINUX'

# Supported configurations:
            # IBM_LSF,
            # LINUX, LINUX_LSF, LINUX_PBS,
            # DARWIN_PGI
# this script can run 4 GSI chem cases
#  1. WRF-Chem GOCART with MODIS AOD observation 
#         bk_core=WRFCHEM_GOCART and obs_type=MODISAOD
#         background= wrfinput_enkf_d01_2012-06-03_18:00:00
#         observations=Aqua_Terra_AOD_BUFR:2012-06-03_00:00:00
#  2. WRF-Chem GOCART with PM25 observation 
#         bk_core=WRFCHEM_GOCART and obs_type=PM25
#         background= wrfinput_enkf_d01_2012-06-03_18:00:00
#         observations=anow.2012060318.bufr
#  3. WRF-Chem PM25 with MP25 observation
#         bk_core=WRFCHEM_PM25 and obs_type=PM25
#         background= wrfinput_enkf_d01_2012-06-03_18:00:00
#         observations=anow.2012060318.bufr
#  4. CMAQ with MP25 observation 
#         bk_core=CMAQ and obs_type=PM25
#         background= cmaq2gsi_4.7_20130621_120000.bin
#         observations=anow.2013062112.bufr
#
#####################################################
# case set up (users should change this part)
#####################################################
#
# ANAL_TIME= analysis time  (YYYYMMDDHH)
# WORK_ROOT= working directory, where GSI runs
# PREPBURF = path of PreBUFR conventional obs
# BK_FILE  = path and name of background file
# OBS_ROOT = path of observations files
# FIX_ROOT = path of fix files
# GSI_EXE  = path and name of the gsi executable 
  ANAL_TIME=2019041818
  JOB_DIR=/public/lapc/changwy/GSI/run01
     #normally you put run scripts here and submit jobs form here, require a copy of gsi.x at this directory
  RUN_NAME=example
  OBS_ROOT=/public/lapc/changwy/GSI/run01/INPUT
  BK_ROOT=/public/lapc/changwy/GSI/run01/INPUT
  GSI_ROOT=/public/lapc/changwy/GSI/comGSIv3.7_EnKFv1.3_cwy1.1_20201123
  CRTM_ROOT=/public/lapc/changwy/GSI/comGSIv3.7_EnKFv1.3_cwy1.1_20201123/CRTM_v2.3.0
  GSI_EXE=${JOB_DIR}/gsi.x  #assume you have a copy of gsi.x here
  WORK_ROOT=${JOB_DIR}/${RUN_NAME}
  FIX_ROOT=${GSI_ROOT}/fix
  GSI_NAMELIST=/public/lapc/changwy/GSI/run01/comgsi_namelist_chem.sh
### -------------------------------------- cwy step 1
  PREPBUFR1=${OBS_ROOT}/pm10bufr
  PREPBUFR2=${OBS_ROOT}/pm25bufr
  PREPBUFR3=${OBS_ROOT}/modisbufr
  PREPBUFR4=${OBS_ROOT}/ce318bufr
  PREPBUFR5=${OBS_ROOT}/srfeabsbufr
  PREPBUFR6=${OBS_ROOT}/srfescabufr
  PREPBUFR7=${OBS_ROOT}/srfeextbufr

  BK_FILE=${BK_ROOT}/wrfout.4bin.kashi/RESULT_KASHI_NODA_2010EMIS/wrfout_d02_2019-04-18_18:00:00
### -------------------------------------- cwy step 1
#
#------------------------------------------------
# bk_core= set background (WRFCHEM_GOCART WRFCHEM_PM25 or CMAQ)
# obs_type= set observation type (MODISAOD or PM25)
# if_clean = clean  : delete temperal files in working directory (default)
#            no     : leave running directory as is (this is for debug only)
### -------------------------------------- cwy step 2
  bk_core=WRFCHEM_MOSAIC
#  bk_core=WRFCHEM_GOCART
#  bk_core=WRFCHEM_PM25
#  obs_type=PM25
### -------------------------------------- cwy step 2
  if_clean=no
#
#
#####################################################
# Users should NOT make changes after this point
#####################################################
#
BYTE_ORDER=Big_Endian
# BYTE_ORDER=Little_Endian

case $ARCH in
   'IBM_LSF')
      ###### IBM LSF (Load Sharing Facility)
      RUN_COMMAND="mpirun.lsf " ;;

   'LINUX')
      if [ $GSIPROC = 1 ]; then
         #### Linux workstation - single processor
         RUN_COMMAND=""
      else
         ###### Linux workstation -  mpi run
        RUN_COMMAND="mpirun -np ${GSIPROC} "
      fi ;;

   'LINUX_LSF')
      ###### LINUX LSF (Load Sharing Facility)
      RUN_COMMAND="mpirun.lsf " ;;

   'LINUX_PBS')
      #### Linux cluster PBS (Portable Batch System)
      RUN_COMMAND="mpirun -np ${GSIPROC} " ;;

   'DARWIN_PGI')
      ### Mac - mpi run
      if [ $GSIPROC = 1 ]; then
         #### Mac workstation - single processor
         RUN_COMMAND=""
      else
         ###### Mac workstation -  mpi run
         RUN_COMMAND="mpirun -np ${GSIPROC} -machinefile ~/mach "
      fi ;;

   * )
     print "error: $ARCH is not a supported platform configuration."
     exit 1 ;;
esac


##################################################################################
# Check GSI needed environment variables are defined and exist
#
 
# Make sure ANAL_TIME is defined and in the correct format
if [ ! "${ANAL_TIME}" ]; then
  echo "ERROR: \$ANAL_TIME is not defined!"
  exit 1
fi

# Make sure WORK_ROOT is defined and exists
if [ ! "${WORK_ROOT}" ]; then
  echo "ERROR: \$WORK_ROOT is not defined!"
  exit 1
fi

# Make sure the background file exists
if [ ! -r "${BK_FILE}" ]; then
  echo "ERROR: ${BK_FILE} does not exist!"
  exit 1
fi

# Make sure OBS_ROOT is defined and exists
if [ ! "${OBS_ROOT}" ]; then
  echo "ERROR: \$OBS_ROOT is not defined!"
  exit 1
fi
if [ ! -d "${OBS_ROOT}" ]; then
  echo "ERROR: OBS_ROOT directory '${OBS_ROOT}' does not exist!"
  exit 1
fi

# Set the path to the GSI static files
if [ ! "${FIX_ROOT}" ]; then
  echo "ERROR: \$FIX_ROOT is not defined!"
  exit 1
fi
if [ ! -d "${FIX_ROOT}" ]; then
  echo "ERROR: fix directory '${FIX_ROOT}' does not exist!"
  exit 1
fi

# Set the path to the CRTM coefficients 
if [ ! "${CRTM_ROOT}" ]; then
  echo "ERROR: \$CRTM_ROOT is not defined!"
  exit 1
fi
if [ ! -d "${CRTM_ROOT}" ]; then
  echo "ERROR: fix directory '${CRTM_ROOT}' does not exist!"
  exit 1
fi


# Make sure the GSI executable exists
if [ ! -x "${GSI_EXE}" ]; then
  echo "ERROR: ${GSI_EXE} does not exist!"
  exit 1
fi

# Check to make sure the number of processors for running GSI was specified
if [ -z "${GSIPROC}" ]; then
  echo "ERROR: The variable $GSIPROC must be set to contain the number of processors to run GSI"
  exit 1
fi

#
##################################################################################
# Create the ram work directory and cd into it

workdir=${WORK_ROOT}
echo " Create working directory:" ${workdir}

if [ -d "${workdir}" ]; then
  rm -rf ${workdir}
fi
mkdir -p ${workdir}
cd ${workdir}

#
##################################################################################

echo " Copy GSI executable, background file, and link observation bufr to working directory"

# Save a copy of the GSI executable in the workdir
cp ${GSI_EXE} gsi.x

# Bring over background field (it's modified by GSI so we can't link to it)

if [ ${bk_core} = WRFCHEM_MOSAIC ] ; then
   cp ${BK_FILE} ./wrf_inout
fi
if [ ${bk_core} = WRFCHEM_GOCART ] ; then
   cp ${BK_FILE} ./wrf_inout
fi
if [ ${bk_core} = WRFCHEM_PM25 ] ; then
   cp ${BK_FILE} ./wrf_inout
fi
if [ ${bk_core} = CMAQ ] ; then
   cp ${BK_FILE} ./cmaq_in.bin
fi

# Link to the observation data
#if [ ${obs_type} = MODISAOD ] ; then
 ln -s ${PREPBUFR3} ./modisbufr
#fi
#if [ ${obs_type} = PM25 ] ; then
 ln -s ${PREPBUFR2} ./pm25bufr
 ln -s ${PREPBUFR1} ./pm10bufr
#fi
 ln -s ${PREPBUFR4} ./ce318bufr
 ln -s ${PREPBUFR5} ./srfeabsbufr
 ln -s ${PREPBUFR6} ./srfescabufr
 ln -s ${PREPBUFR6} ./srfeextbufr
#
##################################################################################

echo " Copy fixed files and link CRTM coefficient files to working directory"

# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (regional only)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

### -------------------------------------- cwy step 3
if [ ${bk_core} = WRFCHEM_MOSAIC ] ; then
  BERROR=${BK_ROOT}/be.gcv
  BERROR_CHEM=${BK_ROOT}/be.gcv
  ANAVINFO=${BK_ROOT}/anavinfo_wrfmosaic
fi
### -------------------------------------- cwy step 3
if [ ${bk_core} = WRFCHEM_GOCART ] ; then
  BERROR=${BK_ROOT}/be.gocart/be.gcv
  BERROR_CHEM=${BK_ROOT}/be.gocart/be.gcv
  ANAVINFO=${BK_ROOT}/anavinfo_wrfgocart_cwy
fi
if [ ${bk_core} = WRFCHEM_PM25 ] ; then
#  BERROR=${FIX_ROOT}/${BYTE_ORDER}/wrf_chem_berror_big_endian
#  BERROR_CHEM=${FIX_ROOT}/${BYTE_ORDER}/wrf_chem_berror_big_endian
#  ANAVINFO=${FIX_ROOT}/anavinfo_wrfchem_pm25
  BERROR=${BK_ROOT}/be.gocart/be.gcv
  BERROR_CHEM=${BK_ROOT}/be.gocart/be.gcv
  ANAVINFO=${BK_ROOT}/anavinfo_wrfpm25_cwy
fi
if [ ${bk_core} = CMAQ ] ; then
  BERROR=${FIX_ROOT}/${BYTE_ORDER}/cmaq_berror_big_endian
  BERROR_CHEM=${FIX_ROOT}/${BYTE_ORDER}/cmaq_berror_big_endian
  ANAVINFO=${FIX_ROOT}/anavinfo_cmaq_pm25
fi

AEROINFO=${FIX_ROOT}/aeroinfo_aod.txt
OBERROR=${FIX_ROOT}/nam_errtable.r3dv
SATANGL=${FIX_ROOT}/global_satangbias.txt
SATINFO=${FIX_ROOT}/global_satinfo.txt
CONVINFO=${FIX_ROOT}/global_convinfo.txt
OZINFO=${FIX_ROOT}/global_ozinfo.txt
PCPINFO=${FIX_ROOT}/global_pcpinfo.txt

#  copy Fixed fields to working directory
 cp $ANAVINFO anavinfo
 cp $BERROR   berror_stats
 cp $BERROR_CHEM   berror_stats_chem
 cp $SATANGL  satbias_angle
 cp $SATINFO  satinfo
 cp $CONVINFO convinfo
 cp $OZINFO   ozinfo
 cp $PCPINFO  pcpinfo
 cp $OBERROR  errtable
 cp $AEROINFO aeroinfo
#
#    # CRTM Spectral and Transmittance coefficients
CRTM_ROOT_ORDER=${CRTM_ROOT}/${BYTE_ORDER}
emiscoef_IRwater=${CRTM_ROOT_ORDER}/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=${CRTM_ROOT_ORDER}/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=${CRTM_ROOT_ORDER}/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=${CRTM_ROOT_ORDER}/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=${CRTM_ROOT_ORDER}/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=${CRTM_ROOT_ORDER}/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=${CRTM_ROOT_ORDER}/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=${CRTM_ROOT_ORDER}/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=${CRTM_ROOT_ORDER}/FASTEM6.MWwater.EmisCoeff.bin
aercoef=${CRTM_ROOT_ORDER}/AerosolCoeff.bin
cldcoef=${CRTM_ROOT_ORDER}/CloudCoeff.bin

ln -s $emiscoef_IRwater ./Nalli.IRwater.EmisCoeff.bin
ln -s $emiscoef_IRice ./NPOESS.IRice.EmisCoeff.bin
ln -s $emiscoef_IRsnow ./NPOESS.IRsnow.EmisCoeff.bin
ln -s $emiscoef_IRland ./NPOESS.IRland.EmisCoeff.bin
ln -s $emiscoef_VISice ./NPOESS.VISice.EmisCoeff.bin
ln -s $emiscoef_VISland ./NPOESS.VISland.EmisCoeff.bin
ln -s $emiscoef_VISsnow ./NPOESS.VISsnow.EmisCoeff.bin
ln -s $emiscoef_VISwater ./NPOESS.VISwater.EmisCoeff.bin
ln -s $emiscoef_MWwater ./FASTEM6.MWwater.EmisCoeff.bin
ln -s $aercoef  ./AerosolCoeff.bin
ln -s $cldcoef  ./CloudCoeff.bin
# Copy CRTM coefficient files based on entries in satinfo file
for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
   ln -s ${CRTM_ROOT_ORDER}/${file}.SpcCoeff.bin ./
   ln -s ${CRTM_ROOT_ORDER}/${file}.TauCoeff.bin ./
done

for file in `awk '{if($1!~"!"){print $1}}' ./aeroinfo | sort | uniq` ;do
   ln -s ${CRTM_ROOT_ORDER}/${file}.SpcCoeff.bin ./
   ln -s ${CRTM_ROOT_ORDER}/${file}.TauCoeff.bin ./
done

# Only need this file for single obs test
 bufrtable=${FIX_ROOT}/prepobs_prep.bufrtable
 cp $bufrtable ./prepobs_prep.bufrtable

# for satellite bias correction
# Users may need to use their own satbias files for correct bias correction
cp ${GSI_ROOT}/fix/comgsi_satbias_in ./satbias_in
cp ${GSI_ROOT}/fix/comgsi_satbias_pc_in ./satbias_pc_in 

#
##################################################################################
# Set some parameters for use by the GSI executable and to build the namelist
echo " Build the namelist "

if [ ${bk_core} = WRFCHEM_MOSAIC ] ; then
 bk_core_arw='.true.'
 bk_if_netcdf='.true.'
 bk_core_cmaq='.false.'
 bk_wrf_pm2_5='.false.'
 bk_laeroana_gocart='.false.'
 bk_laeroana_mosaic='.true.'
fi
if [ ${bk_core} = WRFCHEM_GOCART ] ; then
 bk_core_arw='.true.'
 bk_if_netcdf='.true.'
 bk_core_cmaq='.false.'
 bk_wrf_pm2_5='.false.'
 bk_laeroana_gocart='.true.'
 bk_laeroana_mosaic='.false.'
fi
if [ ${bk_core} = WRFCHEM_PM25 ] ; then
 bk_core_arw='.true.'
 bk_if_netcdf='.true.'
 bk_core_cmaq='.false.'
 bk_wrf_pm2_5='.true.'
 bk_laeroana_gocart='.false.'
 bk_laeroana_mosaic='.false.'
fi
if [ ${bk_core} = CMAQ ] ; then
 bk_core_arw='.false.'
 bk_if_netcdf='.false.'
 bk_core_cmaq='.true.'
 bk_wrf_pm2_5='.false.'
 bk_laeroana_gocart='.false.'
 bk_laeroana_mosaic='.false.'
fi

# Build the GSI namelist on-the-fly
. $GSI_NAMELIST

#
###################################################
#  run  GSI
###################################################
echo ' Run GSI with' ${bk_core} 'background'

case $ARCH in
   'IBM_LSF')
      ${RUN_COMMAND} ./gsi.x < gsiparm.anl > stdout 2>&1  ;;

   * )
      ${RUN_COMMAND} ./gsi.x > stdout 2>&1  ;;
esac

##################################################################
#  run time error check
##################################################################
error=$?

if [ ${error} -ne 0 ]; then
  echo "ERROR: ${GSI} crashed  Exit status=${error}"
  exit ${error}
fi

#
##################################################################
#
#   GSI updating satbias_in
#
# GSI updating satbias_in (only for cycling assimilation)

# Copy the output to more understandable names
ln -s stdout      stdout.anl.${ANAL_TIME}
ln -s wrf_inout   wrfanl.${ANAL_TIME}
ln -s fort.201    fit_p1.${ANAL_TIME}
ln -s fort.202    fit_w1.${ANAL_TIME}
ln -s fort.203    fit_t1.${ANAL_TIME}
ln -s fort.204    fit_q1.${ANAL_TIME}
ln -s fort.207    fit_rad1.${ANAL_TIME}

# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
# NOTE:  Since we set miter=2 in GSI namelist SETUP, outer
#        loop 03 will contain innovations with respect to
#        the analysis.  Creation of o-a innovation files
#        is triggered by write_diag(3)=.true.  The setting
#        write_diag(1)=.true. turns on creation of o-g
#        innovation files.
#

loops="01 03"
for loop in $loops; do

case $loop in
  01) string=ges;;
  03) string=anl;;
   *) string=$loop;;
esac

#  Collect diagnostic files for obs types (groups) below
#   listall="conv amsua_metop-a mhs_metop-a hirs4_metop-a hirs2_n14 msu_n14 \
#          sndr_g08 sndr_g10 sndr_g12 sndr_g08_prep sndr_g10_prep sndr_g12_prep \
#          sndrd1_g08 sndrd2_g08 sndrd3_g08 sndrd4_g08 sndrd1_g10 sndrd2_g10 \
#          sndrd3_g10 sndrd4_g10 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 \
#          hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 \
#          amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua \
#          goes_img_g08 goes_img_g10 goes_img_g11 goes_img_g12 \
#          pcp_ssmi_dmsp pcp_tmi_trmm sbuv2_n16 sbuv2_n17 sbuv2_n18 \
#          omi_aura ssmi_f13 ssmi_f14 ssmi_f15 hirs4_n18 amsua_n18 mhs_n18 \
#          amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_las_f16 \
#          ssmis_uas_f16 ssmis_img_f16 ssmis_env_f16 mhs_metop_b \
#          hirs4_metop_b hirs4_n19 amusa_n19 mhs_n19"
 listall=`ls pe* | cut -f2 -d"." | awk '{print substr($0, 0, length($0)-3)}' | sort | uniq `

   for type in $listall; do
      count=`ls pe*${type}_${loop}* | wc -l`
      if [[ $count -gt 0 ]]; then
         cat pe*${type}_${loop}* > diag_${type}_${string}.${ANAL_TIME}
      fi
   done
done

#  Clean working directory to save only important files 
ls -l * > list_run_directory
if [[ ${if_clean} = clean  ]]; then
  echo ' Clean working directory after GSI run'
  rm -f *Coeff.bin     # all CRTM coefficient files
  rm -f pe0*           # diag files on each processor
  rm -f obs_input.*    # observation middle files
  rm -f siganl sigf0?  # background middle files
  rm -f fsize_*        # delete temperal file for bufr size
fi
#
#
exit 0
