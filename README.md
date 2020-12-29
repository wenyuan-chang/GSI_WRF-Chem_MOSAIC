# GSI_WRF-Chem_MOSAIC
A revised 3D-Var GSI DA code for the WRF-Chem MOSAIC aerosols

## 1. Description
This GSI code is revised based on the official GSI v3.7. It is built and run in the same way as the official version. The user is recommended to read the official user's guide first and  know how to run the official GSI version, create the observational bufr file, and make the background error file with GEN_BE.

This revised GSI version works for the WRF-Chem MOSAIC aerosols with four or eight size bins. It can assimilate the 550nm MODIS AOD, CE318 AOD at four wavelengths (440, 675, 870, 1020 nm), aerosol scattering coefficient at three wavelengths (450, 525, 635 nm) of nephelometer, aerosol absorption coefficient of aethalometer at three wavelengths (470, 520, 660 nm), and the surface concentrations of PM2.5 and PM10. For briefly, this page use [GSI] to denote the main directory of comGSIv3.7_EnKFv1.3_cwy1.1_20201123

User can check the revised codes by
```
cd [GSI]
grep 'cwy mosaic' *.f90
```

## 2. How to build

```
cd [GSI]
mkdir build
cd build
csh
make -j2
```

## 3. How to run

```
cd run01
cp ../[GSI]/build/bin/gsi.x .

vi comgsi_run_chem.ksh
# Edit this file. 
# Make sure that the case-specified absolute directory paths are correct.
# For example, JOB_DIR, OBS_ROOT, BK_ROOT, GSI_ROOT, CRTM_ROOT ...

vi comgsi_run_chem.ksh
# Edit this file. 
# In the section of OBS_INPUT::, only keep the observational bufr files which would be assimilated.

cd INPUT
```
There are three anavinfo_wrfmosaic_chem_xxx files.
<ol>
  <li>anavinfo_wrfmosaic_chem_aer_4bin is for the MOSAIC 4bin aerosols. Aerosol chemical compoistions per size bin is control variable.</li>
  <li>anavinfo_wrfmosaic_chem_mas_4bin is for the MOSAIC 4bin aerosols. Total aerosol mass per size bin is control variable.</li>
  <li>anavinfo_wrfmosaic_chem_mas_8bin is for the MOSAIC 8bin aerosols. Total aerosol mass per size bin is control variable.</li>
</ol>
You do not need to edit any of these anavinfo files. Just link the anavinfo you want to anavinfo_wrfmosaic.

For example

```
ln -sf anavinfo_wrfmosaic_chem_aer_4bin anavinfo_wrfmosaic
```
  
Before running the code, make sure that
<ol>
  <li>The link files are fine, which point to the observations (in bufr format), the WRF-Chem model result, and the background error files (created by GEN_BE). It's the user's responsibility to prepare these data files.</li>
  <li>The background error matchs the control variables in the anavinfor_wrfmosaic that you chose.</li>
  <li>All the links are correctly specified in comgsi_run_chem.ksh.</li>
  <li>CRTM_v2.3.0 is in the correct place where GSI can read it normally. CRTM_v2.3.0 contains a few big files. User can find the CRTM data in the official GSI code package.</li>
</ol>
  
```
qsub comgsi_run_chem.ksh
# run GSI

tail -f example/stdout
#   check the data assimilation. 
```
The analyzed file is wrf_inout

## Acknowledgmnts
This work is supported by the National Key Research and Development Program of China (Grant number 2016YFE0201400).

## References
<ol>
  <li>GSI User’s Guide: https://dtcenter.ucar.edu/com-GSI/users/docs/</li>
  <li>Official GSI download page: https://dtcenter.ucar.edu/com-GSI/users/downloads/index.php </li>
  <li>Chang, W., Zhang, Y., Li, Z., Chen, J., and Li, K.: Improving the Sectional MOSAIC Aerosol models of WRF-Chem with the revised Gridpoint Statistical Interpolation System and multi-wavelength aerosol optical measurements: DAO-K experiment 2019 at Kashi, near the Taklamakan Desert, northwestern China, Atmos. Chem. Phys. Discuss. [preprint], https://doi.org/10.5194/acp-2020-825, in review, 2020.</li>
</ol>

## Feedback
Dr. Wenyuan Chang

changwy@mail.iap.ac.cn

State Key Laboratory of Atmospheric Boundary Layer Physics and Atmospheric Chemistry (LAPC)

Institute of Atmospheric Physics, Chinese Academy of Sciences

Beijing 100029, China.
