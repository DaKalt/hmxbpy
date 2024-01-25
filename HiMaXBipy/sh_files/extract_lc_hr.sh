#!/bin/bash

source @esass_location

PAT=all

ra=@right_ascension
dec=@declination

src=@source_name

main_dir=@main_name
product_dir=@result_dir

region=@region_code

sources='@sources_list'

emin0=@emin0
emax0=@emax0
emin1=@emin1
emax1=@emax1
emin2=@emin2
emax2=@emax2

binning=@binning

#full light curve, 1s bin, TMall
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_${binning}s_HR_" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 5 6 7" \
        srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" lctype="REGULAR-" lcpars="${binning}" lcemin="${emin0} ${emin1} ${emin2}" lcemax="${emax0} ${emax1} ${emax2}" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes" writeinsts="0"

#full light curve, 1s bin, TM12346
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_${binning}s_HR_" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 6" \
        srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" lctype="REGULAR-" lcpars="${binning}" lcemin="${emin0} ${emin1} ${emin2}" lcemax="${emax0} ${emax1} ${emax2}" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes" writeinsts="0"

#full light curve, 1s bin, TM57
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoff_${binning}s_HR_" suffix=".fits" todo="LC LCCORR" insts="5 7" \
        srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" lctype="REGULAR-" lcpars="${binning}" lcemin="${emin0} ${emin1} ${emin2}" lcemax="${emax0} ${emax1} ${emax2}" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes" writeinsts="0"


#all of this is done above anyways
#for j in 1 2 3 4 6
#do
#srctool eventfiles="${sources}" \
#        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_${binning}s_${emin}keV_${emax}keV_${j}20_" suffix=".fits" todo="LC LCCORR" insts="${j}" \
#        srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="${binning}" lcemin="${emin}" lcemax="${emax}" tstep=0.05 xgrid="0.5 1.0" \
#        gtitype="GTI" psftype="2D_PSF" clobber="yes"
#done
#
#for j in 5 7
#do
#srctool eventfiles="${sources}" \
#        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoff_${binning}s_${emin}keV_${emax}keV_${j}20_" suffix=".fits" todo="LC LCCORR" insts="${j}" \
#        srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="${binning}" lcemin="${emin}" lcemax="${emax}" tstep=0.05 xgrid="0.5 1.0" \
#        gtitype="GTI" psftype="2D_PSF" clobber="yes"
#done
