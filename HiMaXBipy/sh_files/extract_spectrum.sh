source /home/erosita/sw/eSASSusers_201009/bin/esass-init.sh

ra=05h29m14.942s
dec=-66d24m46.030s

ra1=05h28m27.120s
dec1=-66d27m36.000s

#change this part with ds9 image
sizeS1=35

sizeBG1=55
#sizeBG2=120

PAT=all

do_lc=true
do_lc_plotting=false
do_spec=true
do_spec_plotting=false

src=BeXRB2_newer
location="$HOME/galaxy/eROSITA/LMC/eRASS/sources/080156"

product_dir="$HOME/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer"

echo "fk5;circle($ra,$dec,$sizeS1\")" > ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/src.reg
echo "fk5;circle($ra1,$dec1,$sizeBG1\")" > ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/bkg.reg

#lim_fracexp=0.15

region=080156
sources="${location}/em00_${region}_020_EventList_c020.fits "
for i in 1 2 3 4
do
	sources+="${location}/em0${i}_${region}_020_EventList_c020.fits "
done
sources+="${location}/em05_${region}_020_EventList_c947.fits"
sources="${sources%%*( )}"
#echo "${product_dir}/src.reg"

if $do_lc
then
#full light curve, 1s bin, TMall
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_1.0s" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 5 6 7" \
        srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="1.0" lcemin="0.2" lcemax="8" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes"

#full light curve, 1s bin, TM12346
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_1.0s" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 6" \
        srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="1.0" lcemin="0.2" lcemax="8" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes"

#both of the above for each eRASS
for i in 1 2 3 4
do
srctool eventfiles="${location}/em0${i}_${region}_020_EventList_c020.fits" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_1.0s" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 5 6 7" \
        srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="1.0" lcemin="0.2" lcemax="8" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes"

#full light curve, 1s bin, TM12346
srctool eventfiles="${location}/em0${i}_${region}_020_EventList_c020.fits" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon_1.0s" suffix=".fits" todo="LC LCCORR" insts="1 2 3 4 6" \
        srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="1.0" lcemin="0.2" lcemax="8" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes"
done

for j in 1 2 3 4 6
do
srctool eventfiles="${sources}" \
        srccoord="J2000;$ra,$dec" prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_1.0s_${j}20" suffix=".fits" todo="LC LCCORR" insts="${j}" \
        srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" lctype="REGULAR+" lcpars="1.0" lcemin="0.2" lcemax="8" tstep=0.05 xgrid="0.5 1.0" \
        gtitype="GTI" psftype="2D_PSF" clobber="yes"
done
fi



if $do_lc_plotting
then
#echo "plot LC here"
#infile=${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_1.0s020_LightCurve_00001.fits
#pfile=${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_1.0s020_LightCurve_00001
#fselect ${infile}[1] ${pfile}_fracexp${lim_fracexp}.fits "FRACEXP>${lim_fracexp}" clobber=yes
#python ${product_dir}/plot_lc*.py ${pfile}_fracexp${lim_fracexp}

#infile=${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_1.0s020_LightCurve_00001.fits 
#pfile=${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon_1.0s020_LightCurve_00001
#fselect ${infile}[1] ${pfile}_fracexp${lim_fracexp}.fits "FRACEXP>${lim_fracexp}" clobber=yes
#python ${product_dir}/plot_lc*.py ${pfile}_fracexp${lim_fracexp}

#for i in 1 2 3 4
#do
#infile=${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_1.0s020_LightCurve_00001.fits 
#pfile=${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_1.0s020_LightCurve_00001
#fselect ${infile}[1] ${pfile}_fracexp${lim_fracexp}.fits "FRACEXP>${lim_fracexp}" clobber=yes
#python ${product_dir}/plot_lc*.py ${pfile}_fracexp${lim_fracexp}

#infile=${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon_1.0s020_LightCurve_00001.fits 
#pfile=${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon_1.0s020_LightCurve_00001
#fselect ${infile}[1] ${pfile}_fracexp${lim_fracexp}.fits "FRACEXP>${lim_fracexp}" clobber=yes
#python ${product_dir}/plot_lc*.py ${pfile}_fracexp${lim_fracexp}
cd ${product_dir}
. ${product_dir}/plot_scans_lc.sh

#done

fi

if $do_spec
then

# spectra on chip detectors
srctool eventfiles="${sources}" \
        prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors
srctool eventfiles="${sources}" \
        prefix="${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for i in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon${i}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMon${i}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for i in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoff${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoffS${i}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoff${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_eROSITA_PAT${PAT}_TMoffS${i}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done


#eRASS0:1
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/em00_${region}_020_EventList_c020.fits ${location}/em01_${region}_020_EventList_c020.fits"\
        prefix="${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/em00_${region}_020_EventList_c020.fits ${location}/em01_${region}_020_EventList_c020.fits"\
        prefix="${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em00_1_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done


#eRASS2:4
for i in 2 3 4
do

# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/em0${i}_${region}_020_EventList_c020.fits"\
        prefix="${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/em0${i}_${region}_020_EventList_c020.fits"\
        prefix="${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_em0${i}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done

done

source ${location}/epoch_spectra.sh

fi


if $do_spec_plotting
then

echo "plot spectra here"

fi
