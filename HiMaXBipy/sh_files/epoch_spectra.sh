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

region=080156

#epoch1
infiles="${location}/em00_${region}_020_EventList_c020.fits ${location}/em01_${region}_020_EventList_c020.fits"
epoch="1"
tmin="6.2E08"
tmax="6.4E08"

evtool eventfiles="${infiles}" outfile="${location}/epoch${epoch}_${region}_020_EventList.fits" gti="${tmin} ${tmax}" clobber=yes
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done


#epoch2
infiles="${location}/em01_${region}_020_EventList_c020.fits ${location}/em02_${region}_020_EventList_c020.fits"
epoch="2"
tmin="6.4E08"
tmax="6.5E08"

evtool eventfiles="${infiles}" outfile="${location}/epoch${epoch}_${region}_020_EventList.fits" gti="${tmin} ${tmax}" clobber=yes
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done


#epoch3
infiles="${location}/em02_${region}_020_EventList_c020.fits ${location}/em03_${region}_020_EventList_c020.fits"
epoch="3"
tmin="6.5E08"
tmax="6.7E08"

evtool eventfiles="${infiles}" outfile="${location}/epoch${epoch}_${region}_020_EventList.fits" gti="${tmin} ${tmax}" clobber=yes
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done

#epoch4
infiles="${location}/em03_${region}_020_EventList_c020.fits ${location}/em04_${region}_020_EventList_c020.fits"
epoch="4"
tmin="6.7E08"
tmax="6.8E08"

evtool eventfiles="${infiles}" outfile="${location}/epoch${epoch}_${region}_020_EventList.fits" gti="${tmin} ${tmax}" clobber=yes
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done


#epoch5
infiles="${location}/em04_${region}_020_EventList_c020.fits ${location}/em05_${region}_020_EventList_c947.fits"
epoch="5"
tmin="6.8E08"
tmax="7.0E08"

evtool eventfiles="${infiles}" outfile="${location}/epoch${epoch}_${region}_020_EventList.fits" gti="${tmin} ${tmax}" clobber=yes
# spectra on chip detectors for each eRASS fun
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors for each eRASS run
srctool eventfiles="${location}/epoch${epoch}_${region}_020_EventList.fits"\
        prefix="${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${product_dir}/src.reg" backreg="${product_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for j in 0 1 2 3 4 6; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMon${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done
for j in 0 5 7; do
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g1.fits comm="group min 1 & exit" clobber=yes
  grppha ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoff${j}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_epoch${epoch}_eROSITA_PAT${PAT}_TMoffS${j}20_SourceSpec_00001_g20.fits comm="group min 20 & exit" clobber=yes
done

mv ${location}/epoch* ${product_dir}/*
