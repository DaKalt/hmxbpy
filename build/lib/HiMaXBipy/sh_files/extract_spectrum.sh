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

group=@group

mode=@mode

epoch=@epoch

cd $product_dir

# spectra on chip detectors
srctool eventfiles="${sources}" \
        prefix="${product_dir}/${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMon" srccoord="fk5;$ra,$dec" srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

# spectra off chip detectors
srctool eventfiles="${sources}" \
        prefix="${product_dir}/${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMoff" srccoord="fk5;$ra,$dec" srcreg="${main_dir}/src.reg" backreg="${main_dir}/bkg.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="5 7" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

for i in 0 1 2 3 4 6; do
  grppha ./${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMon${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMon${i}20_SourceSpec_00001_g${group}.fits comm="group min ${group} & exit" clobber=yes
done
for i in 0 5 7; do
  grppha ./${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMoff${i}20_SourceSpec_00001.fits ${product_dir}/${src}_${region}_${epoch}_eROSITA_${mode}_PAT${PAT}_TMoff${i}20_SourceSpec_00001_g${group}.fits comm="group min ${group} & exit" clobber=yes
done
