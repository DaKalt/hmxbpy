#!/bin/bash

source @esass_location

ra=@right_ascension
dec=@declination

main_dir=@main_name
product_dir=@result_dir

sources='@sources_list'

epoch=@epoch

old_wd=$PWD

cd $product_dir

# spectra on chip detectors
srctool eventfiles="${sources}" \
        prefix="${product_dir}/${epoch}_bxa_bkg_" srccoord="fk5;$ra,$dec" srcreg="${main_dir}/bkg_snr.reg" backreg="${main_dir}/src.reg" exttype="POINT" \
        todo="SPEC ARF RMF EVENTS" insts="1 2 3 4 6" tstep=0.05 xgrid="0.5 1.0" gtitype="GTI" psftype="2D_PSF" clobber=yes \

cd $old_wd
