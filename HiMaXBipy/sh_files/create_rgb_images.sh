#!/bin/bash

source @esass_location

ibin=80
nbinx=2700
nbiny=2700

RA0=@RA0
DE0=@Dec0

r1=@r1
r2=@r2
g1=@g1
g2=@g2
b1=@b1
b2=@b2

ener="${r1}-${r2}"
eneg="${g1}-${g2}"
eneb="${b1}-${b2}"
enet="${r1}-${b2}"

eventfile=@eventfile

smoothing=@smoothing
smooth="smooth${smoothing}"

#######################################
# reproject to new central coordinate #
#######################################
radec2xy file="${eventfile}".fits ra0="${RA0}" dec0="${DE0}"

#####################################
# apply pattern and flag selections #
#####################################
ofile=${eventfile}_cleaned.fits
infile=${eventfile}.fits
evtool $infile $ofile flag=0xC000F000 pattern=15 clobber=yes

###############################################
# create images for the individual telescopes #
###############################################
for itm in 1 2 3 4 5 6 7; do
    infile=$ofile
    file=${eventfile}_${itm}20.fits
    rm -f ima_*_${file}
    evtool $infile ima_${ener}_${file} rawxy_invert=no rawxy="2 383 2 383" emin=${r1} emax=${r2} image=yes telid=${itm} rebin=${ibin} size="${nbinx} ${nbiny}" &
    evtool $infile ima_${eneg}_${file} rawxy_invert=no rawxy="2 383 2 383" emin=${g1} emax=${g2} image=yes telid=${itm} rebin=${ibin} size="${nbinx} ${nbiny}" &
    evtool $infile ima_${eneb}_${file} rawxy_invert=no rawxy="2 383 2 383" emin=${b1} emax=${b2} image=yes telid=${itm} rebin=${ibin} size="${nbinx} ${nbiny}" &
    evtool $infile ima_${enet}_${file} rawxy_invert=no rawxy="2 383 2 383" emin=${r1} emax=${b2} image=yes telid=${itm} rebin=${ibin} size="${nbinx} ${nbiny}" &
    # wait until files are ready
    for ene in $ener $eneg $eneb ; do
        out=ima_${ene}_${file}
        echo ${out}
        while : ; do
            [[ -f "${out}" ]] && break
            echo "Pausing until file exists."
            sleep 60
        done
    done
done

####################################
# sum up images for all telescopes #
####################################
for ene in $ener $eneg $eneb $enet ; do
    itm=1
    file=${eventfile}_${itm}20.fits
    imfile=ima_${ene}_${file}
    cp ${imfile} tmp.fits
    for itm in 2 3 4 5 6 7 ; do
        file=${eventfile}_${itm}20.fits
        imfile=ima_${ene}_${file}
        farith tmp.fits ${imfile} tmp${itm}.fits ADD null=yes clobber=yes
        mv tmp${itm}.fits tmp.fits
    done
    file=${eventfile}_020.fits
    imfile=ima_${ene}_${file}
    mv tmp.fits ${imfile}
    date
    echo "${imfile} created"
done

########################
# create exposure maps #
########################
for itm in 1 2 3 4 5 6 7 ; do
    file=${eventfile}_${itm}20.fits
    rm -f exp_*_${file}
    expmap inputdatasets="ima_${enet}_${file}" \
        emin="${r1} ${g1} ${b1} ${r1}" emax="${r2} ${g2} ${b2} ${b2}" \
        templateimage="ima_${enet}_${file}" withdetmaps=yes \
        mergedmaps="exp_${ener}_${file} exp_${eneg}_${file} exp_${eneb}_${file} exp_${enet}_${file}" &
    sleep 120
done
for itm in 1 2 3 4 5 6 7 ; do
    file=${eventfile}_${itm}20.fits
    # wait until files are ready:
    for ene in $ener $eneg $eneb ; do
        out=exp_${ene}_${file}
        echo ${out}
        while : ; do
            [[ -f "${out}" ]] && break
            echo "Pausing until file exists."
            sleep 600
        done
    done
done

############################################
# sum up exposure maps from all telescopes #
############################################
for ene in $ener $eneg $eneb $enet ; do
    itm=1
    file=${eventfile}_${itm}20.fits
    imfile=exp_${ene}_${file}
    cp ${imfile} tmp.fits
    for itm in 2 3 4 5 6 7 ; do
        file=${eventfile}_${itm}20.fits
        imfile=exp_${ene}_${file}
        farith tmp.fits ${imfile} tmp${itm}.fits ADD null=yes clobber=yes
        mv tmp${itm}.fits tmp.fits
    done
    itm=0
    file=${eventfile}_${itm}20.fits
    imfile=exp_${ene}_${file}
    mv tmp.fits ${imfile}
    date
    echo "${imfile} created"
done

rm -f ima_*_${eventfile}_{1,2,3,4,5,6,7}20.fits exp_*_${eventfile}_{1,2,3,4,5,6,7}20.fits

#######################################
# divide images by their exposure map #
#######################################
itm=0
file=${eventfile}_${itm}20.fits
for ene in $ener $eneg $eneb $enet ; do
    farith ima_${ene}_${file} exp_${ene}_${file} rima_${ene}_${file} DIV null=yes clobber=yes
    echo "farith ima_${ene}_${file} exp_${ene}_${file} rima_${ene}_${file} DIV null=yes clobber=yes"
done

############################
# creating smoothed images #
############################

imat=ima_${enet}_${eventfile}_020.fits
outt=ima_${enet}_${eventfile}_020_${smooth}.fits

# adaptive smoothing:
rm -f ima_convolversset.fits ima_indeximageset.fits
asmooth inset=$imat outset=$outt nconvolvers=50 minwidth=0.0 maxwidth=42.6 widthliststyle=linear desiredsnr=${smoothing} \
        outconvolversset=ima_convolversset.fits outindeximageset=ima_indeximageset.fits writeconvolvers=yes smoothstyle=adaptive
sleep 60

for ene in $ener $eneg $eneb ; do
    ima=ima_${ene}_${eventfile}_020.fits
    out=ima_${ene}_${eventfile}_020_${smooth}.fits
    asmooth inset=$ima outset=$out inconvolversarray=ima_convolversset.fits inindeximagearray=ima_indeximageset.fits \
            writeconvolvers=no smoothstyle=withset &
done

# wait until files are ready:
for ene in $ener $eneg $eneb ; do
    out=ima_${ene}_${eventfile}_020_${smooth}.fits
    while : ; do
        [[ -f "${out}" ]] && break
        echo "Pausing until file ${out} exists."
        sleep 300
        [[ -f "STOP" ]] && exit
    done
done

# divide smoothed images by exposure maps:
itm=0
file=${eventfile}_${itm}20
for ene in $ener $eneg $eneb $enet ; do
    farith ima_${ene}_${file}_${smooth}.fits exp_${ene}_${file}.fits rima_${ene}_${file}_${smooth}.fits DIV null=yes clobber=yes
    echo "farith ima_${ene}_${file}_${smooth}.fits exp_${ene}_${file}.fits rima_${ene}_${file}_${smooth}.fits DIV null=yes clobber=yes"
done

exit
