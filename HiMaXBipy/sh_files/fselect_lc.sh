#!/bin/sh

source @esass_location

infile=@infile
pfile=@pfile
selection=@selection

fselect ${infile}[1] $pfile "${selection}" clobber=yes
