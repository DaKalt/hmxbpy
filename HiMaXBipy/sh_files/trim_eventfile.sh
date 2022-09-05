#!/bin/sh

source @esass_location

infiles=@infiles
outfile=@outfile
tmin=@start
tmax=@stop

evtool eventfiles="${infiles}" outfile="${outfile}" gti="${tmin} ${tmax}" clobber=yes
