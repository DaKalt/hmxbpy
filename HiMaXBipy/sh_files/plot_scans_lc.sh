
fracexp=0.15
fracexp_eRASS4=0.15
fractime=0.57
fracarea=0.05
infile=BeXRB2_newer_080156_eROSITA_PATall_1.0s020_LightCurve_00001.fits
pfile='BeXRB2_newer_080156_LC_TM020'
bin=42360
#binsize=30*1412s (spin period), approximately 3 erodays
mincounts=10

# complete light curve:
fselect ${infile}[1] ${pfile}fracexp${fracexp}_fullLC.fits "FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_fullLC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_UL.py ${pfile}fracexp${fracexp}_fullLC ${mincounts} # "Swift J053041.9-665426"

# eRASS0/1:
fselect ${infile}[1] ${pfile}fracexp${fracexp}_part1aLC.fits "TIME<6.4E08 && FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_part1aLC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp}_part1aLC ${mincounts} # "Swift J053041.9-665426"

# eRASS1/2:
fselect ${infile}[1] ${pfile}fracexp${fracexp}_part1bLC.fits "TIME>6.4E08 && TIME<6.5E08 && FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_part1bLC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp}_part1bLC ${mincounts} # "Swift J053041.9-665426"

# eRASS2/3:
fselect ${infile}[1] ${pfile}fracexp${fracexp}_part2LC.fits "TIME>6.5E08 && TIME<6.7E08 && FRACEXP>${fracexp}" clobber=yes
# use version which marks the time of the xmm observation
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_part2LC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp}_part2LC ${mincounts} # "Swift J053041.9-665426"

# eRASS3/4:
fselect ${infile}[1] ${pfile}fracexp${fracexp_eRASS4}_part3LC.fits "TIME>6.7E08 && TIME<6.8E08 && FRACEXP>${fracexp_eRASS4}" clobber=yes
# use version which marks the time of the xmm observation
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_UL.py ${pfile}fracexp${fracexp_eRASS4}_part3LC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp_eRASS4}_part3LC ${mincounts} # "Swift J053041.9-665426"

## eRASS3/4 tests:
#fselect ${infile}[1] ${pfile}fracexp${fracexp_eRASS4}_part3LC_fractime_fracarea.fits "TIME>6.7E08 && TIME<6.8E08 && FRACTIME>${fractime} && FRACAREA>${fracarea}" clobber=yes
## use version which marks the time of the xmm observation
##python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_UL.py ${pfile}fracexp${fracexp_eRASS4}_part3LC_fractime ${bin} # "Swift J053041.9-665426"
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fractime${fractime}_fracarea${fracarea}_part3LC ${mincounts} # "Swift J053041.9-665426"

#fselect ${infile}[1] ${pfile}fracexp${fracexp_eRASS4}_part3LC_fractime.fits "TIME>6.7E08 && TIME<6.8E08 && FRACTIME>${fractime}" clobber=yes
## use version which marks the time of the xmm observation
##python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_UL.py ${pfile}fracexp${fracexp_eRASS4}_part3LC_fractime ${bin} # "Swift J053041.9-665426"
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fractime${fractime}_part3LC ${mincounts} # "Swift J053041.9-665426"

#fselect ${infile}[1] ${pfile}fracexp${fracexp_eRASS4}_part3LC_fracarea.fits "TIME>6.7E08 && TIME<6.8E08 && FRACAREA>${fracarea}" clobber=yes
## use version which marks the time of the xmm observation
##python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_UL.py ${pfile}fracexp${fracexp_eRASS4}_part3LC_fractime ${bin} # "Swift J053041.9-665426"
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracarea${fracarea}_part3LC ${mincounts} # "Swift J053041.9-665426"

# eRASS4/5:
fselect ${infile}[1] ${pfile}fracexp${fracexp}_part4LC.fits "TIME>6.8E08 && FRACEXP>${fracexp}" clobber=yes
# use version which marks the time of the xmm observation
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_part4LC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp}_part4LC ${mincounts} # "Swift J053041.9-665426"

# checking lightcurve for all TMs
for j in 1 2 3 4 6
do
infile=BeXRB2_newer_080156_eROSITA_PATall_1.0s${j}20_LightCurve_00001.fits
pfile="BeXRB2_newer_080156_LC_TM${j}20"
fselect ${infile}[1] ${pfile}fracexp${fracexp}_fullLC.fits "FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_timebin_UL.py ${pfile}fracexp${fracexp}_fullLC ${bin} # "Swift J053041.9-665426"
python ~/galaxy/eROSITA/LMC/eRASS/sources/080156/src_BeXRB2_newer/plot_lc_BeXRB2_mincounts_UL.py ${pfile}fracexp${fracexp}_fullLC ${mincounts} # "Swift J053041.9-665426"
done

#
exit
#

#infile=BeXRB2_newer_080156_eROSITA_PATall_TMon_1.0s020_LightCurve_00001.fits
#pfile='BeXRB2_newer_080156_LC_TMon'

#fselect ${infile}[1] ${pfile}fracexp${fracexp}_fullLC.fits "FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/plot_scans_lc.py ${pfile}fracexp${fracexp}_fullLC

#fselect ${infile}[1] ${pfile}fracexp${fracexp}.fits "TIME<6.4E08 && FRACEXP>${fracexp}" clobber=yes
#python ~/galaxy/eROSITA/plot_scans_lc.py ${pfile}fracexp${fracexp}
