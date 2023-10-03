import os
from xspec import Spectrum, Model, Fit, Xset, PlotManager, Plot, AllModels, AllData
import json
from astropy.io import fits
import numpy as np
import bxa.xspec as bxa
from bxa.xspec.solver import BXASolver, XSilence
import logging

from HiMaXBipy.spectral_analysis.standard_models_bxa import apl

src_files = []
Z = 0.5
galnh = 0.06
os.environ['WITHOUT_MULTINEST'] = '1'
modelname = 'apl'
Xset.abund = 'wilm'
Fit.statMethod = 'cstat'
Xset.allowPrompting = False  # TODO: maybe make this an option
PlotManager.splashPage = False  # TODO what does this do?
func = apl
working_dir = ''
basename = ''
resume = False

AllData.clear()
AllModels.clear()

backinfo = json.load(
    open(os.environ.get('EROBACK', "erosita_merged_1024.json")))
jlo, jhi = backinfo['ilo'], backinfo['ihi']
ilo, ihi = jlo, jhi
clo = 1
nchan = 1024

pca_models = []
n_srcfiles = len(src_files)

srcs = []
bkgs = []
bkg_files = []

for ispec, specfile in enumerate(src_files):
    if ispec == 0:
        src = Spectrum(specfile)
    else:
        AllData(f'{2*ispec+1}:{2*ispec+1} ' + specfile)
        src = AllData(2*ispec+1)
    # added for proper definition of nchan
    src.ignore("**-**")
    _, nchan = map(int, src.ignoredString().split("-"))
    src.notice("**-**")
    # end of add
    src.ignore("**-0.2 8.0-**")
    for part in src.ignoredString().split():
        if "-" in part:
            klo, khi = map(int, part.split("-"))
        else:
            k = int(part)
            if k < 10:
                klo, khi = 1, k
            else:
                klo, khi = k, nchan

        if klo < 10:
            ilo = max(ilo, khi)
        else:
            ihi = min(ihi, klo)

    src.ignore("1-%d" % ilo)
    src.ignore("%d-%d" % (ihi, nchan))

    bkgfile = src.background.fileName

    src_header = fits.open(src.fileName)['SPECTRUM'].header
    bkg_header = fits.open(bkgfile)['SPECTRUM'].header
    bkg_factor = src_header['BACKSCAL'] / bkg_header['BACKSCAL']

    arf, rmf = src.multiresponse[0].arf, src.multiresponse[0].rmf
    src.background = None

    # copy response to source 2
    for ii in range(1, n_srcfiles+1):
        src.multiresponse[ii] = rmf
        src.multiresponse[ii].arf = arf

    # make unit response for background model:
    src.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo, scaleType="lin",
                 chanOffset=clo, chanWidth=1, sourceNum=1)

    AllData(f"{2*ispec+2}:{2*ispec+2} " + bkgfile)
    bkg = AllData(2*ispec+2)
    bkg.ignore("1-%d" % ilo)
    bkg.ignore("%d-**" % (ihi))

    srcs.append(src)
    bkgs.append(bkgs)
    bkg_files.append(bkgfile)

for ispec in range(n_srcfiles):
    # make two responses, because starting with number 2 does not work.
    src = srcs[ispec]
    bkg = bkgs[ispec]
    bkgfile = bkg_files[ispec]
    for ii in range(2, n_srcfiles+2):
        src.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo,
                     scaleType="lin", chanOffset=clo, chanWidth=1,
                     sourceNum=ii)
        bkg.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo,
                     scaleType="lin", chanOffset=clo, chanWidth=1,
                     sourceNum=ii)
    # delete the first response
    bkg.multiresponse[0] = None

transf_src = func(Model, AllModels, bxa, galnh, Z, n_srcfiles)
for ispec in range(n_srcfiles):
    bkgfile = bkg_files[ispec]
    Model('const*atable{%s_model.fits}' % (bkgfile),
          modName=f'bkgmod{ispec+1}', sourceNum=ispec+2)

bkg_norms = []
for groupid in range(n_srcfiles):
    for bkgid in range(1, n_srcfiles+1):
        fac_src = AllModels(groupNum=2*groupid+1,
                            modName=f'bkgmod{bkgid}').constant.factor
        fac_bkg = AllModels(groupNum=2*groupid+2,
                            modName=f'bkgmod{bkgid}').constant.factor
        AllModels(groupNum=2*groupid+1,
                  modName=f'bkgmod{bkgid}').pcabkg.norm.values = [1, -1]
        AllModels(groupNum=2*groupid+2,
                  modName=f'bkgmod{bkgid}').pcabkg.norm.values = [1, -1]  # does this have to be bkg_factor instead of 1?
        if groupid+1 == bkgid:
            fac_src.values = [1, 0.1, 0.1, 0.1, 10, 10]
            fac_bkg.link = fac_src
            model = AllModels(groupNum=2*groupid+1,
                              modName=f'bkgmod{bkgid}')
            norm = bxa.create_jeffreys_prior_for(model, fac_src)
            norm['name'] = f'back_factor{bkgid}'
            bkg_norms.append(norm)
        else:
            fac_src.values = [0, -1]
            fac_bkg.values = [0, -1]

transformations = transf_src + bkg_norms

outfiles = f'{working_dir}/{basename}'
if not os.path.exists(outfiles):
    os.mkdir(outfiles)

logging.info('running BXA')
analyser = BXASolver(transformations=transformations,
                     outputfiles_basename=outfiles)
analyser.run(resume=resume)
