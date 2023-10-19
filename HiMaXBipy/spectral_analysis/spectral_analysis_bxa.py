#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 17:23:18 2023

@author: David Kaltenbrunner
"""
from astropy.io import fits
from bxa.xspec.solver import BXASolver, XSilence
import bxa.xspec as bxa
import json
import numpy as np
import os


def lum(flux, distance):
    return 4. * np.pi * (distance * 1000 * 3.0857 * 10 ** 18) ** 2 * 10 ** flux


def plot_bxa(Plot, rebinning, src_files, ax_spec, ax_res, colors,
             src_markers, bkg_markers, model_linestyle, epoch_type,
             bkg_factors):
    spec_linestyle = ''
    bkg_linestyle = ''

    energies = {}
    fluxes = {}
    backgrounds = {}
    bkg_energies = {}
    models = {}
    bkg_models = {}
    residuals = {}
    Plot.setRebin(rebinning[0], rebinning[1])
    n_srcfiles = len(src_files)
    for igroup in range(n_srcfiles):
        if n_srcfiles == 1:
            label = f'{epoch_type}'
        else:
            label = f'{epoch_type} {igroup+1}'
        isource = 2*igroup + 1
        Plot.xAxis = 'keV'
        Plot('data')
        EkeV = Plot.x(isource)
        EkeV_err = Plot.xErr(isource)
        data = Plot.y(isource)
        data_err = Plot.yErr(isource)
        ax_spec.errorbar(EkeV, data, xerr=EkeV_err, yerr=data_err,
                         color=colors[igroup], marker=src_markers[igroup],
                         label=label, linestyle=spec_linestyle)
        bkg = (bkg_factors[igroup] * np.array(Plot.y(isource+1))).tolist()
        bkg_err = (bkg_factors[igroup] *
                   np.array(Plot.yErr(isource+1))).tolist()
        EkeV_bkg = Plot.x(isource+1)
        EkeV_bkg_err = Plot.xErr(isource+1)
        ax_spec.errorbar(EkeV_bkg, bkg, xerr=EkeV_bkg_err, yerr=bkg_err,
                         color=colors[igroup], marker=bkg_markers[igroup],
                         label=label, linestyle=spec_linestyle)

        model = Plot.model(isource)
        Esteps = []
        for ii, entry in enumerate(EkeV):
            Esteps.append(entry-EkeV_err[ii])
        Esteps.append(EkeV[-1]+EkeV_err[-1])
        ax_spec.stairs(model, Esteps, color=colors[igroup],
                       linestyle=model_linestyle)
        model_bkg = (bkg_factors[igroup] *
                     np.array(Plot.model(isource+1))).tolist()
        Esteps_bkg = []
        for ii, entry in enumerate(EkeV_bkg):
            Esteps_bkg.append(entry-EkeV_bkg_err[ii])
        Esteps_bkg.append(EkeV_bkg[-1]+EkeV_bkg_err[-1])
        ax_spec.stairs(model_bkg, Esteps_bkg, color=colors[igroup],
                       linestyle=model_linestyle)

        Plot('delchi')
        resid = Plot.y(isource)
        resid_err = Plot.yErr(isource)
        ax_res.errorbar(EkeV, resid, xerr=EkeV_err, yerr=resid_err,
                        color=colors[igroup], linestyle=bkg_linestyle)

        energies[label] = [EkeV, EkeV_err]
        fluxes[label] = [data, data_err]
        backgrounds[label] = [bkg, bkg_err]
        bkg_energies[label] = [EkeV_bkg, EkeV_bkg_err]
        models[label] = [model]
        residuals[label] = [resid, resid_err]
        bkg_models[label] = [model_bkg]

    output = {}
    output['EkeV'] = energies
    output['src_fluxes'] = fluxes
    output['bkg_fluxes'] = backgrounds
    output['EkeV_bkg'] = bkg_energies
    output['model_predictions'] = models
    output['residuals'] = residuals
    output['model_bkg'] = bkg_models

    return output


def fit_bxa(Xset, Fit, PlotManager, AllData, AllModels, Spectrum, Model,
            abund, distance, E_ranges, func, galnh, log, prompting, quantiles,
            src_files, statistic, suffix, resume, working_dir, Z):

    os.environ['WITHOUT_MULTINEST'] = '1'
    Xset.abund = abund
    Fit.statMethod = statistic
    Xset.allowPrompting = prompting
    PlotManager.splashPage = False  # reduces unnecessar output by xspec

    AllData.clear()
    AllModels.clear()

    backinfo = json.load(
        open(os.environ.get('EROBACK', "erosita_merged_1024.json")))
    jlo, jhi = backinfo['ilo'], backinfo['ihi']
    ilo, ihi = jlo, jhi
    clo = 1
    nchan = 1024

    # pca_models = [] #TODO: maybe for future release allow to let the user define pcamodel to use for bkg
    n_srcfiles = len(src_files)

    srcs = []
    bkgs = []
    bkg_files = []
    bkg_factors = []

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
#        src.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo, scaleType="lin",
#                     chanOffset=clo, chanWidth=1, sourceNum=1)
# the above stuff seems counterintuitive, trying out what happens if i leave it out
# pretty sure that it must not be set like that, tried in xspec. Not sure why it worked somehow still
# seems to not have any influence on the fit, no clue why.

        AllData(f"{2*ispec+2}:{2*ispec+2} " + bkgfile)
        bkg = AllData(2*ispec+2)
        bkg.ignore("1-%d" % ilo)
        bkg.ignore("%d-**" % (ihi))

        arf, rmf = bkg.multiresponse[0].arf, bkg.multiresponse[0].rmf
        for ii in range(1, n_srcfiles+1):
            bkg.multiresponse[ii] = rmf
            bkg.multiresponse[ii].arf = arf

        srcs.append(src)
        bkgs.append(bkg)
        bkg_files.append(bkgfile)
        bkg_factors.append(bkg_factor)

    for ispec in range(n_srcfiles):
        # make two responses, because starting with number 2 does not work.
        src = srcs[ispec]
        bkg = bkgs[ispec]
        bkgfile = bkg_files[ispec]
        for ii in range(2, n_srcfiles+2):
            src.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo,
                         scaleType="lin", chanOffset=clo, chanWidth=1.,
                         sourceNum=ii)
            bkg.dummyrsp(lowE=ilo, highE=ihi, nBins=ihi - ilo,
                         scaleType="lin", chanOffset=clo, chanWidth=1.,
                         sourceNum=ii)  # channels need to be set like that to work with the pca model, can't rescale to energies
        # delete the first response
        # bkg.multiresponse[0] = None

    transf_src, nHs_froz, nHs_mod, model_name = func(Model, AllModels, bxa,
                                                     galnh, Z, n_srcfiles)
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
                      modName=f'bkgmod{bkgid}').pcabkg.norm.values = \
                [1./bkg_factors[groupid], -1]
            # has to be set, see https://heasarc.gsfc.nasa.gov/docs/asca/abc_backscal.html
            if groupid+1 == bkgid:
                fac_src.values = [1, 0.1, 0.1, 0.1, 10, 10]
                fac_bkg.link = fac_src
                model = AllModels(groupNum=2*groupid+1,
                                  modName=f'bkgmod{bkgid}')
                norm = bxa.create_jeffreys_prior_for(model, fac_src)
                norm['name'] = f'log(bkg model factor {bkgid})'
                bkg_norms.append(norm)
            else:
                fac_src.values = [0, -1]
                fac_bkg.values = [0, -1]

    with XSilence():
        Fit.nIterations = 1000
        Fit.statMethod = 'cstat'

    transformations = transf_src + bkg_norms

    outfiles = f'{working_dir}/{model_name}{suffix}'
    if not os.path.exists(outfiles):
        os.mkdir(outfiles)

    log.info('running BXA')
    analyser = BXASolver(transformations=transformations,
                         outputfiles_basename=outfiles)
    analyser.run(resume=resume)

    AllData.show()
    AllModels.show()

    absorbed_F = []
    unabsorbed_L = []
    old_posterior = analyser.posterior.copy()
    for ispec in range(n_srcfiles):
        src = srcs[ispec]
        fluxes = []
        lums = []
        for band in E_ranges:
            old_nh = []
            flux = analyser.create_flux_chain(src, erange=f'{band[0]}'
                                              f' {band[1]}')
            fluxes_band = [np.percentile(flux, quantiles[0]),
                           np.percentile(flux, quantiles[1]),
                           np.percentile(flux, quantiles[2])]
            fluxes.append(fluxes_band)
            for nH in nHs_froz:
                old_nh.append(nH.values)
                nH.values = 1e-5
            for nH in nHs_mod:
                nH.values = 1e-5
            flux = analyser.create_flux_chain(src, erange=f'{band[0]}'
                                              f' {band[1]}')
            lums_band = [lum(np.percentile(flux, quantiles[0]), distance),
                         lum(np.percentile(flux, quantiles[1]), distance),
                         lum(np.percentile(flux, quantiles[2]), distance)]
            lums.append(lums_band)
            for inH, nH in enumerate(nHs_froz):
                nH.values = old_nh[inH]
            analyser.posterior[:, 0] = old_posterior[:, 0]
        absorbed_F.append(fluxes)
        unabsorbed_L.append(lums)

    AllData.show()
    AllModels.show()

    return absorbed_F, unabsorbed_L, bkg_factors
