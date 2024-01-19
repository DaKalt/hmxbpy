#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 17:23:18 2023

@author: David Kaltenbrunner
"""
from astropy.io import fits
from bxa.xspec.solver import BXASolver, XSilence
import bxa.xspec as bxa
import corner
import json
import logging
from matplotlib.dates import EPOCH_OFFSET
from matplotlib.figure import Figure
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import numpy as np
import os
from xspec import Xset, Fit, PlotManager, AllData, AllModels, Spectrum, Model,\
    Plot

from HiMaXBipy.spectral_analysis.modded_function import PredictionBand,\
    posterior_predictions_plot


def lum(flux, distance):
    return 4. * np.pi * (distance * 1000 * 3.0857 * 10 ** 18) ** 2 * 10 ** flux


def plot_bxa(rebinning, src_files, ax_spec, ax_res, colors,
             src_markers, bkg_markers, epoch_type, bkg_factors, analyser,
             src_linestyles, bkg_linestyle, hatches, ntransf):
    energies = {}
    fluxes = {}
    backgrounds = {}
    bkg_energies = {}
    models = {}
    bkg_models = {}
    residuals = {}
    Plot.device = '/null'
    Plot.setRebin(rebinning[0], rebinning[1])
    n_srcfiles = len(src_files)
    for igroup in range(n_srcfiles):
        if n_srcfiles > 1:
            label = f'{epoch_type} {igroup+1}'
        else:
            label = f'{epoch_type}'
        isource = 2*igroup + 1
        Plot.xAxis = 'keV'
        Plot('data')
        EkeV = Plot.x(isource).copy()
        EkeV_err = Plot.xErr(isource).copy()
        data = Plot.y(isource).copy()
        data_err = Plot.yErr(isource).copy()
        ax_spec.errorbar(EkeV, data, xerr=EkeV_err, yerr=data_err,
                         color=colors[igroup], marker=src_markers[igroup],
                         label=label, linestyle='', zorder=5)
        bkg = (bkg_factors[igroup] *
               np.array(Plot.y(isource+1).copy())).tolist()
        bkg_err = (bkg_factors[igroup] *
                   np.array(Plot.yErr(isource+1).copy())).tolist()
        EkeV_bkg = Plot.x(isource+1).copy()
        EkeV_bkg_err = Plot.xErr(isource+1).copy()
        ax_spec.errorbar(EkeV_bkg, bkg, xerr=EkeV_bkg_err, yerr=bkg_err,
                         color=colors[igroup], marker=bkg_markers[igroup],
                         linestyle='', zorder=4)

        model = Plot.model(isource).copy()
        Esteps = []
        for ii, entry in enumerate(EkeV):
            Esteps.append(entry-EkeV_err[ii])
        Esteps.append(EkeV[-1]+EkeV_err[-1])
        # ax_spec.stairs(model, Esteps, color='black',
        #                linestyle=model_linestyle)
        model_bkg = (bkg_factors[igroup] *
                     np.array(Plot.model(isource+1).copy())).tolist()
        Esteps_bkg = []
        for ii, entry in enumerate(EkeV_bkg):
            Esteps_bkg.append(entry-EkeV_bkg_err[ii])
        Esteps_bkg.append(EkeV_bkg[-1]+EkeV_bkg_err[-1])
        # ax_spec.stairs(model_bkg, Esteps_bkg, color='red',
        #                linestyle=model_linestyle)

        Plot('delchi')
        resid = Plot.y(isource).copy()
        resid_err = Plot.yErr(isource).copy()
        ax_res.errorbar(EkeV, resid, xerr=EkeV_err, yerr=resid_err,
                        color=colors[igroup], linestyle='',
                        marker=src_markers[igroup])

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

    for igroup in range(n_srcfiles):
        models = []
        bands = []
        Plot.add = False
        Plot.background = False
        for content in posterior_predictions_plot(analyser, plottype='data',
                                                  nsamples=100, group=2*igroup+1):
            xmid = content[:, 0]
            ndata_columns = 6 if Plot.background else 4
            ncomponents = content.shape[1]-ndata_columns
            model_contributions = []
            for component in range(ncomponents):
                y = content[:, ndata_columns+component]
                if component >= len(bands):
                    bands.append(PredictionBand(xmid, ax_spec))
                bands[component].add(y)
                model_contributions.append(y)
            models.append(model_contributions)

        for band, label in zip(bands, 'model prediction'):
            if label == 'ignore':
                continue
            lineargs = dict(drawstyle='steps', color=colors[igroup], zorder=3,
                            linewidth=0.8, linestyle=src_linestyles[igroup])
            shadeargs = dict(color=lineargs['color'], zorder=2)
            band.shade(alpha=0.5, **shadeargs)
            shadeargs = dict(**shadeargs, hatch=hatches[igroup])
            band.shade(q=0.9973/2, alpha=0.2, **shadeargs)
            band.line(label=label, **lineargs)

        models = []
        bands = []
        Plot.add = False
        Plot.background = False
        for content in posterior_predictions_plot(analyser, plottype='data',
                                                  nsamples=100, group=2*igroup+2):
            xmid = content[:, 0]
            ndata_columns = 6 if Plot.background else 4
            ncomponents = content.shape[1]-ndata_columns
            model_contributions = []
            for component in range(ncomponents):
                y = bkg_factors[0] * content[:, ndata_columns+component]
                if component >= len(bands):
                    bands.append(PredictionBand(xmid, ax_spec))
                bands[component].add(y)
                model_contributions.append(y)
            models.append(model_contributions)

        for band, label in zip(bands, 'model prediction'):
            if label == 'ignore':
                continue
            lineargs = dict(drawstyle='steps', color=colors[igroup], zorder=1,
                            linestyle=bkg_linestyle, linewidth=0.8)
            shadeargs = dict(color=lineargs['color'], zorder=0)
            band.shade(alpha=0.5, **shadeargs)
            shadeargs = dict(**shadeargs, hatch=hatches[igroup])
            band.shade(q=0.9973/2., alpha=0.2, **shadeargs)
            band.line(label=label, **lineargs)

    return output


def fit_bxa(abund, distance, E_ranges, func, galnh, log, prompting, quantiles,
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
                norm['name'] = f'log(bkg mod f{bkgid})'
                bkg_norms.append(norm)
            else:
                fac_src.values = [0, -1]
                fac_bkg.values = [0, -1]

    with XSilence():
        Fit.nIterations = 1000
        Fit.statMethod = 'cstat'

    transformations = transf_src + bkg_norms
    ntransf = len(transf_src) #this is for corner plot without background norms

    outfiles = f'{working_dir}/{model_name}{suffix}'
    if not os.path.exists(outfiles):
        os.mkdir(outfiles)

    log.info('running BXA')
    analyser = BXASolver(transformations=transformations,
                         outputfiles_basename=outfiles)
    results = analyser.run(resume=resume)

    AllData.show()
    AllModels.show()

    absorbed_F = []
    unabsorbed_L = []
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
            if band == E_ranges[0] and ispec == 0:
                AllModels.show()
            for inH, nH in enumerate(nHs_froz):
                nH.values = old_nh[inH]
            analyser.set_best_fit()
        absorbed_F.append(fluxes)
        unabsorbed_L.append(lums)

    #this is just to check if analyser.set_best_fit() does its job
    AllData.show()
    AllModels.show()

    return absorbed_F, unabsorbed_L, bkg_factors, analyser, ntransf

def plot_corner(analyser, ntransf, log) -> Figure:
    """Make a corner plot with corner. Altered from ultranest."""
    paramnames = analyser.results['paramnames'][:ntransf]
    data = np.array(analyser.results['weighted_samples']['points']
                    .T[:ntransf].T)
    weights = np.array(analyser.results['weighted_samples']['weights'])
    cumsumweights = np.cumsum(weights)

    mask = cumsumweights > 1e-4

    if mask.sum() == 1:
        if log is not None:
            warn = 'Posterior is still concentrated in a single point:'
            for i, p in enumerate(paramnames):
                v = analyser.results['samples'][mask,i]
                warn += "\n" + '    %-20s: %s' % (p, v)

            log.warning(warn)
            log.info('Try running longer.')
        return plt.figure()

    # monkey patch to disable a useless warning;
    # this is for now turned off to check
    # oldfunc = logging.warning
    # logging.warning = lambda *args, **kwargs: None
    fig = corner.corner(data[mask,:], weights=weights[mask],
                        labels=paramnames, show_titles=True, quiet=True)
    # logging.warning = oldfunc
    return fig

def write_tex(tex_file, tex_info, abs_F, unabs_L, analyser, quantiles):
    tex_file.write('\\begin{{tabular}}{{%s}}\n' % ('c' * (1+len(tex_info)+2)))
    tex_file.write('\\hline\\hline\n')
    line_1 = ''
    for i in range(len(tex_info)):
        line_1 += tex_info[i][0] + ' & '
    line_2 = ''
    for i in range(len(tex_info)):
        line_2 += tex_info[i][1] + ' & '
    tex_file.write(f'Part & {line_1}'
                                '\\mbox{{F$_{{\\rm x}}$}} '
                                '& \\mbox{{L$_{{\\rm x}}$}} \\\\ \n')
    tex_file.write(f'-- & {line_2}'
                                '$\\times$erg cm$^{{-2}}$s$^{{-1}}$ '
                                '& $\\times$erg s$^{{-1}}$ \\\\ \n')
    tex_file.write('\\hline\n')
    for i in range(len(abs_F)):
        tex_file.write('%s\\\\ \n' % ('& '*(1 + len(tex_info) + 1)))
        line = f'{i+1} & '
        for j in range(len(tex_info)):
            if len(tex_info[j][2]) == 1:
                if i > 0:
                    line += ' & '
                    continue
            index = tex_info[j][2][i]
            data = analyser.results['weighted_samples']['points'].T[index]
            weights = np.array(analyser.results['weighted_samples']['weights'])
            cumsumweights = np.cumsum(weights)

            mask = cumsumweights > 1e-4

            lower = corner.quantile(data[mask], quantiles[0]/100.,
                                    weights[mask])[0]
            median = corner.quantile(data[mask], quantiles[1]/100.,
                                     weights[mask])[0]
            upper = corner.quantile(data[mask], quantiles[2]/100.,
                                    weights[mask])[0]
            if tex_info[j][3] == 'log':
                lower = 10 ** lower
                median = 10 ** median
                upper = 10 ** upper 
            line += f'{median}$^{{+{upper-median}}}_{{-{median-lower}}}$ & '
        lower = abs_F[i][0]
        median = abs_F[i][1]
        upper = abs_F[i][2]
        line += f'{median}$^{{+{upper-median}}}_{{-{median-lower}}}$ & '
        lower = unabs_L[i][0]
        median = unabs_L[i][1]
        upper = unabs_L[i][2]
        line += f'{median}$^{{+{upper-median}}}_{{-{median-lower}}}$ \\\\ \n'
        tex_file.write(line)
    tex_file.write('%s\\\\ \n' % ('& '*(1 + len(tex_info) + 1)))