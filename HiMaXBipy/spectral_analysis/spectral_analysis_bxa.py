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
from matplotlib import gridspec
from matplotlib.dates import EPOCH_OFFSET
from matplotlib.figure import Figure
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import numpy as np
import os
from xspec import Xset, Fit, PlotManager, AllData, AllModels, Spectrum, Model,\
    Plot

from HiMaXBipy.io.package_data import num2text
from HiMaXBipy.spectral_analysis.modded_function import PredictionBand,\
    posterior_predictions_plot


def lum(flux, distance):
    return 4. * np.pi * (distance * 1000 * 3.0857 * 10 ** 18) ** 2 * 10 ** flux

def plot_bxa(rebinning, src_files, ax_spec, ax_res, colors,
             src_markers, bkg_markers, epoch_type, bkg_factors, analyser,
             src_linestyles, bkg_linestyle, hatches, ntransf,
             plot_src=False, ):
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
            color_srcbkg = colors[igroup]
            color_bkg = colors[igroup]
        else:
            label = f'{epoch_type}'
            color_srcbkg = colors[0]
            color_bkg = colors[1]
        isource = 2*igroup + 1
        Plot.xAxis = 'keV'
        Plot('data')
        EkeV = Plot.x(isource).copy()
        EkeV_err = Plot.xErr(isource).copy()
        data = Plot.y(isource).copy()
        data_err = Plot.yErr(isource).copy()
        ax_spec.errorbar(EkeV, data, xerr=EkeV_err, yerr=data_err,
                         color=color_srcbkg, marker=src_markers[igroup],
                         label=label, linestyle='', zorder=5)
        bkg = (bkg_factors[igroup] *
               np.array(Plot.y(isource+1).copy())).tolist()
        bkg_err = (bkg_factors[igroup] *
                   np.array(Plot.yErr(isource+1).copy())).tolist()
        EkeV_bkg = Plot.x(isource+1).copy()
        EkeV_bkg_err = Plot.xErr(isource+1).copy()
        ax_spec.errorbar(EkeV_bkg, bkg, xerr=EkeV_bkg_err, yerr=bkg_err,
                         color=color_bkg, marker=bkg_markers[igroup],
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
        if n_srcfiles > 1:
            color_srcbkg = colors[igroup]
            color_src = colors[igroup]
            color_bkg = colors[igroup]
        else:
            color_srcbkg = colors[0]
            color_src = colors[2]
            color_bkg = colors[1]
        #src+bkg
        models = []
        bands = []
        Plot.add = False
        Plot.background = False
        if igroup == 0:
            output['test'] = posterior_predictions_plot(analyser, plottype='data', nsamples=100, group=2*igroup+1)
        for content in posterior_predictions_plot(analyser, plottype='data',
                                                  nsamples=100,
                                                  group=2*igroup+1):
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
            lineargs = dict(drawstyle='steps', color=color_srcbkg, zorder=3,
                            linewidth=0.8, linestyle=src_linestyles[igroup])
            shadeargs = dict(color=lineargs['color'], zorder=2)
            band.shade(alpha=0.5, **shadeargs)
            shadeargs = dict(**shadeargs, hatch=hatches[igroup])
            band.shade(q=0.9973/2, alpha=0.2, **shadeargs)
            band.line(label=label, **lineargs)

        #src
        if plot_src:
            models = []
            bands = []
            Plot.add = False
            Plot.background = False
            posterior_backup = analyser.posterior.copy()
            for i_line in range(analyser.posterior):
                analyser.posterior[i_line][ntransf+igroup] = -10
            for content in posterior_predictions_plot(analyser,
                                                      plottype='data',
                                                      nsamples=100,
                                                      group=2*igroup+1):
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
                lineargs = dict(drawstyle='steps', color=color_srcbkg,
                                zorder=3, linewidth=0.8,
                                linestyle=src_linestyles[igroup])
                shadeargs = dict(color=lineargs['color'], zorder=2)
                band.shade(alpha=0.5, **shadeargs)
                shadeargs = dict(**shadeargs, hatch=hatches[igroup])
                band.shade(q=0.9973/2, alpha=0.2, **shadeargs)
                band.line(label=label, **lineargs)
            print('here')
            AllModels.show()
            analyser.posterior = posterior_backup.copy()
            del posterior_backup
            analyser.set_best_fit()
            AllModels.show()

        #bkg
        models = []
        bands = []
        Plot.add = False
        Plot.background = False
        for content in posterior_predictions_plot(analyser, plottype='data',
                                                  nsamples=100,
                                                  group=2*igroup+2):
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
            lineargs = dict(drawstyle='steps', color=color_bkg, zorder=1,
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
                fac_bkg.link = fac_src #this is correct, more
                # information when fitting both at the same time
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
    luminosity_chains = []
    for ispec in range(n_srcfiles):
        src = srcs[ispec]
        fluxes = []
        lums = []
        for band in E_ranges:
            old_nh = []
            flux = analyser.create_flux_chain(src, erange=f'{band[0]}'
                                              f' {band[1]}')[:,0]
            fluxes_band = [np.percentile(flux, quantiles[0]),
                           np.percentile(flux, quantiles[1]),
                           np.percentile(flux, quantiles[2])]
            fluxes.append(fluxes_band)
            for nH in nHs_froz:
                old_nh.append(nH.values)
                nH.values = 1e-5
            for nH in nHs_mod:
                nH.values = 1e-5
            AllModels.show()
            flux = analyser.create_flux_chain(src, erange=f'{band[0]}'
                                              f' {band[1]}')[:,0]
            lums_band = [lum(np.log10(np.percentile(flux, quantiles[0])),
                             distance),
                         lum(np.log10(np.percentile(flux, quantiles[1])),
                             distance),
                         lum(np.log10(np.percentile(flux, quantiles[2])),
                             distance)]
            lums.append(lums_band)
            luminosity_chains.append(lum(np.log10(flux), distance))
            if band == E_ranges[0] and ispec == 0:
                AllModels.show()
            for inH, nH in enumerate(nHs_froz):
                nH.values = old_nh[inH]
            analyser.set_best_fit()
        absorbed_F.append(fluxes)
        unabsorbed_L.append(lums)

    fig_lum = plot_corner_flux(analyser, luminosity_chains, ntransf)
    fig_lum.savefig(f'{working_dir}/corner_lums.pdf')
        
    # this is just to check if analyser.set_best_fit() does its job
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

def plot_corner_flux(analyser, lum_chains, ntransf) -> Figure:
    '''Make a corner plot with fluxes instead of norms.'''
    paramnames = analyser.results['paramnames'][:ntransf]
    data = np.array(analyser.posterior.T[:ntransf])
    i_norm = 1
    for i_name, entry in enumerate(paramnames):
        if entry.lower().find('norm') >= 0:
            log_scale = int(np.log10(np.median(lum_chains[i_norm-1])))
            paramnames[i_name] = 'L$_{%s}$ (10$^{%s}$ erg/s)' % \
                (i_norm, log_scale)
            data[i_name] = lum_chains[i_norm-1] / 10**log_scale
            i_norm += 1
    data = data.T
    fig = corner.corner(data, labels=paramnames, show_titles=True, quiet=True)
    return fig

def write_tex(tex_file, tex_info, abs_F, unabs_L, analyser, quantiles,
              E_ranges):
    tex_file.write('\\begin{{tabular}}{{%s}}\n' % ('c' * (1+len(tex_info)
                                                          +2*len(E_ranges))))
    tex_file.write('\\hline\\hline\n')
    line_1 = ''
    for i in range(len(tex_info)):
        line_1 += tex_info[i][0]
    line_2 = ''
    for i in range(len(tex_info)):
        line_2 += tex_info[i][1]
    for entry in E_ranges:
        line_1 += (' & \\mbox{{F$_{{\\rm x; %s-%s}}$}} '
                   '& \\mbox{{L$_{{\\rm x; %s-%s}}$}}' % (entry[0], entry[1],
                                                          entry[0], entry[1]))
        line_2 += (' & $\\times$erg cm$^{{-2}}$s$^{{-1}}$ & '
                   '$\\times$erg s$^{{-1}}$')
    tex_file.write(f'Part & {line_1} \\\\ \n')
    tex_file.write(f'-- & {line_2} \\\\ \n')
    tex_file.write('\\hline\n')
    for i in range(len(abs_F)):
        tex_file.write('%s\\\\ \n' % ('& '*(len(tex_info) + 2*len(E_ranges))))
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
            line += f'{median}$^{{+{upper-median}}}_{{-{median-lower}}}$'
        for irange in range(len(E_ranges)):
            lower = abs_F[i][irange][0]
            median = abs_F[i][irange][1]
            upper = abs_F[i][irange][2]
            line += f' & {median}$^{{+{upper-median}}}_{{-{median-lower}}}$'
        for irange in range(len(E_ranges)):
            lower = unabs_L[i][irange][0]
            median = unabs_L[i][irange][1]
            upper = unabs_L[i][irange][2]
            line += f' & {median}$^{{+{upper-median}}}_{{-{median-lower}}}$'
        line += '\\\\ \n'
        tex_file.write(line)
    tex_file.write('%s\\\\ \n' % ('& '*(len(tex_info) + 2*len(E_ranges))))

def setup_axis(Emin, Emax, figsize):
    fig = plt.figure(figsize=(figsize[0], figsize[1]))
    ax_spec = fig.add_subplot(111)
    ax_res = fig.add_subplot(111)
    ax_spec.set_yscale('log')

    xticks = []
    xlabels = []
    for tick in [0.5, 1, 5, 10]:
        if tick>Emin and tick<Emax:
            xticks.append(tick)
            xlabels.append(num2text(tick))

    xminorticks = []
    xminorlabels = []
    xminorticks.append(Emin)
    xminorlabels.append(num2text(Emin))
    for tick in [0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 2, 3, 4, 6, 7, 8, 9]:
        if tick>Emin and tick<Emax:
            xminorticks.append(tick)
            xminorlabels.append('')
    xminorticks.append(Emax)
    xminorlabels.append(num2text(Emax))

    for ax in [ax_spec, ax_res]:
        ax.set_xscale('log')

        ax.grid(which='major', axis='both', zorder=-1,
                color='#111111', linewidth=0.8)
        ax.grid(which='minor', axis='both', zorder=-1,
                color='#222222', linewidth=0.5, linestyle=':')

        ax.tick_params(axis='x', which='major', direction='in',
                        top='on', pad=9, length=5, width=1.5)  # , labelsize=10)
        ax.tick_params(axis='x', which='minor', direction='in',
                        top='on', pad=9, length=3)  # , labelsize=0)
        ax.tick_params(axis='y', which='major', direction='in',
                        right='on', length=5, width=1.5)  # , labelsize=10)
        ax.tick_params(axis='y', which='minor', direction='in',
                        right='on', length=3)  # , labelsize=0)
        ax.set_xticks(xticks, minor=False)
        ax.set_xticks(xminorticks, minor=True)
        ax.set_xticklabels(xminorlabels, minor=True)

    ax_res.set_yticks([0], minor=False)
    ax_res.set_yticks([-4, -2, 2, 4], minor=True)
    ax_res.set_yticklabels(['-4', '-2', '2', '4'], minor=True)
    ax_res.set_xticklabels(xlabels, minor=False)

    ax_spec.set_xticklabels(['']*len(xticks), minor = False)
    ax_spec.set_xticklabels(['']*len(xminorticks), minor=True)

    ax_res_invis = fig.add_subplot(111)
    ax_res_invis.set_frame_on(False)
    ax_res_invis.patch.set_facecolor("none")

    return fig, ax_spec, ax_res, ax_res_invis

def format_axis_pt2(fig, ax_spec, ax_res, ax_res_invis, fig_borders, rescale_F,
                    rescale_chi, E_ranges, src_files, ncols, nrows,
                    height_ratios, width_ratios):
    fig.canvas.draw()
    fig.set_tight_layout(True)
    fig.set_tight_layout(False)
    if len(src_files) > 1:
        ax_spec.legend(bbox_to_anchor=(-0.044, 1.02), loc='upper right',
                        handletextpad=0.1, fontsize=12)
    # hspace = 8.0 / figsize[1] * 0.05
    hspace = 0
    fig.subplots_adjust(
        hspace=hspace, top=fig_borders[0], bottom=fig_borders[1],
        left=fig_borders[2], right=fig_borders[3])

    ax_res_invis.tick_params(left=True, bottom=True,
                                right=True, top=True)
    ax_res_invis.tick_params(axis='x', which='major', direction='in',
                                top='on',   pad=9, length=0)  # , labelsize=10)
    ax_res_invis.tick_params(axis='x', which='minor', direction='in',
                                top='on',   length=0)  # , labelsize=0)
    ax_res_invis.tick_params(axis='y', which='major', direction='in',
                                right='on', length=0)  # , labelsize=10)
    ax_res_invis.tick_params(axis='y', which='minor', direction='in',
                                right='on', length=0)

    ax_res_invis.set_xbound(lower=0, upper=1)
    ax_res_invis.set_xticks([0, 1])
    ax_res_invis.set_xticklabels(['0', '1'], alpha=0)

    ax_res_invis.set_ybound(lower=rescale_F[0], upper=rescale_F[1])
    ax_res_invis.set_yticks(ax_spec.get_yticks())
    ax_res_invis.set_yticklabels(ax_spec.get_yticklabels(), alpha=0)

    gs = gridspec.GridSpec(ncols=ncols,
                            nrows=nrows,
                            height_ratios=height_ratios,
                            width_ratios=width_ratios)

    ax_spec.set_position(gs[0].get_position(fig))
    ax_res.set_position(gs[1].get_position(fig))
    ax_res_invis.set_position(gs[1].get_position(fig))

    ax_spec.set_ylabel('Counts s$^{-1}$ keV$^{-1}$')
    ax_res_invis.set_ylabel('$\\Delta\\chi$')
    ax_res.set_xlabel('Energy (keV)')

    ax_spec.set_ybound(lower=rescale_F[0], upper=rescale_F[1])
    ax_res.set_ybound(lower=rescale_chi[0], upper=rescale_chi[1])

    ax_spec.set_xbound(lower=E_ranges[0][0], upper=E_ranges[0][1])
    ax_res.set_xbound(lower=E_ranges[0][0], upper=E_ranges[0][1])