#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:23:18 2022

@author: David Kaltenbrunner
"""
from astropy.stats import bayesian_blocks
from cmdstanpy import cmdstan_path, install_cmdstan, CmdStanModel
import numpy as np

from HiMaXBipy.io.package_data import round_to_1
import logging

# checking if cmdstan is already installed, otherwise installing it
try:
    cmdstan_path()
except ValueError:
    print('Installing CmdStan')
    install_cmdstan()


def plot_lc_eROday_broken_bayes_old(hdulist, axs, log, mjdref, xflag,
                                    color, obs_periods, short_time, stan_model,
                                    quantiles, time_rel=0, fexp_cut=0.15,
                                    alpha_bg=0.5):
    '''
    Lightcurve rebinned to eROdays with countrates optained with Bayesian fit
    assuming Poissionian distribution for counts and log
    '''
    logger_stan = logging.getLogger('cmdstanpy')
    logger_stan.setLevel(logging.DEBUG)
    logger_stan.handlers = log.handlers.copy()
    pxmin = []
    pxmax = []
    ymin = 0
    pymax = []
    sc_rate = []
    sc_rate_upper = []
    sc_rate_lower = []
    bg_rate = []
    bg_rate_upper = []
    bg_rate_lower = []
    xtime = []
    xtime_d = []
    fexp_full = hdulist[1].data.field('FRACEXP')
    time = hdulist[1].data.field('TIME')[fexp_full > fexp_cut]
    time_mjd = time / 3600. / 24. + mjdref
    delt = hdulist[1].data.field('TIMEDEL')[fexp_full > fexp_cut]
    cnts = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[fexp_full > fexp_cut]
    fexp = hdulist[1].data.field('FRACEXP')[fexp_full > fexp_cut]
    back = np.array(hdulist[1].data.field('BACK_COUNTS'), dtype=int)[
        fexp_full > fexp_cut]
    backrat = hdulist[1].data.field('BACKRATIO')[fexp_full > fexp_cut]
    for i, entry in enumerate(backrat):
        if entry < 0.01:
            backrat[i] = 0.01

    # loading stan model
    model = CmdStanModel(stan_file=stan_model)

    istart = 0
    iend = 0
    tstart = 0
    tend = 0
    nrow = len(time)
    # rebinning in scans and getting sc and bg rates with uncertainties
    # from quantiles
    for i in range(nrow):
        if i == nrow - 1:
            iend = i + 1
        elif time[i + 1] - time[i] > 3600.0:  # elif to avoid error
            iend = i + 1
        else:
            continue
        if istart == 0:
            tstart = time[0]
        else:
            border_low = False
            for period in obs_periods:
                if time_mjd[istart-1] < period[0] and time_mjd[istart] > period[0]:
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow:
            tend = time[nrow - 1]
        else:
            border_high = False
            for period in obs_periods:
                if time_mjd[iend] < period[1] and time_mjd[iend + 1] > period[1]:
                    tend = time[iend] + delt[iend]
                    border_high = True
            if not border_high:
                tend = (time[iend] + time[iend+1]) / 2
        xtime.append((tend+tstart)/2)
        xtime_d.append((tend-tstart)/2)
        data = {}
        data['N'] = iend-istart
        data['dt'] = delt[istart:iend]
        data['sc'] = cnts[istart:iend]
        data['frac_exp'] = fexp[istart:iend]
        data['bg'] = back[istart:iend]
        data['bg_ratio'] = backrat[istart:iend]
        fit = model.sample(data=data, show_progress=False)
        sc_rate_lower.append(np.percentile(fit.stan_variables()['sc_rate'],
                                           quantiles[0]))
        sc_rate.append(np.percentile(fit.stan_variables()['sc_rate'],
                                     quantiles[1]))
        sc_rate_upper.append(np.percentile(fit.stan_variables()['sc_rate'],
                                           quantiles[2]))
        bg_rate_lower.append(np.percentile(fit.stan_variables()['bg_rate'],
                                           quantiles[0]))
        bg_rate.append(np.percentile(fit.stan_variables()['bg_rate'],
                                     quantiles[1]))
        bg_rate_upper.append(np.percentile(fit.stan_variables()['bg_rate'],
                                           quantiles[2]))
        istart = i + 1

    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    xtime = np.array(xtime)
    xtime_d = np.array(xtime_d)
    mjd = xtime * (1. / 24. / 3600.) + mjdref
    mjd_d = xtime_d / 24. / 3600.
    sc_rate = np.array(sc_rate)
    sc_rate_lower = np.array(sc_rate_lower)
    sc_rate_upper = np.array(sc_rate_upper)
    bg_rate = np.array(bg_rate)
    bg_rate_lower = np.array(bg_rate_lower)
    bg_rate_upper = np.array(bg_rate_upper)

    for i_ax, ax in enumerate(axs):
        if xflag == 1:
            xtime_part = xtime[(mjd > obs_periods[i_ax][0]) *
                               (mjd < obs_periods[i_ax][1])]
            xmin = min(xtime_part)
            xmax = max(xtime_part)
        else:
            mjd_part = mjd[(mjd > obs_periods[i_ax][0]) *
                           (mjd < obs_periods[i_ax][1])]
            xmin = min(mjd_part)
            xmax = max(mjd_part)

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(xmin)
            xtime_short = xtime - time_rel
            mjd_short = mjd - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel
        else:
            mjd_short = mjd.copy()
            xtime_short = xtime.copy()

        xm = (xmax-xmin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)

        if xflag == 1:
            ax.errorbar(xtime_short, sc_rate, xerr=xtime_d,
                        yerr=[sc_rate + (-sc_rate_lower),
                              sc_rate_upper + (-sc_rate)],
                        linestyle='None', color=color, fmt='o',
                        zorder=1)
            ax.errorbar(xtime_short, bg_rate, xerr=xtime_d,
                        yerr=[bg_rate + (-bg_rate_lower),
                              bg_rate_upper + (-bg_rate)],
                        linestyle='None', color=color, fmt='x',
                        zorder=1, alpha=alpha_bg)
        else:
            ax.errorbar(mjd_short, sc_rate, xerr=mjd_d,
                        yerr=[sc_rate + (-sc_rate_lower),
                              sc_rate_upper + (-sc_rate)],
                        linestyle='None', color=color, fmt='o',
                        zorder=1)
            ax.errorbar(mjd_short, bg_rate, xerr=mjd_d,
                        yerr=[bg_rate + (-bg_rate_lower),
                              bg_rate_upper + (-bg_rate)],
                        linestyle='None', color=color, fmt='x',
                        zorder=1, alpha=alpha_bg)

    ymax = max([max(sc_rate_upper), max(bg_rate_upper)])
    pymin = ymin - (ymax-ymin)*0.05
    pymax = ymax + (ymax-ymin)*0.05

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_eROday_broken_bayes(hdulist, axs, log, mjdref, xflag,
                                color, obs_periods, short_time, stan_model,
                                quantiles, time_rel=0, fexp_cut=0.15,
                                alpha_bg=0.5, bblocks=False, bbp0=0.01,
                                bbmode='both'):
    '''
    Lightcurve rebinned to eROdays with countrates optained with
    Bayesian fit assuming Poissionian distribution for counts and log;
    fit done for each bin simultaneously
    '''
    logger_stan = logging.getLogger('cmdstanpy')
    logger_stan.setLevel(logging.DEBUG)
    logger_stan.handlers = log.handlers.copy()
    pxmin = []
    pxmax = []
    ymin = 0
    pymax = []
    sc_rate = []
    sc_rate_upper = []
    sc_rate_lower = []
    bg_rate = []
    bg_rate_upper = []
    bg_rate_lower = []
    xtime = []
    xtime_d = []
    fexp_full = hdulist[1].data.field('FRACEXP')
    time = hdulist[1].data.field('TIME')[fexp_full > fexp_cut]
    time_mjd = time / 3600. / 24. + mjdref
    delt = hdulist[1].data.field('TIMEDEL')[fexp_full > fexp_cut]
    cnts = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[fexp_full > fexp_cut]
    fexp = hdulist[1].data.field('FRACEXP')[fexp_full > fexp_cut]
    back = np.array(hdulist[1].data.field('BACK_COUNTS'), dtype=int)[
        fexp_full > fexp_cut]
    backrat = hdulist[1].data.field('BACKRATIO')[fexp_full > fexp_cut]
    for i, entry in enumerate(backrat):
        if entry < 0.01:
            backrat[i] = 0.01

    dts = []
    scs = []
    fexps = []
    bgs = []
    bgrats = []

    istart = 0
    iend = 0
    tstart = 0
    tend = 0
    nrow = len(time)
    # rebinning in scans and getting sc and bg rates with uncertainties
    # from quantiles
    for i in range(nrow):
        if i == nrow - 1:
            iend = i
        elif time[i + 1] - time[i] > 3600.0:  # elif to avoid error
            iend = i
        else:
            continue
        if istart == 0:
            tstart = time[0]
        else:
            border_low = False
            for period in obs_periods:
                if time_mjd[istart-1] < period[0] and time_mjd[istart] > period[0]:
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow - 1:
            tend = time[iend] + delt[iend]
        else:
            border_high = False
            for period in obs_periods:
                if time_mjd[iend] < period[1] and time_mjd[iend + 1] > period[1]:
                    tend = time[iend] + delt[iend]
                    border_high = True
            if not border_high:
                tend = (time[iend] + time[iend+1]) / 2
        xtime.append((tend+tstart)/2)
        xtime_d.append((tend-tstart)/2)
        # data['N'] = iend-istart
        dts.append(delt[istart:iend+1].tolist())
        scs.append(cnts[istart:iend+1].tolist())
        fexps.append(fexp[istart:iend+1].tolist())
        bgs.append(back[istart:iend+1].tolist())
        bgrats.append(backrat[istart:iend+1].tolist())
        istart = i + 1

    N = 0
    for k in range(len(scs)):
        if len(scs[k]) > N:
            N = len(scs[k])
    for k in range(len(scs)):
        while len(scs[k]) < N:
            dts[k].append(0)
            scs[k].append(0)
            fexps[k].append(0)
            bgs[k].append(0)
            bgrats[k].append(0)

    data = {}
    data['N'] = N
    data['M'] = len(scs)
    data['dt'] = np.array(dts)
    data['sc'] = np.array(scs)
    data['frac_exp'] = np.array(fexps)
    data['bg'] = np.array(bgs)
    data['bg_ratio'] = np.array(bgrats)

    # loading stan model
    model = CmdStanModel(stan_file=stan_model)
    fit = model.sample(data=data, show_progress=False)
    sc_rate_lower = np.percentile(fit.stan_variables()['sc_rate'],
                                  quantiles[0], axis=0)
    sc_rate = np.percentile(fit.stan_variables()['sc_rate'],
                            quantiles[1], axis=0)
    sc_rate_upper = np.percentile(fit.stan_variables()['sc_rate'],
                                  quantiles[2], axis=0)
    bg_rate_lower = np.percentile(fit.stan_variables()['bg_rate'],
                                  quantiles[0], axis=0)
    bg_rate = np.percentile(fit.stan_variables()['bg_rate'],
                            quantiles[1], axis=0)
    bg_rate_upper = np.percentile(fit.stan_variables()['bg_rate'],
                                  quantiles[2], axis=0)
    sc_bg_rate_lower = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                     quantiles[0], axis=0)
    sc_bg_rate = np.percentile(fit.stan_variables()['sc_bg_rate'],
                               quantiles[1], axis=0)
    sc_bg_rate_upper = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                     quantiles[2], axis=0)

    amp_var = np.percentile(fit.stan_variables()['amp_dev'], quantiles[1])
    amp_var_low = np.percentile(fit.stan_variables()['amp_dev'], quantiles[0])
    amp_var_up = np.percentile(fit.stan_variables()['amp_dev'], quantiles[2])
    amp_frac = np.percentile(fit.stan_variables()['amp_frac'], quantiles[1])
    amp_frac_low = np.percentile(
        fit.stan_variables()['amp_frac'], quantiles[0])
    amp_frac_up = np.percentile(
        fit.stan_variables()['amp_frac'], quantiles[2])
    logger_stan.warning(f'eROday AmpVar={amp_var}+{amp_var_up-amp_var}'
                        f'-{amp_var-amp_var_low}\n')
    logger_stan.warning(f'AmpFrac={amp_frac}'
                        f'+{amp_frac_up-amp_frac}'
                        f'-{amp_frac-amp_frac_low}\n')
    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    xtime = np.array(xtime)
    xtime_d = np.array(xtime_d)
    mjd = xtime * (1. / 24. / 3600.) + mjdref
    mjd_d = xtime_d / 24. / 3600.
    sc_rate = np.array(sc_rate)
    sc_rate_lower = np.array(sc_rate_lower)
    sc_rate_upper = np.array(sc_rate_upper)
    bg_rate = np.array(bg_rate)
    bg_rate_lower = np.array(bg_rate_lower)
    bg_rate_upper = np.array(bg_rate_upper)
    sc_bg_rate = np.array(sc_bg_rate)
    sc_bg_rate_lower = np.array(sc_bg_rate_lower)
    sc_bg_rate_upper = np.array(sc_bg_rate_upper)
    sc_rate_err = [sc_rate + (-sc_rate_lower),
                   sc_rate_upper + (-sc_rate)]
    bg_rate_err = [bg_rate + (-bg_rate_lower),
                   bg_rate_upper + (-bg_rate)]
    sc_bg_rate_err = [sc_bg_rate + (-sc_bg_rate_lower),
                      sc_bg_rate_upper + (-sc_bg_rate)]
    t = xtime.copy()

    for i_ax, ax in enumerate(axs):
        if xflag == 1:
            xtime_part = xtime[(mjd > obs_periods[i_ax][0]) *
                               (mjd < obs_periods[i_ax][1])]
            xtime_d_part = xtime_d[(mjd > obs_periods[i_ax][0]) *
                                   (mjd < obs_periods[i_ax][1])]
            xmin = min(xtime_part - xtime_d_part)
            xmax = max(xtime_part + xtime_d_part)
        else:
            mjd_part = mjd[(mjd > obs_periods[i_ax][0]) *
                           (mjd < obs_periods[i_ax][1])]
            mjd_d_part = mjd_d[(mjd > obs_periods[i_ax][0]) *
                               (mjd < obs_periods[i_ax][1])]
            xmin = min(mjd_part - mjd_d_part)
            xmax = max(mjd_part + mjd_d_part)

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(xmin)
            xtime_short = xtime - time_rel
            mjd_short = mjd - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel
        else:
            mjd_short = mjd.copy()
            xtime_short = xtime.copy()

        xm = (xmax-xmin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)

        if xflag == 1:
            t = xtime_short
            terr = xtime_d
        else:
            t = mjd_short
            terr = mjd_d
        ax.errorbar(t, sc_rate, xerr=terr,
                    yerr=sc_rate_err,
                    linestyle='None', color=color, fmt='o',
                    zorder=4)
        ax.errorbar(t, bg_rate, xerr=terr,
                    yerr=bg_rate_err,
                    linestyle='None', color=color, fmt='x',
                    zorder=3, alpha=alpha_bg)

    ymax = max([max(sc_rate_upper), max(bg_rate_upper)])
    pymin = ymin - (ymax-ymin)*0.05
    pymax = ymax + (ymax-ymin)*0.05

    if bblocks:
        if bbmode == 'sc' or bbmode == 'both':
            err = np.max(sc_rate_err, axis=0)
            err = np.min([sc_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            sc_rate_bb = np.percentile(fit.stan_variables()['sc_rate'],
                                       quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(sc_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')
        elif bbmode == 'sum':
            err = np.max(sc_bg_rate_err, axis=0)
            err = np.min([sc_bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_bg_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            sc_bg_rate_bb = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                          quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(sc_bg_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')
        if bbmode == 'both':
            err = np.max(bg_rate_err, axis=0)
            err = np.min([bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, bg_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            bg_rate_bb = np.percentile(fit.stan_variables()['bg_rate'],
                                       quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(bg_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_mincounts_broken_bayes(hdulist, axs, log, mjdref, xflag,
                                   mincounts, color, obs_periods,
                                   short_time, stan_model, quantiles,
                                   time_rel=0, fexp_cut=0.15,
                                   alpha_bg=0.5, bblocks=False, bbp0=0.01,
                                   bbmode='both'):
    '''
    Lightcurve rebinned to eROdays with countrates optained with
    Bayesian fit assuming Poissionian distribution for counts and log;
    fit done for each bin simultaneously
    '''
    logger_stan = logging.getLogger('cmdstanpy')
    logger_stan.setLevel(logging.DEBUG)
    logger_stan.handlers = log.handlers.copy()
    pxmin = []
    pxmax = []
    ymin = 0
    pymax = []
    sc_rate = []
    sc_rate_upper = []
    sc_rate_lower = []
    bg_rate = []
    bg_rate_upper = []
    bg_rate_lower = []
    xtime = []
    xtime_d = []
    fexp_full = hdulist[1].data.field('FRACEXP')
    time = hdulist[1].data.field('TIME')[fexp_full > fexp_cut]
    time_mjd = time / 3600. / 24. + mjdref
    delt = hdulist[1].data.field('TIMEDEL')[fexp_full > fexp_cut]
    cnts = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[fexp_full > fexp_cut]
    fexp = hdulist[1].data.field('FRACEXP')[fexp_full > fexp_cut]
    back = np.array(hdulist[1].data.field('BACK_COUNTS'), dtype=int)[
        fexp_full > fexp_cut]
    backrat = hdulist[1].data.field('BACKRATIO')[fexp_full > fexp_cut]
    for i, entry in enumerate(backrat):
        if entry < 0.01:
            backrat[i] = 0.01

    dts = []
    scs = []
    fexps = []
    bgs = []
    bgrats = []

    istart = 0
    iend = 0
    tstart = 0
    tend = 0
    istart_old = 0
    nrow = len(time)
    # rebinning in scans and getting sc and bg rates with uncertainties
    # from quantiles
    for i in range(nrow):
        if i == nrow - 1:
            iend = i
            if (cnts[istart:nrow].sum() < mincounts and
                    not time_mjd[istart_old] < obs_periods[-1][0]):
                istart = istart_old
                del (xtime[-1], xtime_d[-1], dts[-1], scs[-1], fexps[-1],
                     bgs[-1], bgrats[-1])
        elif (cnts[istart:i].sum() >= mincounts and
              time[i + 1] - time[i] > 3600.0):  # elif to avoid error
            iend = i
        else:
            period_last = False
            period_first = False
            for period in obs_periods:
                if time_mjd[i] < period[1] and time_mjd[i+1] > period[1]:
                    period_last = True
                    if time_mjd[istart_old] < period[0] or istart == 0:
                        period_first = True
            if period_last:
                iend = i
                if not period_first:
                    istart = istart_old
                    del (xtime[-1], xtime_d[-1], dts[-1], scs[-1], fexps[-1],
                         bgs[-1], bgrats[-1])
            else:
                continue

        if istart == 0:
            tstart = time[0]
        else:
            border_low = False
            for period in obs_periods:
                if time_mjd[istart-1] < period[0] and time_mjd[istart] > period[0]:
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow - 1:
            tend = time[iend] + delt[iend]
        else:
            border_high = False
            for period in obs_periods:
                if time_mjd[iend] < period[1] and time_mjd[iend + 1] > period[1]:
                    tend = time[iend] + delt[iend]
                    border_high = True
            if not border_high:
                tend = (time[iend] + time[iend+1]) / 2
        xtime.append((tend+tstart)/2)
        xtime_d.append((tend-tstart)/2)
        # data['N'] = iend-istart
        dts.append(delt[istart:iend+1].tolist())
        scs.append(cnts[istart:iend+1].tolist())
        fexps.append(fexp[istart:iend+1].tolist())
        bgs.append(back[istart:iend+1].tolist())
        bgrats.append(backrat[istart:iend+1].tolist())
        istart_old = istart
        istart = i + 1

    N = 0
    for k in range(len(scs)):
        if len(scs[k]) > N:
            N = len(scs[k])
    for k in range(len(scs)):
        while len(scs[k]) < N:
            dts[k].append(0)
            scs[k].append(0)
            fexps[k].append(0)
            bgs[k].append(0)
            bgrats[k].append(0)

    data = {}
    data['N'] = N
    data['M'] = len(scs)
    data['dt'] = np.array(dts)
    data['sc'] = np.array(scs)
    data['frac_exp'] = np.array(fexps)
    data['bg'] = np.array(bgs)
    data['bg_ratio'] = np.array(bgrats)

    # loading stan model
    model = CmdStanModel(stan_file=stan_model)
    fit = model.sample(data=data, show_progress=False)
    sc_rate_lower = np.percentile(fit.stan_variables()['sc_rate'],
                                  quantiles[0], axis=0)
    sc_rate = np.percentile(fit.stan_variables()['sc_rate'],
                            quantiles[1], axis=0)
    sc_rate_upper = np.percentile(fit.stan_variables()['sc_rate'],
                                  quantiles[2], axis=0)
    bg_rate_lower = np.percentile(fit.stan_variables()['bg_rate'],
                                  quantiles[0], axis=0)
    bg_rate = np.percentile(fit.stan_variables()['bg_rate'],
                            quantiles[1], axis=0)
    bg_rate_upper = np.percentile(fit.stan_variables()['bg_rate'],
                                  quantiles[2], axis=0)
    sc_bg_rate_lower = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                     quantiles[0], axis=0)
    sc_bg_rate = np.percentile(fit.stan_variables()['sc_bg_rate'],
                               quantiles[1], axis=0)
    sc_bg_rate_upper = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                     quantiles[2], axis=0)
    t = xtime.copy()

    amp_var = np.percentile(fit.stan_variables()['amp_dev'], quantiles[1])
    amp_var_low = np.percentile(fit.stan_variables()['amp_dev'], quantiles[0])
    amp_var_up = np.percentile(fit.stan_variables()['amp_dev'], quantiles[2])
    amp_frac = np.percentile(fit.stan_variables()['amp_frac'], quantiles[1])
    amp_frac_low = np.percentile(
        fit.stan_variables()['amp_frac'], quantiles[0])
    amp_frac_up = np.percentile(
        fit.stan_variables()['amp_frac'], quantiles[2])
    logger_stan.warning(f'mincounts {mincounts} AmpVar={amp_var}'
                        f'+{amp_var_up-amp_var}'
                        f'-{amp_var-amp_var_low}\n')
    logger_stan.warning(f'AmpFrac={amp_frac}'
                        f'+{amp_frac_up-amp_frac}'
                        f'-{amp_frac-amp_frac_low}\n')

    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    xtime = np.array(xtime)
    xtime_d = np.array(xtime_d)
    mjd = xtime * (1. / 24. / 3600.) + mjdref
    mjd_d = xtime_d / 24. / 3600.
    sc_rate = np.array(sc_rate)
    sc_rate_lower = np.array(sc_rate_lower)
    sc_rate_upper = np.array(sc_rate_upper)
    bg_rate = np.array(bg_rate)
    bg_rate_lower = np.array(bg_rate_lower)
    bg_rate_upper = np.array(bg_rate_upper)
    sc_bg_rate = np.array(sc_bg_rate)
    sc_bg_rate_lower = np.array(sc_bg_rate_lower)
    sc_bg_rate_upper = np.array(sc_bg_rate_upper)
    sc_rate_err = [sc_rate + (-sc_rate_lower),
                   sc_rate_upper + (-sc_rate)]
    bg_rate_err = [bg_rate + (-bg_rate_lower),
                   bg_rate_upper + (-bg_rate)]
    sc_bg_rate_err = [sc_bg_rate + (-sc_bg_rate_lower),
                      sc_bg_rate_upper + (-sc_bg_rate)]

    for i_ax, ax in enumerate(axs):
        if xflag == 1:
            xtime_part = xtime[(mjd > obs_periods[i_ax][0]) *
                               (mjd < obs_periods[i_ax][1])]
            xtime_d_part = xtime_d[(mjd > obs_periods[i_ax][0]) *
                                   (mjd < obs_periods[i_ax][1])]
            xmin = min(xtime_part - xtime_d_part)
            xmax = max(xtime_part + xtime_d_part)
        else:
            mjd_part = mjd[(mjd > obs_periods[i_ax][0]) *
                           (mjd < obs_periods[i_ax][1])]
            mjd_d_part = mjd_d[(mjd > obs_periods[i_ax][0]) *
                               (mjd < obs_periods[i_ax][1])]
            xmin = min(mjd_part - mjd_d_part)
            xmax = max(mjd_part + mjd_d_part)

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(xmin)
            xtime_short = xtime - time_rel
            mjd_short = mjd - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel
        else:
            mjd_short = mjd.copy()
            xtime_short = xtime.copy()

        xm = (xmax-xmin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)

        if xflag == 1:
            t = xtime_short
            terr = xtime_d
        else:
            t = mjd_short
            terr = mjd_d
        ax.errorbar(t, sc_rate, xerr=terr,
                    yerr=sc_rate_err,
                    linestyle='None', color=color, fmt='o',
                    zorder=4)
        ax.errorbar(t, bg_rate, xerr=terr,
                    yerr=bg_rate_err,
                    linestyle='None', color=color, fmt='x',
                    zorder=3, alpha=alpha_bg)

    ymax = max([max(sc_rate_upper), max(bg_rate_upper)])
    pymin = ymin - (ymax-ymin)*0.05
    pymax = ymax + (ymax-ymin)*0.05

    if bblocks:
        if bbmode == 'sc' or bbmode == 'both':
            err = np.max(sc_rate_err, axis=0)
            err = np.min([sc_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            sc_rate_bb = np.percentile(fit.stan_variables()['sc_rate'],
                                       quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(sc_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')
        elif bbmode == 'sum':
            err = np.max(sc_bg_rate_err, axis=0)
            err = np.min([sc_bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_bg_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            sc_bg_rate_bb = np.percentile(fit.stan_variables()['sc_bg_rate'],
                                          quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(sc_bg_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')
        if bbmode == 'both':
            err = np.max(bg_rate_err, axis=0)
            err = np.min([bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, bg_rate, sigma=err,
                                        fitness='measure', p0=bbp0)
            dts_bb = []
            scs_bb = []
            fexps_bb = []
            bgs_bb = []
            bgrats_bb = []
            for i in range(len(t_blocked) - 1):
                if xflag == 1:
                    tstart = t_blocked[i] + time_rel
                    tend = t_blocked[i+1] + time_rel
                else:
                    tstart = (t_blocked[i] + time_rel - mjdref) * 24. * 3600.
                    tend = (t_blocked[i+1] + time_rel - mjdref) * 24. * 3600.
                dts_bb.append(delt[(time >= tstart)*(time <= tend)].tolist())
                scs_bb.append(cnts[(time >= tstart)*(time <= tend)].tolist())
                fexps_bb.append(fexp[(time >= tstart)*(time <= tend)].tolist())
                bgs_bb.append(back[(time >= tstart)*(time <= tend)].tolist())
                bgrats_bb.append(
                    backrat[(time >= tstart)*(time <= tend)].tolist())

            N_bb = 0
            for k in range(len(scs_bb)):
                if len(scs_bb[k]) > N_bb:
                    N_bb = len(scs_bb[k])
            for k in range(len(scs_bb)):
                while len(scs_bb[k]) < N_bb:
                    dts_bb[k].append(0)
                    scs_bb[k].append(0)
                    fexps_bb[k].append(0)
                    bgs_bb[k].append(0)
                    bgrats_bb[k].append(0)

            data = {}
            data['N'] = N_bb
            data['M'] = len(scs_bb)
            data['dt'] = np.array(dts_bb)
            data['sc'] = np.array(scs_bb)
            data['frac_exp'] = np.array(fexps_bb)
            data['bg'] = np.array(bgs_bb)
            data['bg_ratio'] = np.array(bgrats_bb)

            # loading stan model
            model = CmdStanModel(stan_file=stan_model)
            fit = model.sample(data=data, show_progress=False)
            bg_rate_bb = np.percentile(fit.stan_variables()['bg_rate'],
                                       quantiles[1], axis=0)
            for ax in axs:
                ax.stairs(bg_rate_bb, t_blocked, color=color, zorder=2,
                          linestyle='--')

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel
