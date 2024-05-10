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
                                    alpha_bg=0.3):
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
                if (time_mjd[istart-1] < period[0]
                        and time_mjd[istart] > period[0]):
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow:
            tend = time[nrow - 1]
        else:
            border_high = False
            for period in obs_periods:
                if (time_mjd[iend] < period[1]
                        and time_mjd[iend + 1] > period[1]):
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
                                alpha_bg=0.3, bblocks=False, bbp0=0.003,
                                bbmode='both', yscale='linear'):
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
    ftime = hdulist[1].data.field('FRACTIME')[fexp_full > fexp_cut]
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
                if (time_mjd[istart-1] < period[0]
                        and time_mjd[istart] > period[0]):
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow - 1:
            tend = time[iend] + delt[iend]
        else:
            border_high = False
            for period in obs_periods:
                if (time_mjd[iend] < period[1]
                        and time_mjd[iend + 1] > period[1]):
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

    amp_dev_min = np.percentile(
        fit.stan_variables()['amp_dev_min'], quantiles[1])
    amp_dev_min_low = np.percentile(
        fit.stan_variables()['amp_dev_min'], quantiles[0])
    amp_dev_min_up = np.percentile(
        fit.stan_variables()['amp_dev_min'], quantiles[2])
    amp_frac_min = np.percentile(
        fit.stan_variables()['amp_frac_min'], quantiles[1])
    amp_frac_min_low = np.percentile(
        fit.stan_variables()['amp_frac_min'], quantiles[0])
    amp_frac_min_up = np.percentile(
        fit.stan_variables()['amp_frac_min'], quantiles[2])
    logger_stan.warning(f'eROday AmpVar(min)={amp_dev_min}+'
                        f'{amp_dev_min_up-amp_dev_min}'
                        f'-{amp_dev_min-amp_dev_min_low}\n')
    logger_stan.warning(f'eROday AmpFrac(min)={amp_frac_min}'
                        f'+{amp_frac_min_up-amp_frac_min}'
                        f'-{amp_frac_min-amp_frac_min_low}\n')

    amp_dev_med = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[1])
    amp_dev_med_low = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[0])
    amp_dev_med_up = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[2])
    amp_frac_med = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[1])
    amp_frac_med_low = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[0])
    amp_frac_med_up = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[2])
    logger_stan.warning(f'eROday AmpVar(med)={amp_dev_med}+'
                        f'{amp_dev_med_up-amp_dev_med}'
                        f'-{amp_dev_med-amp_dev_med_low}\n')
    logger_stan.warning(f'eROday AmpFrac(med)={amp_frac_med}'
                        f'+{amp_frac_med_up-amp_frac_med}'
                        f'-{amp_frac_med-amp_frac_med_low}\n')

    yrate = np.array(sc_rate)
    yrate_lower = np.array(sc_rate_lower)
    yrate_upper = np.array(sc_rate_upper)
    yrate_e = np.sqrt(((yrate-yrate_lower)**2+(yrate_upper-yrate)**2)/2)

    i_max = np.argmax(yrate_lower)
    i_min = np.argmin(yrate_upper)

    ampl_max2 = yrate_lower[i_max] - yrate_upper[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    log.info(f'AMPL_MAX scan: {ampl_max}\n')
    log.info(f'Variability V scan = {variability}\n')
    log.info(f'AMPL_SIG scan: {ampl_sig}\n')
    log.info(f'AMPL_MAX conservative scan: {ampl_max2}\n')
    log.info(f'AMPL_SIG2 scan: {ampl_sig2}\n')

    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    # average rate:
    data_av = {}
    data_av['N'] = len(delt)
    data_av['M'] = 1
    data_av['dt'] = np.array([delt])
    data_av['sc'] = np.array([cnts], dtype=int)
    data_av['frac_exp'] = np.array([fexp])
    data_av['bg'] = np.array([back], dtype=int)
    data_av['bg_ratio'] = np.array([backrat])

    fit_av = model.sample(data=data_av, show_progress=False)
    sc_rate_av_lower = np.percentile(fit_av.stan_variables()['sc_rate'],
                                     quantiles[0], axis=0)[0]
    sc_rate_av = np.percentile(fit_av.stan_variables()['sc_rate'],
                               quantiles[1], axis=0)[0]
    sc_rate_av_upper = np.percentile(fit_av.stan_variables()['sc_rate'],
                                     quantiles[2], axis=0)[0]
    bg_rate_av_lower = np.percentile(fit_av.stan_variables()['bg_rate'],
                                     quantiles[0], axis=0)[0]
    bg_rate_av = np.percentile(fit_av.stan_variables()['bg_rate'],
                               quantiles[1], axis=0)[0]
    bg_rate_av_upper = np.percentile(fit_av.stan_variables()['bg_rate'],
                                     quantiles[2], axis=0)[0]
    logger_stan.warning(f'Average Source Rate={sc_rate_av}'
                        f'+{sc_rate_av_upper-sc_rate_av}'
                        f'-{sc_rate_av-sc_rate_av_lower} cts/s')
    logger_stan.warning(f'Average Background Rate={bg_rate_av}'
                        f'+{bg_rate_av_upper-bg_rate_av}'
                        f'-{bg_rate_av-bg_rate_av_lower} cts/s')
    tot_exp = ftime.sum()
    tot_fexp = (delt * fexp).sum()
    logger_stan.warning(f'Total Exposure={tot_exp}s')
    logger_stan.warning(f'Total Fractional Exposure={tot_fexp}s')

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
    if yscale == 'linear':
        ymin = 0
        pymin = ymin - (ymax-ymin)*0.05
        pymax = ymax + (ymax-ymin)*0.05
    elif yscale == 'log':
        if max(bg_rate_lower) > min(sc_rate_lower) * 1e-1 and alpha_bg > 0:
            ymin = min([min(sc_rate_lower), min(bg_rate_lower)])
        else:
            ymin = min(sc_rate_lower)
        pymin = ymin / ((ymax/ymin) ** 0.05)
        pymax = ymax * ((ymax/ymin) ** 0.05)
    else:
        pymin = 0
        pymax = 100

    if bblocks:
        if bbmode == 'sc' or bbmode == 'both':
            err = np.max(sc_rate_err, axis=0)
            err = np.min([sc_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(sc_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None)
        elif bbmode == 'sum':
            err = np.max(sc_bg_rate_err, axis=0)
            err = np.min([sc_bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_bg_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(sc_bg_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None)
        if bbmode == 'both':
            err = np.max(bg_rate_err, axis=0)
            err = np.min([bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, bg_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(bg_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None,
                          alpha=alpha_bg)

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_mincounts_broken_bayes(hdulist, axs, log, mjdref, xflag,
                                   mincounts, color, obs_periods,
                                   short_time, stan_model, quantiles,
                                   time_rel=0, fexp_cut=0.15,
                                   alpha_bg=0.3, bblocks=False, bbp0=0.003,
                                   bbmode='both', yscale='linear'):
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
    ftime = hdulist[1].data.field('FRACTIME')[fexp_full > fexp_cut]
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
                if (time_mjd[istart-1] < period[0]
                        and time_mjd[istart] > period[0]):
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow - 1:
            tend = time[iend] + delt[iend]
        else:
            border_high = False
            for period in obs_periods:
                if (time_mjd[iend] < period[1]
                        and time_mjd[iend + 1] > period[1]):
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

    amp_dev_min = np.percentile(fit.stan_variables()['amp_dev_min'],
                                quantiles[1])
    amp_dev_min_low = np.percentile(fit.stan_variables()['amp_dev_min'],
                                    quantiles[0])
    amp_dev_min_up = np.percentile(fit.stan_variables()['amp_dev_min'],
                                   quantiles[2])
    amp_frac_min = np.percentile(fit.stan_variables()['amp_frac_min'],
                                 quantiles[1])
    amp_frac_min_low = np.percentile(
        fit.stan_variables()['amp_frac_min'], quantiles[0])
    amp_frac_min_up = np.percentile(
        fit.stan_variables()['amp_frac_min'], quantiles[2])
    logger_stan.warning(f'mincounts {mincounts} AmpVar(min)={amp_dev_min}'
                        f'+{amp_dev_min_up-amp_dev_min}'
                        f'-{amp_dev_min-amp_dev_min_low}\n')
    logger_stan.warning(f'mincounts {mincounts} AmpFrac(min)={amp_frac_min}'
                        f'+{amp_frac_min_up-amp_frac_min}'
                        f'-{amp_frac_min-amp_frac_min_low}\n')

    amp_dev_med = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[1])
    amp_dev_med_low = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[0])
    amp_dev_med_up = np.percentile(
        fit.stan_variables()['amp_dev_med'], quantiles[2])
    amp_frac_med = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[1])
    amp_frac_med_low = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[0])
    amp_frac_med_up = np.percentile(
        fit.stan_variables()['amp_frac_med'], quantiles[2])
    logger_stan.warning(f'mincounts {mincounts} AmpVar(med)={amp_dev_med}+'
                        f'{amp_dev_med_up-amp_dev_med}'
                        f'-{amp_dev_med-amp_dev_med_low}\n')
    logger_stan.warning(f'mincounts {mincounts} AmpFrac(med)={amp_frac_med}'
                        f'+{amp_frac_med_up-amp_frac_med}'
                        f'-{amp_frac_med-amp_frac_med_low}\n')

    yrate = np.array(sc_rate)
    yrate_lower = np.array(sc_rate_lower)
    yrate_upper = np.array(sc_rate_upper)
    yrate_e = np.sqrt(((yrate-yrate_lower)**2+(yrate_upper-yrate)**2)/2)

    i_max = np.argmax(yrate_lower)
    i_min = np.argmin(yrate_upper)

    ampl_max2 = yrate_lower[i_max] - yrate_upper[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    log.info(f'AMPL_MAX mincounts {mincounts}: {ampl_max}\n')
    log.info(f'Variability V mincounts {mincounts} = {variability}\n')
    log.info(f'AMPL_SIG mincounts {mincounts}: {ampl_sig}\n')
    log.info(f'AMPL_MAX conservative mincounts {mincounts}: {ampl_max2}\n')
    log.info(f'AMPL_SIG2 mincounts {mincounts}: {ampl_sig2}\n')

    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    # average rate:
    data_av = {}
    data_av['N'] = len(delt)
    data_av['M'] = 1
    data_av['dt'] = np.array([delt])
    data_av['sc'] = np.array([cnts], dtype=int)
    data_av['frac_exp'] = np.array([fexp])
    data_av['bg'] = np.array([back], dtype=int)
    data_av['bg_ratio'] = np.array([backrat])

    fit_av = model.sample(data=data_av, show_progress=False)
    sc_rate_av_lower = np.percentile(fit_av.stan_variables()['sc_rate'],
                                     quantiles[0], axis=0)[0]
    sc_rate_av = np.percentile(fit_av.stan_variables()['sc_rate'],
                               quantiles[1], axis=0)[0]
    sc_rate_av_upper = np.percentile(fit_av.stan_variables()['sc_rate'],
                                     quantiles[2], axis=0)[0]
    bg_rate_av_lower = np.percentile(fit_av.stan_variables()['bg_rate'],
                                     quantiles[0], axis=0)[0]
    bg_rate_av = np.percentile(fit_av.stan_variables()['bg_rate'],
                               quantiles[1], axis=0)[0]
    bg_rate_av_upper = np.percentile(fit_av.stan_variables()['bg_rate'],
                                     quantiles[2], axis=0)[0]
    logger_stan.warning(f'Average Source Rate={sc_rate_av}'
                        f'+{sc_rate_av_upper-sc_rate_av}'
                        f'-{sc_rate_av-sc_rate_av_lower} cts/s')
    logger_stan.warning(f'Average Background Rate={bg_rate_av}'
                        f'+{bg_rate_av_upper-bg_rate_av}'
                        f'-{bg_rate_av-bg_rate_av_lower} cts/s')
    tot_exp = ftime.sum()
    tot_fexp = (delt * fexp).sum()
    logger_stan.warning(f'Total Exposure={tot_exp}s')
    logger_stan.warning(f'Total Fractional Exposure={tot_fexp}s')

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
                    zorder=6)
        ax.errorbar(t, bg_rate, xerr=terr,
                    yerr=bg_rate_err,
                    linestyle='None', color=color, fmt='x',
                    zorder=5, alpha=alpha_bg)

    ymax = max([max(sc_rate_upper), max(bg_rate_upper)])
    if yscale == 'linear':
        ymin = 0
        pymin = ymin - (ymax-ymin)*0.05
        pymax = ymax + (ymax-ymin)*0.05
    elif yscale == 'log':
        if max(bg_rate_lower) > min(sc_rate_lower) * 1e-1 and alpha_bg > 0:
            ymin = min([min(sc_rate_lower), min(bg_rate_lower)])
        else:
            ymin = min(sc_rate_lower)
        pymin = ymin / ((ymax/ymin) ** 0.05)
        pymax = ymax * ((ymax/ymin) ** 0.05)
    else:
        pymin = 0
        pymax = 100

    if bblocks:
        if bbmode == 'sc' or bbmode == 'both':
            err = np.max(sc_rate_err, axis=0)
            err = np.min([sc_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(sc_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None)
        elif bbmode == 'sum':
            err = np.max(sc_bg_rate_err, axis=0)
            err = np.min([sc_bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, sc_bg_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(sc_bg_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None)
        if bbmode == 'both':
            err = np.max(bg_rate_err, axis=0)
            err = np.min([bg_rate, err], axis=0)
            t_blocked = bayesian_blocks(t, bg_rate, sigma=err,
                                        fitness='measures', p0=bbp0)
            for i in range(1, len(t_blocked) - 1):
                entry_tb = t_blocked[i]
                index_tb = len(t[t < entry_tb]) - 1
                t_blocked[i] = t[index_tb] + terr[index_tb]
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
            t_blocked_plot = t_blocked.copy()
            t_blocked_plot[0] = t[0] - terr[0]
            t_blocked_plot[-1] = t[-1] + terr[-1]
            for ax in axs:
                ax.stairs(bg_rate_bb, t_blocked_plot, color=color, zorder=2,
                          linestyle='--', baseline = None,
                          alpha=alpha_bg)

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel

def plot_hr_mincounts_broken_bayes(hdulist, axs, log, mjdref, xflag,
                                   mincounts_full, mincounts_ind, color,
                                   obs_periods, short_time, stan_model,
                                   quantiles, time_rel=0, fexp_cut=0.15,
                                   alpha_bg=0.3, bblocks=False, bbp0=0.003,
                                   bbmode='both', yscale='linear'):
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
    pymax = []
    xtime = []
    xtime_d = []
    fexp_full = hdulist[1].data.field('FRACEXP')[:,0]
    time = hdulist[1].data.field('TIME')[fexp_full > fexp_cut]
    time_mjd = time / 3600. / 24. + mjdref
    delt = hdulist[1].data.field('TIMEDEL')[fexp_full > fexp_cut]
    cnts0 = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[:,0][fexp_full > fexp_cut]
    cnts1 = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[:,1][fexp_full > fexp_cut]
    cnts2 = np.array(hdulist[1].data.field('COUNTS'),
                    dtype=int)[:,2][fexp_full > fexp_cut]
    fexp1 = hdulist[1].data.field('FRACEXP')[:,1][fexp_full > fexp_cut]
    fexp2 = hdulist[1].data.field('FRACEXP')[:,2][fexp_full > fexp_cut]
    back1 = np.array(hdulist[1].data.field('BACK_COUNTS'), dtype=int)[:,1][
        fexp_full > fexp_cut]
    back2 = np.array(hdulist[1].data.field('BACK_COUNTS'), dtype=int)[:,2][
        fexp_full > fexp_cut]
    backrat = hdulist[1].data.field('BACKRATIO')[fexp_full > fexp_cut]
    ftime = hdulist[1].data.field('FRACTIME')[fexp_full > fexp_cut]
    for i, entry in enumerate(backrat):
        if entry < 0.01:
            backrat[i] = 0.01

    dts = []

    scs1 = []
    fexps1 = []
    bgs1 = []
    scs2 = []
    fexps2 = []
    bgs2 = []

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
            if ((cnts0[istart:nrow].sum() < mincounts_full
                 or cnts1[istart:nrow].sum() < mincounts_ind
                 or cnts2[istart:nrow].sum() < mincounts_ind)
                 and not time_mjd[istart_old] < obs_periods[-1][0]):
                if istart_old > obs_periods[-1][0]:
                    istart = istart_old
                    del (xtime[-1], xtime_d[-1], dts[-1], bgrats[-1], scs1[-1],
                        fexps1[-1], bgs1[-1], scs2[-1], fexps2[-1], bgs2[-1])
        elif (cnts0[istart:i].sum() >= mincounts_full
              and cnts1[istart:i].sum() >= mincounts_ind
              and cnts2[istart:i].sum() >= mincounts_ind
              and time[i + 1] - time[i] > 3600.0):  # elif to avoid error
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
                    del (xtime[-1], xtime_d[-1], dts[-1], bgrats[-1], scs1[-1],
                         fexps1[-1], bgs1[-1],
                         scs2[-1], fexps2[-1], bgs2[-1])
            else:
                continue

        if istart == 0:
            tstart = time[0]
        else:
            border_low = False
            for period in obs_periods:
                if (time_mjd[istart-1] < period[0]
                        and time_mjd[istart] > period[0]):
                    tstart = time[istart] - delt[istart]
                    border_low = True
            if not border_low:
                tstart = (time[istart] + time[istart-1]) / 2
        if iend == nrow - 1:
            tend = time[iend] + delt[iend]
        else:
            border_high = False
            for period in obs_periods:
                if (time_mjd[iend] < period[1]
                        and time_mjd[iend + 1] > period[1]):
                    tend = time[iend] + delt[iend]
                    border_high = True
            if not border_high:
                tend = (time[iend] + time[iend+1]) / 2
        xtime.append((tend+tstart)/2)
        xtime_d.append((tend-tstart)/2)
        # data['N'] = iend-istart
        dts.append(delt[istart:iend+1].tolist())
        scs1.append(cnts1[istart:iend+1].tolist())
        fexps1.append(fexp1[istart:iend+1].tolist())
        bgs1.append(back1[istart:iend+1].tolist())
        scs2.append(cnts2[istart:iend+1].tolist())
        fexps2.append(fexp2[istart:iend+1].tolist())
        bgs2.append(back2[istart:iend+1].tolist())
        bgrats.append(backrat[istart:iend+1].tolist())
        istart_old = istart
        istart = i + 1

    N = 0
    for k in range(len(scs1)):
        if len(scs1[k]) > N:
            N = len(scs1[k])
    for k in range(len(scs1)):
        while len(scs1[k]) < N:
            dts[k].append(0)
            scs1[k].append(0)
            fexps1[k].append(0)
            bgs1[k].append(0)
            bgrats[k].append(0)
    for k in range(len(scs2)):
        while len(scs2[k]) < N:
            scs2[k].append(0)
            fexps2[k].append(0)
            bgs2[k].append(0)

    data = {}
    data['M'] = len(scs1)
    data['dt'] = np.array(dts)
    data['bg_ratio'] = np.array(bgrats)
    data['N'] = N
    data['sc1'] = np.array(scs1)
    data['frac_exp1'] = np.array(fexps1)
    data['bg1'] = np.array(bgs1)
    data['sc2'] = np.array(scs2)
    data['frac_exp2'] = np.array(fexps2)
    data['bg2'] = np.array(bgs2)

    # loading stan model
    model = CmdStanModel(stan_file=stan_model)
    fit = model.sample(data=data, show_progress=False)
    sc_rate_lower1 = np.percentile(fit.stan_variables()['sc_rate1'],
                                   quantiles[0], axis=0)
    sc_rate1 = np.percentile(fit.stan_variables()['sc_rate1'],
                             quantiles[1], axis=0)
    sc_rate_upper1 = np.percentile(fit.stan_variables()['sc_rate1'],
                                   quantiles[2], axis=0)
    bg_rate_lower1 = np.percentile(fit.stan_variables()['bg_rate1'],
                                   quantiles[0], axis=0)
    bg_rate1 = np.percentile(fit.stan_variables()['bg_rate1'],
                             quantiles[1], axis=0)
    bg_rate_upper1 = np.percentile(fit.stan_variables()['bg_rate1'],
                                   quantiles[2], axis=0)
    sc_bg_rate_lower1 = np.percentile(fit.stan_variables()['sc_bg_rate1'],
                                      quantiles[0], axis=0)
    sc_bg_rate1 = np.percentile(fit.stan_variables()['sc_bg_rate1'],
                                quantiles[1], axis=0)
    sc_bg_rate_upper1 = np.percentile(fit.stan_variables()['sc_bg_rate1'],
                                      quantiles[2], axis=0)
    sc_rate_lower2 = np.percentile(fit.stan_variables()['sc_rate2'],
                                   quantiles[0], axis=0)
    sc_rate2 = np.percentile(fit.stan_variables()['sc_rate2'],
                             quantiles[1], axis=0)
    sc_rate_upper2 = np.percentile(fit.stan_variables()['sc_rate2'],
                                   quantiles[2], axis=0)
    bg_rate_lower2 = np.percentile(fit.stan_variables()['bg_rate2'],
                                   quantiles[0], axis=0)
    bg_rate2 = np.percentile(fit.stan_variables()['bg_rate2'],
                             quantiles[1], axis=0)
    bg_rate_upper2 = np.percentile(fit.stan_variables()['bg_rate2'],
                                   quantiles[2], axis=0)
    sc_bg_rate_lower2 = np.percentile(fit.stan_variables()['sc_bg_rate2'],
                                      quantiles[0], axis=0)
    sc_bg_rate2 = np.percentile(fit.stan_variables()['sc_bg_rate2'],
                                quantiles[1], axis=0)
    sc_bg_rate_upper2 = np.percentile(fit.stan_variables()['sc_bg_rate2'],
                                      quantiles[2], axis=0)
    frac = np.percentile(fit.stan_variables()['frac'], quantiles[1], axis=0)
    frac_upper = np.percentile(fit.stan_variables()['frac'],
                               quantiles[2], axis=0)
    frac_lower = np.percentile(fit.stan_variables()['frac'],
                               quantiles[0], axis=0)
    
    t = xtime.copy()

    yrate = np.array(sc_rate1)
    yrate_lower = np.array(sc_rate_lower1)
    yrate_upper = np.array(sc_rate_upper1)
    yrate_e = np.sqrt(((yrate-yrate_lower)**2+(yrate_upper-yrate)**2)/2)

    i_max = np.argmax(yrate_lower)
    i_min = np.argmin(yrate_upper)

    ampl_max2 = yrate_lower[i_max] - yrate_upper[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    log.info(f'AMPL_MAX hr band 1 {mincounts_full}: {ampl_max}\n')
    log.info(f'Variability V mincounts band 1 {mincounts_full} = '
             f'{variability}\n')
    log.info(f'AMPL_SIG mincounts band 1 {mincounts_full}: {ampl_sig}\n')
    log.info(f'AMPL_MAX conservative mincounts band 1 {mincounts_full}: '
             f'{ampl_max2}\n')
    log.info(f'AMPL_SIG2 mincounts band 1 {mincounts_full}: {ampl_sig2}\n')
    
    yrate = np.array(sc_rate2)
    yrate_lower = np.array(sc_rate_lower2)
    yrate_upper = np.array(sc_rate_upper2)
    yrate_e = np.sqrt(((yrate-yrate_lower)**2+(yrate_upper-yrate)**2)/2)

    i_max = np.argmax(yrate_lower)
    i_min = np.argmin(yrate_upper)

    ampl_max2 = yrate_lower[i_max] - yrate_upper[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    log.info(f'AMPL_MAX hr band 2 {mincounts_ind}: {ampl_max}\n')
    log.info(f'Variability V mincounts band 2 {mincounts_ind} = '
             f'{variability}\n')
    log.info(f'AMPL_SIG mincounts band 2 {mincounts_ind}: {ampl_sig}\n')
    log.info(f'AMPL_MAX conservative mincounts band 2 {mincounts_ind}: '
             f'{ampl_max2}\n')
    log.info(f'AMPL_SIG2 mincounts band 2 {mincounts_ind}: {ampl_sig2}\n')


    if istart != nrow:
        raise Exception('Something went wrong in last bin.')

    xtime = np.array(xtime)
    xtime_d = np.array(xtime_d)
    mjd = xtime * (1. / 24. / 3600.) + mjdref
    mjd_d = xtime_d / 24. / 3600.
    sc_rate1 = np.array(sc_rate1)
    sc_rate_lower1 = np.array(sc_rate_lower1)
    sc_rate_upper1 = np.array(sc_rate_upper1)
    bg_rate1 = np.array(bg_rate1)
    bg_rate_lower1 = np.array(bg_rate_lower1)
    bg_rate_upper1 = np.array(bg_rate_upper1)
    sc_bg_rate1 = np.array(sc_bg_rate1)
    sc_bg_rate_lower1 = np.array(sc_bg_rate_lower1)
    sc_bg_rate_upper1 = np.array(sc_bg_rate_upper1)
    sc_rate_err1 = [sc_rate1 + (-sc_rate_lower1),
                    sc_rate_upper1 + (-sc_rate1)]
    bg_rate_err1 = [bg_rate1 + (-bg_rate_lower1),
                    bg_rate_upper1 + (-bg_rate1)]
    sc_bg_rate_err1 = [sc_bg_rate1 + (-sc_bg_rate_lower1),
                       sc_bg_rate_upper1 + (-sc_bg_rate1)]
    sc_rate_lower2 = np.array(sc_rate_lower2)
    sc_rate_upper2 = np.array(sc_rate_upper2)
    bg_rate2 = np.array(bg_rate2)
    bg_rate_lower2 = np.array(bg_rate_lower2)
    bg_rate_upper2 = np.array(bg_rate_upper2)
    sc_bg_rate2 = np.array(sc_bg_rate2)
    sc_bg_rate_lower2 = np.array(sc_bg_rate_lower2)
    sc_bg_rate_upper2 = np.array(sc_bg_rate_upper2)
    sc_rate_err2 = [sc_rate2 + (-sc_rate_lower2),
                    sc_rate_upper2 + (-sc_rate2)]
    bg_rate_err2 = [bg_rate2 + (-bg_rate_lower2),
                    bg_rate_upper2 + (-bg_rate2)]
    sc_bg_rate_err2 = [sc_bg_rate2 + (-sc_bg_rate_lower2),
                       sc_bg_rate_upper2 + (-sc_bg_rate2)]
    frac = np.array(frac)
    frac_lower = np.array(frac_lower)
    frac_upper = np.array(frac_upper)
    frac_err = [frac + (-frac_lower),
                frac_upper + (-frac)]
                    

    for i_ax, group in enumerate(axs):
        for i_en, ax in enumerate(group):
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
            if i_en == 0:
                pxmin.append(xmin - xm)
                pxmax.append(xmax + xm)

            if xflag == 1:
                t = xtime_short
                terr = xtime_d
            else:
                t = mjd_short
                terr = mjd_d
            if i_en == 0:
                ax.errorbar(t, sc_rate1, xerr=terr,
                            yerr=sc_rate_err1,
                            linestyle='None', color=color, fmt='o',
                            zorder=6)
                ax.errorbar(t, bg_rate1, xerr=terr,
                            yerr=bg_rate_err1,
                            linestyle='None', color=color, fmt='x',
                            zorder=5, alpha=alpha_bg)
            elif i_en == 1:
                ax.errorbar(t, sc_rate2, xerr=terr,
                            yerr=sc_rate_err2,
                            linestyle='None', color=color, fmt='o',
                            zorder=6)
                ax.errorbar(t, bg_rate2, xerr=terr,
                            yerr=bg_rate_err2,
                            linestyle='None', color=color, fmt='x',
                            zorder=5, alpha=alpha_bg)
            else:
                ax.errorbar(t, frac, xerr=terr,
                            yerr=frac_err,
                            linestyle='None', color=color, fmt='o',
                            zorder=6)

    if alpha_bg > 0:
        ymax1 = max([max(sc_rate_upper1), max(bg_rate_upper1)])
        ymax2 = max([max(sc_rate_upper2), max(bg_rate_upper2)])
    else:
        ymax1 = max(sc_rate_upper1)
        ymax2 = max(sc_rate_upper2)
    if yscale == 'linear':
        ymin = 0
        pymin = [ymin - (ymax1-ymin)*0.05,
                 ymin - (ymax2-ymin)*0.05,
                 -1]
        pymax = [ymax1 + (ymax1-ymin)*0.05,
                 ymax2 + (ymax2-ymin)*0.05,
                 1]
    elif yscale == 'log':
        if max(bg_rate_lower1) > min(sc_rate_lower1) * 1e-1 and alpha_bg > 0:
            ymin1 = min([min(sc_rate_lower1), min(bg_rate_lower1)])
        else:
            ymin1 = min(sc_rate_lower1)
        if max(bg_rate_lower2) > min(sc_rate_lower2) * 1e-1 and alpha_bg > 0:
            ymin2 = min([min(sc_rate_lower2), min(bg_rate_lower2)])
        else:
            ymin2 = min(sc_rate_lower2)
        pymin = [ymin1 / ((ymax1/ymin1) ** 0.05),
                 ymin2 / ((ymax2/ymin2) ** 0.05),
                 -1]
        pymax = [ymax1 * ((ymax1/ymin1) ** 0.05),
                 ymax2 * ((ymax2/ymin2) ** 0.05),
                 1]
    else:
        pymin = 0
        pymax = 100

    logger_stan.handlers = []

    return pxmin, pxmax, pymin, pymax, time_rel
