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


def plot_lc_eROday_broken_bayes_old(hdulist, axs, logfile, mjdref, xflag,
                                    color, obs_periods, short_time, stan_model,
                                    quantiles, time_rel=0, fexp_cut=0.15,
                                    alpha_bg=0.5):
    '''
    Lightcurve rebinned to eROdays with countrates optained with Bayesian fit
    assuming Poissionian distribution for counts and log
    '''
    logger = logging.getLogger('cmdstanpy')
    handler = logging.FileHandler(filename=logfile, mode='w')
    logger.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setLevel(logging.WARNING)
    logger.addHandler(handler)
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
    logging.basicConfig()

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_eROday_broken_bayes(hdulist, axs, logfile, mjdref, xflag,
                                color, obs_periods, short_time, stan_model,
                                quantiles, time_rel=0, fexp_cut=0.15,
                                alpha_bg=0.5):
    '''
    Lightcurve rebinned to eROdays with countrates optained with
    Bayesian fit assuming Poissionian distribution for counts and log;
    fit done for each bin simultaneously
    '''
    logger = logging.getLogger('cmdstanpy')
    handler = logging.FileHandler(filename=logfile, mode='w')
    logger.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setLevel(logging.WARNING)
    logger.addHandler(handler)
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
        # data['N'] = iend-istart
        dts.append(delt[istart:iend].tolist())
        scs.append(cnts[istart:iend].tolist())
        fexps.append(fexp[istart:iend].tolist())
        bgs.append(back[istart:iend].tolist())
        bgrats.append(backrat[istart:iend].tolist())
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
    logging.basicConfig()

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_mincounts_broken_bayes(hdulist, axs, logfile, mjdref, xflag,
                                   mincounts, color, obs_periods,
                                   short_time, stan_model, quantiles,
                                   time_rel=0, fexp_cut=0.15,
                                   alpha_bg=0.5):
    '''
    Lightcurve rebinned to eROdays with countrates optained with
    Bayesian fit assuming Poissionian distribution for counts and log;
    fit done for each bin simultaneously
    '''
    logger = logging.getLogger('cmdstanpy')
    handler = logging.FileHandler(filename=logfile, mode='w')
    logger.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setLevel(logging.WARNING)
    logger.addHandler(handler)
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
            iend = i + 1
            if (cnts[istart:nrow].sum() < mincounts and
                    not time_mjd[istart_old] < obs_periods[-1][0]):
                istart = istart_old
                del (xtime[-1], xtime_d[-1], dts[-1], scs[-1], fexps[-1],
                     bgs[-1], bgrats[-1])
        elif (cnts[istart:i].sum() >= mincounts and
              time[i + 1] - time[i] > 3600.0):  # elif to avoid error
            iend = i + 1
        else:  # TODO: case when at the end of an obs period
            period_last = False
            period_first = False
            for period in obs_periods:
                if time_mjd[i] < period[1] and time_mjd[i+1] > period[1]:
                    period_last = True
                    if time_mjd[istart_old] < period[0]:
                        period_first = True
            if period_last:
                iend = i + 1
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
        # data['N'] = iend-istart
        dts.append(delt[istart:iend].tolist())
        scs.append(cnts[istart:iend].tolist())
        fexps.append(fexp[istart:iend].tolist())
        bgs.append(back[istart:iend].tolist())
        bgrats.append(backrat[istart:iend].tolist())
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

    print(mjd)
    print(obs_periods)

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
    logging.basicConfig()

    return pxmin, pxmax, pymin, pymax, time_rel
