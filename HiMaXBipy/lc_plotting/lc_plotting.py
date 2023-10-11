#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Kaltenbrunner
"""
import numpy as np
from HiMaXBipy.io.package_data import round_to_1
import matplotlib.ticker as plticker
from matplotlib import ticker, rcParams
import matplotlib.pyplot as plt


def plot_lc_UL(hdulist, ax, log, mjdref, xflag, mincounts, color):
    time = hdulist[1].data.field('TIME')
    delt = hdulist[1].data.field('TIMEDEL')
    cnts = hdulist[1].data.field('COUNTS')
    fexp = hdulist[1].data.field('FRACEXP')
    ftim = hdulist[1].data.field('FRACTIME')
    farea = hdulist[1].data.field('FRACAREA')
    back = hdulist[1].data.field('BACK_COUNTS')
    backrat = hdulist[1].data.field('BACKRATIO')
    # delta = delt[0]
    # most representative value of background ratio
    backrat_med = np.median(backrat)

    nrow = len(time)
    log.info(f'Total number of time bins in table: {nrow}\n')
    start_time = time[0]
    end_time = time[-1]+delt[-1]
    time0 = time - start_time
    start_time0 = time0[0]
    end_time0 = time0[-1]+delt[-1]
    mjds = mjdref + start_time/24./3600.
    mjde = mjdref + end_time/24./3600.

    log.info(
        f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
    log.info(f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

    # calculate new rate by error propagation assuming there are enough
    # counts in each bin to assume Gaussian statistics usinf
    # infirmationon total counts background exposure and backrat
    uplimit = []

    # rate = []
    # rate_e = []
    xtime = []
    xtime_d = []
    yrate = []
    yrate_e = []
    start_bin = time0[0]
    end_bin = start_bin + delt[0]
    nbin = 1
    tcounts = 0
    tback = 0
    trate = 0.
    texp = 0.
    trate_e = 0.
    tcount_e = 0.
    tback_e = 0.
    ttim = 0.
    istart = 0
    iend = -1  # kald:for intended functionality
    counts = 0
    nrate = 0.
    nrate_e = 0.
    nexp = 0.
    bkg = 0.
    ncount_e = 0.
    nback_e = 0.
    narea = 0.
    ntime = 0.
    netcounts = -1.
    netcounts_e = -1.
    trate_ee = -1.

    for i in range(nrow):
        tmp = time0[i] - start_bin
        tcounts = tcounts + cnts[i]
        texp = texp + fexp[i]
        ttim = ttim + ftim[i]
        tback = tback + back[i]
        netcounts = (tcounts - tback*backrat_med)
        trate = (tcounts - tback*backrat_med)/(texp*delt[0])
        tcount_e = tcounts**0.5
        tback_e = tback**0.5
        trate_e = (tcount_e**2) + ((tback_e**2)*(backrat_med**2))
        netcounts_e = trate_e**0.5
        trate_ee = trate_e**0.5/(texp*delt[0])
        if tmp > 3600.0:
            # bin finished
            xtime.append(0.5*(time0[istart]+time0[iend]))
            xtime_d.append(0.5*(time0[iend]-time0[istart]))
            end_bin = time0[iend]
            counts = 0
            nrate = 0.
            nrate_e = 0.
            nexp = 0.
            bkg = 0.
            ncount_e = 0.
            nback_e = 0.
            narea = 0.
            ntime = 0.

            for j in range(istart, (iend+1)):
                nexp = nexp + fexp[j]
                counts = counts + cnts[j]
                bkg = bkg + back[j]
                nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                ncount_e = counts**0.5
                nback_e = bkg**0.5
                nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
                narea = narea + farea[j]
                ntime = ntime + ftim[j]
            if istart == iend+1:
                nexp = nexp + fexp[istart]
                counts = counts + cnts[istart]
                bkg = bkg + back[istart]
                nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                ncount_e = counts**0.5
                nback_e = bkg**0.5
                nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
            if counts < mincounts:
                uplimit.append(1)
            else:
                uplimit.append(0)
            log.debug(
                r'bin number, Start,end time of bin, duration of bin, ' +
                r'counts, rate, error,fracexp, fractime, fracarea, loop ' +
                f'start, loop end {nbin} {start_bin} {end_bin} ' +
                f'{end_bin-start_bin} {counts} {nrate} ' +
                f'{nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} ' +
                f'{istart} {iend}\n')
            # start next bin
            nbin += 1
            istart = i
            iend += 1
            start_bin = time0[istart]
        else:
            iend += 1
    end_bin = time0[iend]
    xtime.append(0.5*(time0[istart]+time0[iend]))
    xtime_d.append(0.5*(time0[iend]-time0[istart]))
    counts = 0
    nrate = 0.
    nrate_e = 0.
    nexp = 0.
    for j in range(istart, (iend+1)):
        nexp = nexp + fexp[j]
        counts = counts + cnts[j]
        bkg = bkg + back[j]
        nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
        ncount_e = counts**0.5
        nback_e = bkg**0.5
        nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
    if istart == iend+1:
        nexp = nexp + fexp[istart]
        counts = counts + cnts[istart]
        bkg = bkg + back[istart]
        nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
        ncount_e = counts**0.5
        nback_e = bkg**0.5
        nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))

    yrate.append(nrate)
    yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
    if counts < mincounts:
        uplimit.append(1)
    else:
        uplimit.append(0)

    log.debug(
        f'Start/end time of last bin {nbin} {start_bin} {end_bin} ' +
        f'{end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} ' +
        f'{nexp} {istart} {iend}\n')
    # start of first bin at 0:
    xtime = xtime-xtime[0]+xtime_d[0]
    log.debug(f'{xtime} {xtime_d}\n')
    mjd = np.array(xtime)/3600./24. + mjds
    mjd_d = np.array(xtime_d)/3600./24.

    log.info(f'Total counts: {tcounts} +/- {tcount_e}\n')
    log.info(f'Net counts {netcounts} +/- {netcounts_e}\n')
    log.info(f'binsize is {delt[0]}\n')
    log.info(f'Average rate: {trate} +/- {trate_ee}\n')

    log.info(f'Total exposure: {ttim}\n')
    log.info(f'Total fract. exposure: {texp}\n')

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    if xflag == 1:
        xmin = min(xtime)
        xmax = max(xtime)
    else:
        xmin = min(mjd)
        xmax = max(mjd)

    ymin = min(yrate + (-yrate_e))
    ymax = max(yrate + yrate_e * np.logical_not(uplimit))
    xm = (xmax-xmin)*0.05
    ym = (ymax-ymin)*0.05
    pxmin = xmin - xm
    pxmax = xmax + xm
    pymin = ymin - ym
    pymax = ymax + ym

    # Plot limits
    log.debug(f"Plot limits: {pxmin} {pxmax} {pymin} {pymax}\n")

    if xflag == 1:
        ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e, uplims=uplimit,
                    linestyle='None', color=color, fmt='o', zorder=1)
    else:
        ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e, uplims=uplimit,
                    linestyle='None', color=color, fmt='o', zorder=1)

    return pxmin, pxmax, pymin, pymax


def plot_lc_mincounts(hdulist, ax, log, mjdref, xflag, mincounts, color):
    time = hdulist[1].data.field('TIME')
    delt = hdulist[1].data.field('TIMEDEL')
    cnts = hdulist[1].data.field('COUNTS')
    fexp = hdulist[1].data.field('FRACEXP')
    ftim = hdulist[1].data.field('FRACTIME')
    farea = hdulist[1].data.field('FRACAREA')
    back = hdulist[1].data.field('BACK_COUNTS')
    backrat = hdulist[1].data.field('BACKRATIO')
    # delta = delt[0]
    backrat_med = np.median(backrat)

    nrow = len(time)
    log.info(f'Total number of time bins in table: {nrow}\n')
    start_time = time[0]
    end_time = time[-1]+delt[-1]
    time0 = time - start_time
    start_time0 = time0[0]
    end_time0 = time0[-1]+delt[-1]
    mjds = mjdref + start_time/24./3600.
    mjde = mjdref + end_time/24./3600.

    log.info(
        f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
    log.info(f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

    # calculate new rate by error propagation assuming there are enough
    # counts in each bin to assume Gaussian statistics usinf infirmation
    # on total counts background exposure and backrat
    # rate = []
    # rate_e = []
    xtime = []
    xtime_d = []
    yrate = []
    yrate_e = []
    start_bin = time0[0]
    end_bin = start_bin + delt[0]
    nbin = 1
    tcounts = 0
    tcounts_c = 0
    tback = 0
    trate = 0.
    texp = 0.
    trate_e = 0.
    tcount_e = 0.
    # tcount_e_c = 0.
    tback_e = 0.
    ttim = 0.
    istart = 0
    istart_tmp = 0
    iend = -1  # kald: for intended functionality
    ncount_e = 0.
    nback_e = 0.
    narea = 0.
    ntime = 0.
    nrat = []
    netcounts = -1
    netcounts_e = -1
    trate_ee = -1
    # mean_count = 0
    # sn = 2.0
    # xtime_tmp = 0.0
    # xtime_d_start_tmp = 0.0
    log.info(f'Minimum number of counts is {mincounts}\n')
    for i in range(nrow):
        # tmp = time0[i] - start_bin
        tcounts = tcounts + cnts[i]
        tcounts_c = tcounts_c + cnts[i]
        texp = texp + fexp[i]
        ttim = ttim + ftim[i]
        tback = tback + back[i]
        netcounts = (tcounts - tback*backrat_med)
        trate = (tcounts - tback*backrat_med)/(texp*delt[0])
        tcount_e = tcounts**0.5
        # tcount_e_c = tcounts_c**0.5
        tback_e = tback**0.5
        # if tcounts_c > 0:
        #     sn = tcounts_c/tcount_e_c
        # else:
        #     sn = 0.0
        trate_e = (tcount_e**2) + ((tback_e**2)*(backrat_med**2))
        netcounts_e = trate_e**0.5
        trate_ee = trate_e**0.5/(texp*delt[0])
        if tcounts_c > mincounts:
            tcounts_c = 0
            # tcount_e_c = 0.0
            if iend < nrow - 2:
                iend += 1  # kald: in this case one more bin needs to be
            # included compared to the usual lightcurves, else <10
            # countsfor all newly defined bins
            # without the +1 before but this way sections are directly
            # connected #kald
            xtime.append(0.5*(time0[istart]+time0[iend+1]))
            xtime_d.append(0.5*(time0[iend+1]-time0[istart]))
            # xtime_tmp = 0.5*(time0[istart]+time0[iend])
            # xtime_d_start_tmp = time0[istart]
            end_bin = time0[iend]
            counts = 0
            nrate = 0.
            nrate_e = 0.
            nexp = 0.
            bkg = 0.
            ncount_e = 0.
            nback_e = 0.
            narea = 0.
            ntime = 0.
            nrat = []
            # nrat_med = 0.
            for j in range(istart, iend + 1):
                nexp = nexp + fexp[j]
                counts = counts + cnts[j]
                bkg = bkg + back[j]
                nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                ncount_e = counts**0.5
                nback_e = bkg**0.5
                nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
                narea = narea + farea[j]
                ntime = ntime + ftim[j]
                nrat.append(backrat[j])
            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
            # nrat_med = np.median(nrat)

            log.debug(
                r'bin number, Start,end time of bin, duration of bin, ' +
                r'counts, rate, error,fracexp, fractime, fracarea, loop ' +
                f'start, loop end {nbin} {start_bin} {end_bin} ' +
                f'{end_bin-start_bin} {counts} {nrate} ' +
                f'{nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} ' +
                f'{istart} {iend}\n')
            # start next bin
            nbin += 1
            istart_tmp = istart
            if i < nrow - 1:
                istart = i + 1  # was only i #kald
                start_bin = time0[istart]
        else:
            iend += 1
    end_bin = time0[iend]
    # ycount_mean = 0
    # yerror_mean = 0.0
    # xtime_mean = 0.0
    counts = 0
    nrate = 0.
    nrate_e = 0.
    nexp = 0.
    bkg = 0.  # added #kald
    for j in range(istart, (iend+1)):
        nexp = nexp + fexp[j]
        counts = counts + cnts[j]
        bkg = bkg + back[j]
        nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
        ncount_e = counts**0.5
        nback_e = bkg**0.5
        nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
    if istart == iend+1:
        nexp = nexp + fexp[istart]
        counts = counts + cnts[istart]
        bkg = bkg + back[istart]
        nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
        ncount_e = counts**0.5
        nback_e = bkg**0.5
        nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
    log.debug(
        f'Start/end time of last bin {nbin} {start_bin} {end_bin} ' +
        f'{end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} ' +
        f'{nexp} {istart} {iend}\n')
    if counts > mincounts:  # this was tcounts_c which should be 0 all
        # the time at this point #kald
        log.debug(f'Condition satisfied with {counts} {bkg}\n')

        yrate.append(nrate)
        yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
        # see above comment about connected #kald
        xtime.append(0.5*(time0[istart]+time0[iend]))
        xtime_d.append(0.5*(time0[iend]-time0[istart]))
    else:
        log.debug('averaging\n')
        for j in range(istart_tmp, iend + 1):
            nexp = nexp + fexp[j]
            counts = counts + cnts[j]
            bkg = bkg + back[j]
            nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
            ncount_e = counts**0.5
            nback_e = bkg**0.5
            nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
            narea = narea + farea[j]
            ntime = ntime + ftim[j]
            nrat.append(backrat[j])
        if len(yrate) > 0:
            yrate[-1] = nrate
            yrate_e[-1] = nrate_e**0.5/(nexp*delt[0])

            #xtime_mean = (0.5*(time0[istart]+time0[iend]))
            xtime[-1] = (time0[istart_tmp] + time0[iend]) * \
                0.5  # changed to used istart_tmp #kald
            # taken out xerrer as a variable #kald
            xtime_d[-1] = 0.5*(time0[iend]-time0[istart_tmp])
        else:
            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))

            #xtime_mean = (0.5*(time0[istart]+time0[iend]))
            xtime.append((time0[istart_tmp] + time0[iend]) *
                         0.5)  # changed to used istart_tmp #kald
            # taken out xerrer as a variable #kald
            xtime_d.append(0.5*(time0[iend]-time0[istart_tmp]))

    # start of first bin at 0:
    mjd = np.array(xtime) / 3600. / 24. + start_time / \
        3600. / 24. + mjdref  # adjusted here #kald
    xtime = xtime-xtime[0]+xtime_d[0]
    mjd_d = np.array(xtime_d)/3600./24.

    log.info(f'Total counts: {tcounts} +/- {tcount_e}\n')
    log.info(f'Net counts {netcounts} +/- {netcounts_e}\n')
    log.info(f'binsize is {delt[0]}\n')
    log.info(f'Average rate: {trate} +/- {trate_ee}\n')

    log.info(f'Total exposure: {ttim}\n')
    log.info(f'Total fract. exposure: {texp}\n')

    if xflag == 1:
        xmin = xtime[0] - xtime_d[0]
        xmax = xtime[-1] + xtime_d[-1]
    else:
        xmin = mjd[0] - mjd_d[0]
        xmax = mjd[-1] + mjd_d[-1]

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    ymin = min(yrate + (-yrate_e))
    ymax = max(yrate + yrate_e)
    xm = (xmax-xmin)*0.05
    ym = (ymax-ymin)*0.05
    pxmin = xmin - xm
    pxmax = xmax + xm
    pymin = ymin - ym
    pymax = ymax + ym

    # Plot limits
    log.debug(f"Plot limits: {pxmin} {pxmax} {pymin} {pymax}\n")

    if xflag == 1:
        ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                    linestyle='None', color=color, fmt='o', zorder=2)
    else:
        ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                    linestyle='None', color=color, fmt='o', zorder=2)

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    i_max = np.argmax(yrate + (-yrate_e))
    i_min = np.argmax(yrate_e + (-yrate))

    ampl_max2 = yrate[i_max] - yrate_e[i_max] - yrate[i_min] - yrate_e[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    log.info(f'AMPL_MAX: {ampl_max}\n')
    log.info(f'Variability V = {variability}\n')
    log.info(f'AMPL_SIG: {ampl_sig}\n')
    log.info(f'AMPL_MAX conservative: {ampl_max2}\n')
    log.info(f'AMPL_SIG2: {ampl_sig2}\n')

    return pxmin, pxmax, pymin, pymax


def format_axis(ax, pxmin, pxmax, pymin, pymax, ticknumber_x, ticknumber_y):
    ax.set_xlim(pxmin, pxmax)
    ax.set_ylim(pymin, pymax)
    # this locator puts ticks at regular intervals
    loc = plticker.MultipleLocator(base=1.0)
    ax.yaxis.set_major_locator(loc)
    x_formatter = plticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)
    ax.tick_params(axis='x', which='major', direction='in',
                   top='on',   pad=5, length=5)  # , labelsize=10)
    ax.tick_params(axis='x', which='minor', direction='in',
                   top='on',   length=3)  # , labelsize=0)
    ax.tick_params(axis='y', which='major', direction='in',
                   right='on', length=5)  # , labelsize=10)
    ax.tick_params(axis='y', which='minor', direction='in',
                   right='on', length=3)  # , labelsize=0)

    tick_size_y = round_to_1((pymax - pymin) / ticknumber_y)
    yticks = []
    for i in range(-int(ticknumber_y), int(ticknumber_y)):
        if i * tick_size_y > pymin and i * tick_size_y < pymax:
            yticks.append(i * tick_size_y)
    ax.set_yticks(yticks)

    tick_size_x = round_to_1((pxmax - pxmin) / ticknumber_x)
    xticks = []
    centre_x = np.round((pxmax + pxmin) / 2., -int(np.log10(tick_size_x)))
    for i in range(-int(ticknumber_x), int(ticknumber_x)):
        if (i * tick_size_x + centre_x > pxmin and
                i * tick_size_x + centre_x < pxmax):
            xticks.append(i * tick_size_x + centre_x)
    ax.set_xticks(xticks)


def get_boundaries(in_bounds):
    out_bounds = np.zeros(4)
    for i in range(4):
        if i % 2 == 0:
            out_bounds[i] = np.min(np.array(in_bounds), axis=0)[i]
        else:
            out_bounds[i] = np.max(np.array(in_bounds), axis=0)[i]
    return out_bounds[0], out_bounds[1], out_bounds[2], out_bounds[3]


def get_boundaries_broken(in_bounds):
    out_bounds = []
    out_bounds.append(np.min(in_bounds, axis=0)[0])
    out_bounds.append(np.max(in_bounds, axis=0)[1])
    out_bounds.append(min(np.min(in_bounds, axis=0)[2]))
    out_bounds.append(max(np.max(in_bounds, axis=0)[3]))
    return out_bounds[0], out_bounds[1], out_bounds[2], out_bounds[3]


def plot_lc_mincounts_broken_new(hdulist, axs, log, mjdref, xflag,
                                 mincounts, color, obs_periods, short_time,
                                 time_rel=0):
    pxmax = []
    pxmin = []
    pymax = []
    pymin = []
    for i_ax, ax in enumerate(axs):
        time_full = hdulist[1].data.field('TIME')
        time = hdulist[1].data.field('TIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt = hdulist[1].data.field('TIMEDEL')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts = hdulist[1].data.field('COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp = hdulist[1].data.field('FRACEXP')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim = hdulist[1].data.field('FRACTIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea = hdulist[1].data.field('FRACAREA')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        back = hdulist[1].data.field('BACK_COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat = hdulist[1].data.field('BACKRATIO')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        backrat_med = np.median(backrat)

        nrow = len(time)
        log.info(f'Total number of time bins in table: {nrow}\n')
        start_time = time[0]
        end_time = time[-1]+delt[-1]
        time0 = time - start_time
        start_time0 = time0[0]
        end_time0 = time0[-1]+delt[-1]
        mjds = mjdref + start_time/24./3600.
        mjde = mjdref + end_time/24./3600.

        log.info(
            f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
        log.info(
            f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

        # calculate new rate by error propagation assuming there are
        # enough counts in each bin to assume Gaussian statistics using
        # infirmation on total counts background exposure and backrat
        # rate = []
        # rate_e = []
        xtime = []
        xtime_d = []
        yrate = []
        yrate_e = []
        start_bin = time0[0]
        end_bin = start_bin + delt[0]
        nbin = 1
        tcounts = 0
        tcounts_c = 0
        tback = 0
        trate = 0.
        texp = 0.
        trate_e = 0.
        tcount_e = 0.
        # tcount_e_c = 0.
        tback_e = 0.
        ttim = 0.
        istart = 0
        istart_tmp = 0
        iend = -1  # kald: for intended functionality
        ncount_e = 0.
        nback_e = 0.
        narea = 0.
        ntime = 0.
        nrat = []
        netcounts = -1.
        netcounts_e = -1.
        trate_ee = -1.
        # mean_count = 0
        # sn = 2.0
        # xtime_tmp = 0.0
        # xtime_d_start_tmp = 0.0
        log.info(f'Minimum number of counts is {mincounts}\n')
        for i in range(nrow):
            # tmp = time0[i] - start_bin
            tcounts = tcounts + cnts[i]
            tcounts_c = tcounts_c + cnts[i]
            texp = texp + fexp[i]
            ttim = ttim + ftim[i]
            tback = tback + back[i]
            netcounts = (tcounts - tback*backrat_med)
            trate = (tcounts - tback*backrat_med)/(texp*delt[0])
            tcount_e = tcounts**0.5
            # tcount_e_c = tcounts_c**0.5
            tback_e = tback**0.5
            # if tcounts_c > 0:
            #     sn = tcounts_c/tcount_e_c
            # else:
            #     sn = 0.0
            trate_e = (tcount_e**2) + ((tback_e**2)*(backrat_med**2))
            netcounts_e = trate_e**0.5
            trate_ee = trate_e**0.5/(texp*delt[0])
            if tcounts_c > mincounts:
                tcounts_c = 0
                # tcount_e_c = 0.0
                if iend < nrow - 2:
                    iend += 1  # kald: in this case one more bin needs to be
                # included compared to the usual lightcurves, else <10
                # countsfor all newly defined bins
                # without the +1 before but this way sections are
                # directly connected #kald
                xtime.append(0.5*(time0[istart]+time0[iend+1]))
                xtime_d.append(0.5*(time0[iend+1]-time0[istart]))
                # xtime_tmp = 0.5*(time0[istart]+time0[iend])
                # xtime_d_start_tmp = time0[istart]
                end_bin = time0[iend]
                counts = 0
                nrate = 0.
                nrate_e = 0.
                nexp = 0.
                bkg = 0.
                ncount_e = 0.
                nback_e = 0.
                narea = 0.
                ntime = 0.
                nrat = []
                # nrat_med = 0.
                for j in range(istart, iend + 1):
                    nexp = nexp + fexp[j]
                    counts = counts + cnts[j]
                    bkg = bkg + back[j]
                    nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                    ncount_e = counts**0.5
                    nback_e = bkg**0.5
                    nrate_e = (ncount_e**2) + ((nback_e**2) *
                                               (backrat_med**2))
                    narea = narea + farea[j]
                    ntime = ntime + ftim[j]
                    nrat.append(backrat[j])
                yrate.append(nrate)
                yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
                # nrat_med = np.median(nrat)

                log.debug(
                    r'bin number, Start,end time of bin, duration of bin, ' +
                    r'counts , rate, error, fracexp, fractime, fracarea, ' +
                    f'loop start, loop end {nbin} {start_bin} {end_bin} ' +
                    f'{end_bin-start_bin} {counts} {nrate} ' +
                    f'{nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} ' +
                    f'{istart} {iend}\n')
                # start next bin
                nbin += 1
                istart_tmp = istart
                if i < nrow - 1:
                    istart = i + 1  # was only i #kald
                    start_bin = time0[istart]
            else:
                iend += 1
        end_bin = time0[iend]
        # ycount_mean = 0
        # yerror_mean = 0.0
        # xtime_mean = 0.0
        counts = 0
        nrate = 0.
        nrate_e = 0.
        nexp = 0.
        bkg = 0.  # added #kald
        for j in range(istart, (iend+1)):
            nexp = nexp + fexp[j]
            counts = counts + cnts[j]
            bkg = bkg + back[j]
            nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
            ncount_e = counts**0.5
            nback_e = bkg**0.5
            nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
        if istart == iend+1:
            nexp = nexp + fexp[istart]
            counts = counts + cnts[istart]
            bkg = bkg + back[istart]
            nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
            ncount_e = counts**0.5
            nback_e = bkg**0.5
            nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
        log.info(
            f'Start/end time of last bin {nbin} {start_bin} {end_bin} ' +
            f'{end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} ' +
            f'{nexp} {istart} {iend}\n')
        if counts > mincounts:  # this was tcounts_c which should be 0
            # all the time at this point #kald
            log.debug(f'Condition satisfied with {counts} {bkg}\n')

            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
            # see above comment about connected #kald
            xtime.append(0.5*(time0[istart]+time0[iend]))
            xtime_d.append(0.5*(time0[iend]-time0[istart]))
        else:
            log.debug('averaging\n')
            for j in range(istart_tmp, iend + 1):
                nexp = nexp + fexp[j]
                counts = counts + cnts[j]
                bkg = bkg + back[j]
                nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                ncount_e = counts**0.5
                nback_e = bkg**0.5
                nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
                narea = narea + farea[j]
                ntime = ntime + ftim[j]
                nrat.append(backrat[j])
            if len(yrate) > 0:
                yrate[-1] = nrate
                yrate_e[-1] = nrate_e**0.5/(nexp*delt[0])

                #xtime_mean = (0.5*(time0[istart]+time0[iend]))
                xtime[-1] = (time0[istart_tmp] + time0[iend]) * \
                    0.5  # changed to used istart_tmp #kald
                # taken out xerrer as a variable #kald
                xtime_d[-1] = 0.5*(time0[iend]-time0[istart_tmp])
            else:
                yrate.append(nrate)
                yrate_e.append(nrate_e**0.5/(nexp*delt[0]))

                #xtime_mean = (0.5*(time0[istart]+time0[iend]))
                xtime.append((time0[istart_tmp] + time0[iend]) *
                             0.5)  # changed to used istart_tmp #kald
                # taken out xerrer as a variable #kald
                xtime_d.append(0.5*(time0[iend]-time0[istart_tmp]))

        # start of first bin at 0:
        mjd = np.array(xtime) / 3600. / 24. + start_time / \
            3600. / 24. + mjdref  # adjusted here #kald
        xtime = xtime-xtime[0]+xtime_d[0]
        mjd_d = np.array(xtime_d)/3600./24.

        log.info(f'Total counts: {tcounts} +/- {tcount_e}\n')
        log.info(f'Net counts {netcounts} +/- {netcounts_e}\n')
        log.info(f'binsize is {delt[0]}\n')
        log.info(f'Average rate: {trate} +/- {trate_ee}\n')

        log.info(f'Total exposure: {ttim}\n')
        log.info(f'Total fract. exposure: {texp}\n')

        if xflag == 1:
            xmin = xtime[0] - xtime_d[0]
            xmax = xtime[-1] + xtime_d[-1]
        else:
            xmin = mjd[0] - mjd_d[0]
            xmax = mjd[-1] + mjd_d[-1]

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(np.round(xmin))
            xtime = xtime - time_rel
            mjd = mjd - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel

        yrate = np.array(yrate)
        yrate_e = np.array(yrate_e)

        ymin = min(yrate + (-yrate_e))
        ymax = max(yrate + yrate_e)
        xm = (xmax-xmin)*0.05
        ym = (ymax-ymin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)
        pymin.append(ymin - ym)
        pymax.append(ymax + ym)

        # Plot limits
        log.debug(
            f"Plot limits: {pxmin[i_ax]} {pxmax[i_ax]} {pymin[i_ax]} " +
            f"{pymax[i_ax]}\n")

        if xflag == 1:
            ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                        linestyle='None', color=color, fmt='o', zorder=2)
        else:
            ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                        linestyle='None', color=color, fmt='o', zorder=2)

        yrate = np.array(yrate)
        yrate_e = np.array(yrate_e)

        i_max = np.argmax(yrate + (-yrate_e))
        i_min = np.argmax(yrate_e + (-yrate))

        ampl_max2 = yrate[i_max] - yrate_e[i_max] - \
            yrate[i_min] - yrate_e[i_min]
        ampl_max = yrate[i_max] - yrate[i_min]
        if yrate[i_min] > 0:
            variability = yrate[i_max] / yrate[i_min]
        else:
            variability = -1
        ampl_sig = ampl_max / \
            np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
        ampl_sig2 = ampl_max2 / \
            np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

        log.info(f'AMPL_MAX: {ampl_max}\n')
        log.info(f'Variability V = {variability}\n')
        log.info(f'AMPL_SIG: {ampl_sig}\n')
        log.info(f'AMPL_MAX conservative: {ampl_max2}\n')
        log.info(f'AMPL_SIG2: {ampl_sig2}\n')

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_UL_broken_new(hdulist, axs, log, mjdref, xflag, mincounts,
                          color, obs_periods, short_time, time_rel=0):
    pxmax = []
    pxmin = []
    pymax = []
    pymin = []
    for i_ax, ax in enumerate(axs):
        time_full = hdulist[1].data.field('TIME')
        time = hdulist[1].data.field('TIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt = hdulist[1].data.field('TIMEDEL')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts = hdulist[1].data.field('COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp = hdulist[1].data.field('FRACEXP')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim = hdulist[1].data.field('FRACTIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea = hdulist[1].data.field('FRACAREA')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        back = hdulist[1].data.field('BACK_COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat = hdulist[1].data.field('BACKRATIO')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        # most representative value of background ratio
        backrat_med = np.median(backrat)

        nrow = len(time)
        log.info(f'Total number of time bins in table: {nrow}\n')
        start_time = time[0]
        end_time = time[-1]+delt[-1]
        time0 = time - start_time
        start_time0 = time0[0]
        end_time0 = time0[-1]+delt[-1]
        mjds = mjdref + start_time/24./3600.
        mjde = mjdref + end_time/24./3600.

        log.info(
            f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
        log.info(
            f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

        # calculate new rate by error propagation assuming there are
        # enough counts in each bin to assume Gaussian statistics using
        # infirmation on total counts background exposure and backrat
        uplimit = []

        # rate = []
        # rate_e = []
        xtime = []
        xtime_d = []
        yrate = []
        yrate_e = []
        start_bin = time0[0]
        end_bin = start_bin + delt[0]
        nbin = 1
        tcounts = 0
        tback = 0
        trate = 0.
        texp = 0.
        trate_e = 0.
        tcount_e = 0.
        tback_e = 0.
        ttim = 0.
        istart = 0
        iend = -1  # kald:for intended functionality
        counts = 0
        nrate = 0.
        nrate_e = 0.
        nexp = 0.
        bkg = 0.
        ncount_e = 0.
        nback_e = 0.
        narea = 0.
        ntime = 0.
        netcounts = -1.
        netcounts_e = -1.
        trate_ee = -1.

        for i in range(nrow):
            tmp = time0[i] - start_bin
            tcounts = tcounts + cnts[i]
            texp = texp + fexp[i]
            ttim = ttim + ftim[i]
            tback = tback + back[i]
            netcounts = (tcounts - tback*backrat_med)
            trate = (tcounts - tback*backrat_med)/(texp*delt[0])
            tcount_e = tcounts**0.5
            tback_e = tback**0.5
            trate_e = (tcount_e**2) + ((tback_e**2)*(backrat_med**2))
            netcounts_e = trate_e**0.5
            trate_ee = trate_e**0.5/(texp*delt[0])
            if tmp > 3600.0:
                # bin finished
                xtime.append(0.5*(time0[istart]+time0[iend]))
                xtime_d.append(0.5*(time0[iend]-time0[istart]))
                end_bin = time0[iend]
                counts = 0
                nrate = 0.
                nrate_e = 0.
                nexp = 0.
                bkg = 0.
                ncount_e = 0.
                nback_e = 0.
                narea = 0.
                ntime = 0.

                for j in range(istart, (iend+1)):
                    nexp = nexp + fexp[j]
                    counts = counts + cnts[j]
                    bkg = bkg + back[j]
                    nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                    ncount_e = counts**0.5
                    nback_e = bkg**0.5
                    nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
                    narea = narea + farea[j]
                    ntime = ntime + ftim[j]
                if istart == iend+1:
                    nexp = nexp + fexp[istart]
                    counts = counts + cnts[istart]
                    bkg = bkg + back[istart]
                    nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
                    ncount_e = counts**0.5
                    nback_e = bkg**0.5
                    nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
                yrate.append(nrate)
                yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
                if counts < mincounts:
                    uplimit.append(1)
                else:
                    uplimit.append(0)
                log.debug(
                    r'bin number, Start,end time of bin, duration of bin, ' +
                    r'counts , rate, error,fracexp, fractime, fracarea, ' +
                    f'loop start, loop end {nbin} {start_bin} {end_bin} ' +
                    f'{end_bin-start_bin} {counts} {nrate} ' +
                    f'{nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} ' +
                    f'{istart} {iend}\n')
                # start next bin
                nbin += 1
                iend += 1
                if i <= nrow - 1:
                    istart = i
                    start_bin = time0[istart]
            else:
                iend += 1
        end_bin = time0[iend]
        xtime.append(0.5*(time0[istart]+time0[iend]))
        xtime_d.append(0.5*(time0[iend]-time0[istart]))
        counts = 0
        nrate = 0.
        nrate_e = 0.
        nexp = 0.
        for j in range(istart, (iend+1)):
            nexp = nexp + fexp[j]
            counts = counts + cnts[j]
            bkg = bkg + back[j]
            nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
            ncount_e = counts**0.5
            nback_e = bkg**0.5
            nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))
        if istart == iend+1:
            nexp = nexp + fexp[istart]
            counts = counts + cnts[istart]
            bkg = bkg + back[istart]
            nrate = (counts - bkg*backrat_med)/(nexp*delt[0])
            ncount_e = counts**0.5
            nback_e = bkg**0.5
            nrate_e = (ncount_e**2) + ((nback_e**2)*(backrat_med**2))

        yrate.append(nrate)
        yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
        if counts < mincounts:
            uplimit.append(1)
        else:
            uplimit.append(0)

        log.debug(
            f'Start/end time of last bin {nbin} {start_bin} {end_bin} ' +
            f'{end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} ' +
            f'{nexp} {istart} {iend}\n')
        # start of first bin at 0:
        xtime = xtime-xtime[0]+xtime_d[0]
        log.debug(f'{xtime} {xtime_d}\n')
        mjd = np.array(xtime)/3600./24. + mjds
        mjd_d = np.array(xtime_d)/3600./24.

        log.info(f'Total counts: {tcounts} +/- {tcount_e}\n')
        log.info(f'Net counts {netcounts} +/- {netcounts_e}\n')
        log.info(f'binsize is {delt[0]}\n')
        log.info(f'Average rate: {trate} +/- {trate_ee}\n')

        log.info(f'Total exposure: {ttim}\n')
        log.info(f'Total fract. exposure: {texp}\n')

        yrate = np.array(yrate)
        yrate_e = np.array(yrate_e)

        if xflag == 1:
            xmin = min(xtime)
            xmax = max(xtime)
        else:
            xmin = min(mjd)
            xmax = max(mjd)

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(np.round(xmin))
            xtime = xtime - time_rel
            mjd = mjd - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel

        ymin = min(yrate + (-yrate_e))
        ymax = max(yrate + yrate_e * np.logical_not(uplimit))
        xm = (xmax-xmin)*0.05
        ym = (ymax-ymin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)
        pymin.append(ymin - ym)
        pymax.append(ymax + ym)

        # Plot limits
        log.debug(
            f"Plot limits: {pxmin[i_ax]} {pxmax[i_ax]} {pymin[i_ax]} " +
            f"{pymax[i_ax]}\n")

        if xflag == 1:
            ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                        uplims=uplimit, linestyle='None', color=color, fmt='o',
                        zorder=1)
        else:
            ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                        uplims=uplimit, linestyle='None', color=color, fmt='o',
                        zorder=1)

    return pxmin, pxmax, pymin, pymax, time_rel


def format_axis_broken_new(fig, axs, pxmins, pxmaxs, pymin, pymax,
                           ticknumber_x, ticknumber_y, ncols, nrows, d, tilt,
                           diag_color, big_ax, yscale):
    # proportion of vertical to horizontal extent of the slanted line
    d = np.tan(tilt)
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none",
                  color=diag_color, mec=diag_color, mew=1, clip_on=False)

    start_x, end_x = 0, 0

    big_ax.tick_params(left=True, bottom=True,
                       right=True, top=True)
    loc = plticker.MultipleLocator(base=1.0)
    big_ax.yaxis.set_major_locator(loc)
    x_formatter = plticker.ScalarFormatter(useOffset=False)
    big_ax.xaxis.set_major_formatter(x_formatter)
    big_ax.tick_params(axis='x', which='major', direction='in',
                       top='on',   pad=9, length=0)  # , labelsize=10)
    big_ax.tick_params(axis='x', which='minor', direction='in',
                       top='on',   length=0)  # , labelsize=0)
    big_ax.tick_params(axis='y', which='major', direction='in',
                       right='on', length=0)  # , labelsize=10)
    big_ax.tick_params(axis='y', which='minor', direction='in',
                       right='on', length=0)

    for i_ax, ax in enumerate(axs):
        loc = plticker.MultipleLocator(base=1.0)
        ax.yaxis.set_major_locator(loc)
        x_formatter = plticker.ScalarFormatter(useOffset=False)
        ax.xaxis.set_major_formatter(x_formatter)
        ax.tick_params(axis='x', which='major', direction='in',
                       top='on',   pad=9, length=5)  # , labelsize=10)
        ax.tick_params(axis='x', which='minor', direction='in',
                       top='on',   length=3)  # , labelsize=0)
        ax.tick_params(axis='y', which='major', direction='in',
                       right='on', length=5)  # , labelsize=10)
        ax.tick_params(axis='y', which='minor', direction='in',
                       right='on', length=3)  # , labelsize=0)

        if i_ax == 0:
            ax.plot([1, 1], [0, 1], transform=ax.transAxes, **kwargs)
            ax.spines.right.set_visible(False)
            ax.tick_params(which='major', right=False, labelright=False)
            ax.tick_params(which='minor', right=False, labelright=False)
        elif i_ax == len(axs) - 1:
            ax.plot([0, 0], [0, 1], transform=ax.transAxes, **kwargs)
            ax.spines.left.set_visible(False)
            ax.tick_params(which='major', left=False, labelleft=False)
            ax.tick_params(which='minor', left=False, labelleft=False)
        else:
            ax.plot([0, 0, 1, 1], [0, 1, 0, 1],
                    transform=ax.transAxes, **kwargs)
            ax.spines.right.set_visible(False)
            ax.spines.left.set_visible(False)
            ax.tick_params(which='major', left=False, labelleft=False,
                           right=False, labelright=False)
            ax.tick_params(which='minor', left=False, labelleft=False,
                           right=False, labelright=False)

        if i_ax == 0 or i_ax == len(axs) - 1:
            if yscale == 'linear':
                tick_size_y = round_to_1((pymax - pymin) / ticknumber_y)
                yticks = []
                for j in range(-2*int(ticknumber_y), 2*int(ticknumber_y) + 1):
                    if j * tick_size_y > pymin and j * tick_size_y < pymax:
                        yticks.append(j * tick_size_y)
                ax.set_yticks(yticks)
                big_ax.set_yticks(yticks)
            elif yscale == 'log':
                yticks = []
                for power in range(int(np.floor(np.log10(pymin))),
                                   int(np.ceil(np.log10(pymax)))):
                    if (np.log10(pymax/pymin) >= 2 * ticknumber_y
                            and power % 2 == 1):
                        continue
                    yticks.append(10**power)
                    if np.log10(pymax/pymin) * 2 <= ticknumber_y:
                        yticks.append(2 * 10 ** power)
                    if np.log10(pymax/pymin) <= ticknumber_y:
                        yticks.append(5 * 10 ** power)
                ax.set_yticks(yticks)
                big_ax.set_yticks(yticks)

            # longest_y = ''
            # for entry in ax.get_yticks():
            #     entry = round_to_1(entry)
            #     if len(str(entry)) > len(longest_y):
            #         longest_y = str(entry)

        tick_size_x = np.round(
            (pxmaxs[i_ax] - pxmins[i_ax]) / ticknumber_x)
        if ticknumber_x % 2 == 0:
            shift_x = np.round(tick_size_x / 2)
        else:
            shift_x = 0

        xticks = []
        centre_x = np.round(
            (pxmaxs[i_ax] + pxmins[i_ax]) / 2.)
        for j in range(-int(ticknumber_x), int(ticknumber_x) + 1):
            if (j * tick_size_x + centre_x - shift_x > pxmins[i_ax] and
                    j * tick_size_x + centre_x - shift_x < pxmaxs[i_ax]):
                xticks.append(j * tick_size_x + centre_x - shift_x)
        if len(xticks) <= 1:
            xticks = [centre_x - np.round(tick_size_x / 2),
                      centre_x + np.round(tick_size_x / 2)]
        ax.set_xticks(xticks)
        if i_ax == 0:
            start_x = xticks[0]
        elif i_ax == len(axs) - 1:
            end_x = xticks[-1]

        ax.set_xbound(lower=pxmins[i_ax], upper=pxmaxs[i_ax])
        ax.set_ybound(lower=pymin, upper=pymax)
    big_ax.set_ybound(lower=pymin, upper=pymax)

    big_ax.set_xbound(lower=0, upper=1)
    big_ax.set_xticks([0, 1])

    big_ax.set_xticklabels([start_x, end_x], alpha=0)
    big_ax.yaxis.set_ticklabels(big_ax.yaxis.get_ticklabels(), alpha=0)


def plot_lc_mincounts_hr(hdulist_1, hdulist_2, axs, log, mjdref, xflag,
                         mincounts, color, obs_periods, short_time,
                         time_rel=0):
    return [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]


def plot_lc_UL_hr(hdulist_1, hdulist_2, axs, log, mjdref, xflag, mincounts,
                  color, obs_periods, short_time, time_rel=0):
    pxmax = []
    pxmin = []
    pymax = []
    pymin = []
    for i_ax, ax in enumerate(axs):
        time_full_1 = hdulist_1[1].data.field('TIME')
        time_1 = hdulist_1[1].data.field('TIME')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt_1 = hdulist_1[1].data.field('TIMEDEL')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts_1 = hdulist_1[1].data.field('COUNTS')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp_1 = hdulist_1[1].data.field('FRACEXP')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim_1 = hdulist_1[1].data.field('FRACTIME')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea_1 = hdulist_1[1].data.field('FRACAREA')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        back_1 = hdulist_1[1].data.field('BACK_COUNTS')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat_1 = hdulist_1[1].data.field('BACKRATIO')[(
            time_full_1/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_1/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        # most representative value of background ratio
        backrat_med_1 = np.median(backrat_1)

        nrow_1 = len(time_1)
        start_time_1 = time_1[0]
        end_time_1 = time_1[-1]+delt_1[-1]
        time0_1 = time_1 - start_time_1
        start_time0_1 = time0_1[0]
        end_time0_1 = time0_1[-1]+delt_1[-1]
        mjds_1 = mjdref + start_time_1/24./3600.
        mjde_1 = mjdref + end_time_1/24./3600.

        # calculate new rate by error propagation assuming there are
        # enough counts in each bin to assume Gaussian statistics using
        # infirmation on total counts background exposure and backrat
        uplimit_1 = []

        # rate = []
        # rate_e = []
        xtime_1 = []
        xtime_d_1 = []
        yrate_1 = []
        yrate_e_1 = []
        start_bin_1 = time0_1[0]
        end_bin_1 = start_bin_1 + delt_1[0]
        nbin_1 = 1
        tcounts_1 = 0
        tback_1 = 0
        trate_1 = 0.
        texp_1 = 0.
        trate_e_1 = 0.
        tcount_e_1 = 0.
        tback_e_1 = 0.
        ttim_1 = 0.
        istart_1 = 0
        iend_1 = -1  # kald:for intended functionality
        counts_1 = 0
        nrate_1 = 0.
        nrate_e_1 = 0.
        nexp_1 = 0.
        bkg_1 = 0.
        ncount_e_1 = 0.
        nback_e_1 = 0.
        narea_1 = 0.
        ntime_1 = 0.
        netcounts_1 = -1.
        netcounts_e_1 = -1.
        trate_ee_1 = -1.

        for i in range(nrow_1):
            tmp_1 = time0_1[i] - start_bin_1
            tcounts_1 = tcounts_1 + cnts_1[i]
            texp_1 = texp_1 + fexp_1[i]
            ttim_1 = ttim_1 + ftim_1[i]
            tback_1 = tback_1 + back_1[i]
            netcounts_1 = (tcounts_1 - tback_1*backrat_med_1)
            trate_1 = (tcounts_1 - tback_1*backrat_med_1)/(texp_1*delt_1[0])
            tcount_e_1 = tcounts_1**0.5
            tback_e_1 = tback_1**0.5
            trate_e_1 = (tcount_e_1**2) + ((tback_e_1**2)*(backrat_med_1**2))
            netcounts_e_1 = trate_e_1**0.5
            trate_ee_1 = trate_e_1**0.5/(texp_1*delt_1[0])
            if tmp_1 > 3600.0:
                # bin finished
                xtime_1.append(0.5*(time0_1[istart_1]+time0_1[iend_1]))
                xtime_d_1.append(0.5*(time0_1[iend_1]-time0_1[istart_1]))
                end_bin_1 = time0_1[iend_1]
                counts_1 = 0
                nrate_1 = 0.
                nrate_e_1 = 0.
                nexp_1 = 0.
                bkg_1 = 0.
                ncount_e_1 = 0.
                nback_e_1 = 0.
                narea_1 = 0.
                ntime_1 = 0.

                for j in range(istart_1, (iend_1+1)):
                    nexp_1 = nexp_1 + fexp_1[j]
                    counts_1 = counts_1 + cnts_1[j]
                    bkg_1 = bkg_1 + back_1[j]
                    nrate_1 = (counts_1 - bkg_1*backrat_med_1) / \
                        (nexp_1*delt_1[0])
                    ncount_e_1 = counts_1**0.5
                    nback_e_1 = bkg_1**0.5
                    nrate_e_1 = (ncount_e_1**2) + \
                        ((nback_e_1**2)*(backrat_med_1**2))
                    narea_1 = narea_1 + farea_1[j]
                    ntime_1 = ntime_1 + ftim_1[j]
                if istart_1 == iend_1+1:
                    nexp_1 = nexp_1 + fexp_1[istart_1]
                    counts_1 = counts_1 + cnts_1[istart_1]
                    bkg_1 = bkg_1 + back_1[istart_1]
                    nrate_1 = (counts_1 - bkg_1*backrat_med_1) / \
                        (nexp_1*delt_1[0])
                    ncount_e_1 = counts_1**0.5
                    nback_e_1 = bkg_1**0.5
                    nrate_e_1 = (ncount_e_1**2) + \
                        ((nback_e_1**2)*(backrat_med_1**2))
                yrate_1.append(nrate_1)
                yrate_e_1.append(nrate_e_1**0.5/(nexp_1*delt_1[0]))
                if counts_1 < mincounts:
                    uplimit_1.append(1)
                else:
                    uplimit_1.append(0)
                # start next bin
                nbin_1 += 1
                iend_1 += 1
                if i <= nrow_1 - 1:
                    istart_1 = i
                    start_bin_1 = time0_1[istart_1]
            else:
                iend_1 += 1
        end_bin_1 = time0_1[iend_1]
        xtime_1.append(0.5*(time0_1[istart_1]+time0_1[iend_1]))
        xtime_d_1.append(0.5*(time0_1[iend_1]-time0_1[istart_1]))
        counts_1 = 0
        nrate_1 = 0.
        nrate_e_1 = 0.
        nexp_1 = 0.
        for j in range(istart_1, (iend_1+1)):
            nexp_1 = nexp_1 + fexp_1[j]
            counts_1 = counts_1 + cnts_1[j]
            bkg_1 = bkg_1 + back_1[j]
            nrate_1 = (counts_1 - bkg_1*backrat_med_1)/(nexp_1*delt_1[0])
            ncount_e_1 = counts_1**0.5
            nback_e_1 = bkg_1**0.5
            nrate_e_1 = (ncount_e_1**2) + ((nback_e_1**2)*(backrat_med_1**2))
        if istart_1 == iend_1+1:
            nexp_1 = nexp_1 + fexp_1[istart_1]
            counts_1 = counts_1 + cnts_1[istart_1]
            bkg_1 = bkg_1 + back_1[istart_1]
            nrate_1 = (counts_1 - bkg_1*backrat_med_1)/(nexp_1*delt_1[0])
            ncount_e_1 = counts_1**0.5
            nback_e_1 = bkg_1**0.5
            nrate_e_1 = (ncount_e_1**2) + ((nback_e_1**2)*(backrat_med_1**2))

        yrate_1.append(nrate_1)
        yrate_e_1.append(nrate_e_1**0.5/(nexp_1*delt_1[0]))
        if counts_1 < mincounts:
            uplimit_1.append(1)
        else:
            uplimit_1.append(0)

        # start of first bin at 0:
        xtime_1 = xtime_1-xtime_1[0]+xtime_d_1[0]
        mjd_1 = np.array(xtime_1)/3600./24. + mjds_1
        mjd_d_1 = np.array(xtime_d_1)/3600./24.
        xtime_d_1 = np.array(xtime_d_1)

        yrate_1 = np.array(yrate_1)
        yrate_e_1 = np.array(yrate_e_1)

        if xflag == 1:
            xmin = min(xtime_1)
            xmax = max(xtime_1)
        else:
            xmin = min(mjd_1)
            xmax = max(mjd_1)

        if short_time:
            if i_ax == 0 and time_rel == 0:
                time_rel = int(np.round(xmin))
            xtime_1 = xtime_1 - time_rel
            mjd_1 = mjd_1 - time_rel
            xmin = xmin - time_rel
            xmax = xmax - time_rel

        time_full_2 = hdulist_2[1].data.field('TIME')
        time_2 = hdulist_2[1].data.field('TIME')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt_2 = hdulist_2[1].data.field('TIMEDEL')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts_2 = hdulist_2[1].data.field('COUNTS')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp_2 = hdulist_2[1].data.field('FRACEXP')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim_2 = hdulist_2[1].data.field('FRACTIME')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea_2 = hdulist_2[1].data.field('FRACAREA')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        back_2 = hdulist_2[1].data.field('BACK_COUNTS')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat_2 = hdulist_2[1].data.field('BACKRATIO')[(
            time_full_2/3600./24.+mjdref > obs_periods[i_ax][0]) *
            (time_full_2/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        # most representative value of background ratio
        backrat_med_2 = np.median(backrat_2)

        nrow_2 = len(time_2)
        start_time_2 = time_2[0]
        end_time_2 = time_2[-1]+delt_2[-1]
        time0_2 = time_2 - start_time_2
        start_time0_2 = time0_2[0]
        end_time0_2 = time0_2[-1]+delt_2[-1]
        mjds_2 = mjdref + start_time_2/24./3600.
        mjde_2 = mjdref + end_time_2/24./3600.

        # calculate new rate by error propagation assuming there are
        # enough counts in each bin to assume Gaussian statistics using
        # infirmation on total counts background exposure and backrat
        uplimit_2 = []

        # rate = []
        # rate_e = []
        xtime_2 = []
        xtime_d_2 = []
        yrate_2 = []
        yrate_e_2 = []
        start_bin_2 = time0_2[0]
        end_bin_2 = start_bin_2 + delt_2[0]
        nbin_2 = 1
        tcounts_2 = 0
        tback_2 = 0
        trate_2 = 0.
        texp_2 = 0.
        trate_e_2 = 0.
        tcount_e_2 = 0.
        tback_e_2 = 0.
        ttim_2 = 0.
        istart_2 = 0
        iend_2 = -1  # kald:for intended functionality
        counts_2 = 0
        nrate_2 = 0.
        nrate_e_2 = 0.
        nexp_2 = 0.
        bkg_2 = 0.
        ncount_e_2 = 0.
        nback_e_2 = 0.
        narea_2 = 0.
        ntime_2 = 0.
        netcounts_2 = -1.
        netcounts_e_2 = -1.
        trate_ee_2 = -1.

        for i in range(nrow_2):
            tmp_2 = time0_2[i] - start_bin_2
            tcounts_2 = tcounts_2 + cnts_2[i]
            texp_2 = texp_2 + fexp_2[i]
            ttim_2 = ttim_2 + ftim_2[i]
            tback_2 = tback_2 + back_2[i]
            netcounts_2 = (tcounts_2 - tback_2*backrat_med_2)
            trate_2 = (tcounts_2 - tback_2*backrat_med_2)/(texp_2*delt_2[0])
            tcount_e_2 = tcounts_2**0.5
            tback_e_2 = tback_2**0.5
            trate_e_2 = (tcount_e_2**2) + ((tback_e_2**2)*(backrat_med_2**2))
            netcounts_e_2 = trate_e_2**0.5
            trate_ee_2 = trate_e_2**0.5/(texp_2*delt_2[0])
            if tmp_2 > 3600.0:
                # bin finished
                xtime_2.append(0.5*(time0_2[istart_2]+time0_2[iend_2]))
                xtime_d_2.append(0.5*(time0_2[iend_2]-time0_2[istart_2]))
                end_bin_2 = time0_2[iend_2]
                counts_2 = 0
                nrate_2 = 0.
                nrate_e_2 = 0.
                nexp_2 = 0.
                bkg_2 = 0.
                ncount_e_2 = 0.
                nback_e_2 = 0.
                narea_2 = 0.
                ntime_2 = 0.

                for j in range(istart_2, (iend_2+1)):
                    nexp_2 = nexp_2 + fexp_2[j]
                    counts_2 = counts_2 + cnts_2[j]
                    bkg_2 = bkg_2 + back_2[j]
                    nrate_2 = (counts_2 - bkg_2*backrat_med_2) / \
                        (nexp_2*delt_2[0])
                    ncount_e_2 = counts_2**0.5
                    nback_e_2 = bkg_2**0.5
                    nrate_e_2 = (ncount_e_2**2) + \
                        ((nback_e_2**2)*(backrat_med_2**2))
                    narea_2 = narea_2 + farea_2[j]
                    ntime_2 = ntime_2 + ftim_2[j]
                if istart_2 == iend_2+1:
                    nexp_2 = nexp_2 + fexp_2[istart_2]
                    counts_2 = counts_2 + cnts_2[istart_2]
                    bkg_2 = bkg_2 + back_2[istart_2]
                    nrate_2 = (counts_2 - bkg_2*backrat_med_2) / \
                        (nexp_2*delt_2[0])
                    ncount_e_2 = counts_2**0.5
                    nback_e_2 = bkg_2**0.5
                    nrate_e_2 = (ncount_e_2**2) + \
                        ((nback_e_2**2)*(backrat_med_2**2))
                yrate_2.append(nrate_2)
                yrate_e_2.append(nrate_e_2**0.5/(nexp_2*delt_2[0]))
                if counts_2 < mincounts:
                    uplimit_2.append(1)
                else:
                    uplimit_2.append(0)
                # start next bin
                nbin_2 += 1
                iend_2 += 1
                if i <= nrow_2 - 1:
                    istart_2 = i
                    start_bin_2 = time0_2[istart_2]
            else:
                iend_2 += 1
        end_bin_2 = time0_2[iend_2]
        xtime_2.append(0.5*(time0_2[istart_2]+time0_2[iend_2]))
        xtime_d_2.append(0.5*(time0_2[iend_2]-time0_2[istart_2]))
        counts_2 = 0
        nrate_2 = 0.
        nrate_e_2 = 0.
        nexp_2 = 0.
        for j in range(istart_2, (iend_2+1)):
            nexp_2 = nexp_2 + fexp_2[j]
            counts_2 = counts_2 + cnts_2[j]
            bkg_2 = bkg_2 + back_2[j]
            nrate_2 = (counts_2 - bkg_2*backrat_med_2)/(nexp_2*delt_2[0])
            ncount_e_2 = counts_2**0.5
            nback_e_2 = bkg_2**0.5
            nrate_e_2 = (ncount_e_2**2) + ((nback_e_2**2)*(backrat_med_2**2))
        if istart_2 == iend_2+1:
            nexp_2 = nexp_2 + fexp_2[istart_2]
            counts_2 = counts_2 + cnts_2[istart_2]
            bkg_2 = bkg_2 + back_2[istart_2]
            nrate_2 = (counts_2 - bkg_2*backrat_med_2)/(nexp_2*delt_2[0])
            ncount_e_2 = counts_2**0.5
            nback_e_2 = bkg_2**0.5
            nrate_e_2 = (ncount_e_2**2) + ((nback_e_2**2)*(backrat_med_2**2))

        yrate_2.append(nrate_2)
        yrate_e_2.append(nrate_e_2**0.5/(nexp_2*delt_2[0]))
        if counts_2 < mincounts:
            uplimit_2.append(1)
        else:
            uplimit_2.append(0)

        # start of first bin at 0:
        xtime_2 = xtime_2-xtime_2[0]+xtime_d_2[0]
        mjd_2 = np.array(xtime_2)/3600./24. + mjds_2
        mjd_d_2 = np.array(xtime_d_2)/3600./24.
        xtime_d_2 = np.array(xtime_d_2)

        yrate_2 = np.array(yrate_2)
        yrate_e_2 = np.array(yrate_e_2)

        if short_time:
            xtime_2 = xtime_2 - time_rel
            mjd_2 = mjd_2 - time_rel

        filter_1 = [entry in mjd_2 for entry in mjd_1]
        filter_2 = [entry in mjd_1 for entry in mjd_2]

        filter = (yrate_1[filter_1]+yrate_2[filter_2] != 0)
        uplimit_1 = np.array(uplimit_1)
        uplimit_2 = np.array(uplimit_2)

        yrate_1 = yrate_1[filter_1][filter]
        yrate_e_1 = yrate_e_1[filter_1][filter]
        xtime_1 = xtime_1[filter_1][filter]
        xtime_d_1 = xtime_d_1[filter_1][filter]
        mjd_1 = mjd_1[filter_1][filter]
        mjd_d_1 = mjd_d_1[filter_1][filter]
        uplimit_1 = uplimit_1[filter_1][filter]

        yrate_2 = yrate_2[filter_2][filter]
        yrate_e_2 = yrate_e_2[filter_2][filter]
        xtime_2 = xtime_2[filter_2][filter]
        xtime_d_2 = xtime_d_2[filter_2][filter]
        mjd_2 = mjd_2[filter_2][filter]
        mjd_d_2 = mjd_d_2[filter_2][filter]
        uplimit_2 = uplimit_2[filter_2][filter]

        yrate = (yrate_1 + (-yrate_2)) / (yrate_1 + yrate_2)
        yrate_e = np.sqrt((2*yrate_2/(yrate_1+yrate_2))**2 * (yrate_e_1)**2 +
                          (2*yrate_1/(yrate_1+yrate_2))**2 * (yrate_e_2)**2)

        ###########
        ymin = min(yrate + (-yrate_e))
        ymax = max(yrate + yrate_e * np.logical_not(uplimit_1))
        xm = (xmax-xmin)*0.05
        ym = (ymax-ymin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)
        pymin.append(ymin - ym)
        pymax.append(ymax + ym)

        # Plot limits
        log.debug(
            f"Plot limits: {pxmin[i_ax]} {pxmax[i_ax]} {pymin[i_ax]} " +
            f"{pymax[i_ax]}\n")

        if xflag == 1:
            ax.errorbar(xtime_1, yrate, xerr=xtime_d_1, yerr=yrate_e,
                        uplims=uplimit_1, linestyle='None', color=color, fmt='o',
                        zorder=1, lolims=uplimit_2)
        else:
            ax.errorbar(mjd_1, yrate, xerr=mjd_d_1, yerr=yrate_e,
                        uplims=uplimit_1, linestyle='None', color=color, fmt='o',
                        zorder=1, lolims=uplimit_2)

    return pxmin, pxmax, pymin, pymax, time_rel
