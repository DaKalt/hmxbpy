import numpy as np
from HiMaXBipy.io.package_data import round_to_1
import matplotlib.ticker as plticker
from matplotlib import ticker, rcParams
import matplotlib.pyplot as plt


def plot_lc_UL(hdulist, ax, logfile, mjdref, xflag, mincounts, color):
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
    logfile.write(f'Total number of time bins in table: {nrow}\n')
    start_time = time[0]
    end_time = time[-1]+delt[-1]
    time0 = time - start_time
    start_time0 = time0[0]
    end_time0 = time0[-1]+delt[-1]
    mjds = mjdref + start_time/24./3600.
    mjde = mjdref + end_time/24./3600.

    logfile.write(
        f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
    logfile.write(f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

    # calculate new rate by error propagation assuming there are enough counts in each bin to assume Gaussian statistics usinf infirmation on total counts background exposure and backrat
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
            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
            if counts < mincounts:
                uplimit.append(1)
            else:
                uplimit.append(0)
            logfile.write(
                f'bin number, Start,end time of bin, duration of bin, counts , rate, error,fracexp, fractime, fracarea, loop start, loop end {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate} {nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} {istart} {iend}\n')
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

    yrate.append(nrate)
    yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
    if counts < mincounts:
        uplimit.append(1)
    else:
        uplimit.append(0)

    logfile.write(
        f'Start/end time of last bin {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} {nexp} {istart} {iend}\n')
    # start of first bin at 0:
    xtime = xtime-xtime[0]+xtime_d[0]
    logfile.write(f'{xtime} {xtime_d}\n')
    mjd = np.array(xtime)/3600./24. + mjds
    mjd_d = np.array(xtime_d)/3600./24.

    logfile.write(f'Total counts: {tcounts} +/- {tcount_e}\n')
    logfile.write(f'Net counts {netcounts} +/- {netcounts_e}\n')
    logfile.write(f'binsize is {delt[0]}\n')
    logfile.write(f'Average rate: {trate} +/- {trate_ee}\n')

    logfile.write(f'Total exposure: {ttim}\n')
    logfile.write(f'Total fract. exposure: {texp}\n')

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    if xflag == 1:
        xmin = min(xtime)
        xmax = max(xtime)
    else:
        xmin = min(mjd)
        xmax = max(mjd)

    ymin = min(yrate - yrate_e)
    ymax = max(yrate + yrate_e * np.logical_not(uplimit))
    xm = (xmax-xmin)*0.05
    ym = (ymax-ymin)*0.05
    pxmin = xmin - xm
    pxmax = xmax + xm
    pymin = ymin - ym
    pymax = ymax + ym

    # Plot limits
    logfile.write(f"Plot limits: {pxmin} {pxmax} {pymin} {pymax}\n")

    if xflag == 1:
        ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                    uplims=uplimit, linestyle='None', color=color, fmt='o', zorder = 1)
    else:
        ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                    uplims=uplimit, linestyle='None', color=color, fmt='o', zorder=1)

    return pxmin, pxmax, pymin, pymax


def plot_lc_mincounts(hdulist, ax, logfile, mjdref, xflag, mincounts, color):
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
    logfile.write(f'Total number of time bins in table: {nrow}\n')
    start_time = time[0]
    end_time = time[-1]+delt[-1]
    time0 = time - start_time
    start_time0 = time0[0]
    end_time0 = time0[-1]+delt[-1]
    mjds = mjdref + start_time/24./3600.
    mjde = mjdref + end_time/24./3600.

    logfile.write(
        f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
    logfile.write(f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

    # calculate new rate by error propagation assuming there are enough counts in each bin to assume Gaussian statistics usinf infirmation on total counts background exposure and backrat
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
    # mean_count = 0
    # sn = 2.0
    # xtime_tmp = 0.0
    # xtime_d_start_tmp = 0.0
    logfile.write(f'Minimum number of counts is {mincounts}\n')
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
            iend += 1  # kald: in this case one more bin needs to be included compared to the usual lightcurves, else <10 countsfor all newly defined bins
            # without the +1 before but this way sections are directly connected #kald
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

            logfile.write(
                f'bin number, Start,end time of bin, duration of bin, counts , rate, error,fracexp, fractime, fracarea, loop start, loop end {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate} {nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} {istart} {iend}\n')
            # start next bin
            nbin += 1
            istart_tmp = istart
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
    logfile.write(
        f'Start/end time of last bin {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} {nexp} {istart} {iend}\n')
    if counts > mincounts:  # this was tcounts_c which should be 0 all the time at this point #kald
        logfile.write(f'Condition satisfied with {counts} {bkg}\n')

        yrate.append(nrate)
        yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
        # see above comment about connected #kald
        xtime.append(0.5*(time0[istart]+time0[iend]))
        xtime_d.append(0.5*(time0[iend]-time0[istart]))
    else:
        logfile.write('averaging\n')
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

    logfile.write(f'Total counts: {tcounts} +/- {tcount_e}\n')
    logfile.write(f'Net counts {netcounts} +/- {netcounts_e}\n')
    logfile.write(f'binsize is {delt[0]}\n')
    logfile.write(f'Average rate: {trate} +/- {trate_ee}\n')

    logfile.write(f'Total exposure: {ttim}\n')
    logfile.write(f'Total fract. exposure: {texp}\n')

    if xflag == 1:
        xmin = xtime[0] - xtime_d[0]
        xmax = xtime[-1] + xtime_d[-1]
    else:
        xmin = mjd[0] - mjd_d[0]
        xmax = mjd[-1] + mjd_d[-1]

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    ymin = min(yrate - yrate_e)
    ymax = max(yrate + yrate_e)
    xm = (xmax-xmin)*0.05
    ym = (ymax-ymin)*0.05
    pxmin = xmin - xm
    pxmax = xmax + xm
    pymin = ymin - ym
    pymax = ymax + ym

    # Plot limits
    logfile.write(f"Plot limits: {pxmin} {pxmax} {pymin} {pymax}\n")

    if xflag == 1:
        ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                    linestyle='None', color=color, fmt='o', zorder = 2)
    else:
        ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                    linestyle='None', color=color, fmt='o', zorder=2)

    yrate = np.array(yrate)
    yrate_e = np.array(yrate_e)

    i_max = np.argmax(yrate - yrate_e)
    i_min = np.argmax(yrate_e - yrate)

    ampl_max2 = yrate[i_max] - yrate_e[i_max] - yrate[i_min] - yrate_e[i_min]
    ampl_max = yrate[i_max] - yrate[i_min]
    if yrate[i_min] > 0:
        variability = yrate[i_max] / yrate[i_min]
    else:
        variability = -1
    ampl_sig = ampl_max / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)
    ampl_sig2 = ampl_max2 / np.sqrt(yrate_e[i_max] ** 2 + yrate_e[i_min] ** 2)

    logfile.write(f'AMPL_MAX: {ampl_max}\n')
    logfile.write(f'Variability V = {variability}\n')
    logfile.write(f'AMPL_SIG: {ampl_sig}\n')
    logfile.write(f'AMPL_MAX conservative: {ampl_max2}\n')
    logfile.write(f'AMPL_SIG2: {ampl_sig2}\n')

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
        if i * tick_size_x + centre_x > pxmin and i * tick_size_x + centre_x < pxmax:
            xticks.append(i * tick_size_x + centre_x)
    ax.set_xticks(xticks)


def standardize_ticks(axs, xbase=None, ybase=None):
    """Make all of the internal axes share tick bases
    Parameters
    ----------
    xbase, ybase: (optional) None or float
        If `xbase` or `ybase` is a float, manually set all tick locators to
        this base. Otherwise, use the largest base across internal subplots
        for that axis.
    """
    if xbase is None:
        if axs[0].xaxis.get_scale() == "log":
            xbase = max(
                ax.xaxis.get_ticklocs()[1] / ax.xaxis.get_ticklocs()[0]
                for ax in axs
                if ax.get_subplotspec().is_last_row()
            )
        else:
            xbase = max(
                ax.xaxis.get_ticklocs()[1] - ax.xaxis.get_ticklocs()[0]
                for ax in axs
                if ax.get_subplotspec().is_last_row()
            )
    # if ybase is None:
    #     if axs[0].yaxis.get_scale() == "log":
    #         ybase = max(
    #             ax.yaxis.get_ticklocs()[1] / ax.yaxis.get_ticklocs()[0]
    #             for ax in axs
    #             if ax.get_subplotspec().is_first_col()
    #         )
    #     else:
    #         ybase = max(
    #             ax.yaxis.get_ticklocs()[1] - ax.yaxis.get_ticklocs()[0]
    #             for ax in axs
    #             if ax.get_subplotspec().is_first_col()
    #         )

    for ax in axs:
        # if ax.get_subplotspec().is_first_col():
        #     if ax.yaxis.get_scale() == "log":
        #         ax.yaxis.set_major_locator(ticker.LogLocator(ybase))
        #     else:
        #         ax.yaxis.set_major_locator(ticker.MultipleLocator(ybase))
        if ax.get_subplotspec().is_last_row():
            if ax.xaxis.get_scale() == "log":
                ax.xaxis.set_major_locator(ticker.LogLocator(xbase))
            else:
                ax.xaxis.set_major_locator(ticker.MultipleLocator(xbase))


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


def plot_lc_mincounts_broken_new(hdulist, axs, logfile, mjdref, xflag, mincounts, color, obs_periods, short_time, time_rel=0):
    pxmax = []
    pxmin = []
    pymax = []
    pymin = []
    for i_ax, ax in enumerate(axs):
        time_full = hdulist[1].data.field('TIME')
        time = hdulist[1].data.field('TIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt = hdulist[1].data.field('TIMEDEL')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts = hdulist[1].data.field('COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp = hdulist[1].data.field('FRACEXP')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim = hdulist[1].data.field('FRACTIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea = hdulist[1].data.field('FRACAREA')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        back = hdulist[1].data.field('BACK_COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat = hdulist[1].data.field('BACKRATIO')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        backrat_med = np.median(backrat)

        nrow = len(time)
        logfile.write(f'Total number of time bins in table: {nrow}\n')
        start_time = time[0]
        end_time = time[-1]+delt[-1]
        time0 = time - start_time
        start_time0 = time0[0]
        end_time0 = time0[-1]+delt[-1]
        mjds = mjdref + start_time/24./3600.
        mjde = mjdref + end_time/24./3600.

        logfile.write(
            f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
        logfile.write(
            f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

        # calculate new rate by error propagation assuming there are enough counts in each bin to assume Gaussian statistics usinf infirmation on total counts background exposure and backrat
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
        # mean_count = 0
        # sn = 2.0
        # xtime_tmp = 0.0
        # xtime_d_start_tmp = 0.0
        logfile.write(f'Minimum number of counts is {mincounts}\n')
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
                iend += 1  # kald: in this case one more bin needs to be included compared to the usual lightcurves, else <10 countsfor all newly defined bins
                # without the +1 before but this way sections are directly connected #kald
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

                logfile.write(
                    f'bin number, Start,end time of bin, duration of bin, counts , rate, error,fracexp, fractime, fracarea, loop start, loop end {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate} {nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} {istart} {iend}\n')
                # start next bin
                nbin += 1
                istart_tmp = istart
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
        logfile.write(
            f'Start/end time of last bin {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} {nexp} {istart} {iend}\n')
        if counts > mincounts:  # this was tcounts_c which should be 0 all the time at this point #kald
            logfile.write(f'Condition satisfied with {counts} {bkg}\n')

            yrate.append(nrate)
            yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
            # see above comment about connected #kald
            xtime.append(0.5*(time0[istart]+time0[iend]))
            xtime_d.append(0.5*(time0[iend]-time0[istart]))
        else:
            logfile.write('averaging\n')
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

        logfile.write(f'Total counts: {tcounts} +/- {tcount_e}\n')
        logfile.write(f'Net counts {netcounts} +/- {netcounts_e}\n')
        logfile.write(f'binsize is {delt[0]}\n')
        logfile.write(f'Average rate: {trate} +/- {trate_ee}\n')

        logfile.write(f'Total exposure: {ttim}\n')
        logfile.write(f'Total fract. exposure: {texp}\n')

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

        ymin = min(yrate - yrate_e)
        ymax = max(yrate + yrate_e)
        xm = (xmax-xmin)*0.05
        ym = (ymax-ymin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)
        pymin.append(ymin - ym)
        pymax.append(ymax + ym)

        # Plot limits
        logfile.write(
            f"Plot limits: {pxmin[i_ax]} {pxmax[i_ax]} {pymin[i_ax]} {pymax[i_ax]}\n")

        if xflag == 1:
            ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                        linestyle='None', color=color, fmt='o', zorder=2)
        else:
            ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                        linestyle='None', color=color, fmt='o', zorder=2)

        yrate = np.array(yrate)
        yrate_e = np.array(yrate_e)

        i_max = np.argmax(yrate - yrate_e)
        i_min = np.argmax(yrate_e - yrate)

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

        logfile.write(f'AMPL_MAX: {ampl_max}\n')
        logfile.write(f'Variability V = {variability}\n')
        logfile.write(f'AMPL_SIG: {ampl_sig}\n')
        logfile.write(f'AMPL_MAX conservative: {ampl_max2}\n')
        logfile.write(f'AMPL_SIG2: {ampl_sig2}\n')

    return pxmin, pxmax, pymin, pymax, time_rel


def plot_lc_UL_broken_new(hdulist, axs, logfile, mjdref, xflag, mincounts, color, obs_periods, short_time, time_rel=0):
    pxmax = []
    pxmin = []
    pymax = []
    pymin = []
    for i_ax, ax in enumerate(axs):
        time_full = hdulist[1].data.field('TIME')
        time = hdulist[1].data.field('TIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        delt = hdulist[1].data.field('TIMEDEL')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        cnts = hdulist[1].data.field('COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        fexp = hdulist[1].data.field('FRACEXP')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        ftim = hdulist[1].data.field('FRACTIME')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        farea = hdulist[1].data.field('FRACAREA')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        back = hdulist[1].data.field('BACK_COUNTS')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        backrat = hdulist[1].data.field('BACKRATIO')[(
            time_full/3600./24.+mjdref > obs_periods[i_ax][0]) * (time_full/3600./24.+mjdref < obs_periods[i_ax][1])]
        # delta = delt[0]
        # most representative value of background ratio
        backrat_med = np.median(backrat)

        nrow = len(time)
        logfile.write(f'Total number of time bins in table: {nrow}\n')
        start_time = time[0]
        end_time = time[-1]+delt[-1]
        time0 = time - start_time
        start_time0 = time0[0]
        end_time0 = time0[-1]+delt[-1]
        mjds = mjdref + start_time/24./3600.
        mjde = mjdref + end_time/24./3600.

        logfile.write(
            f'Start of 1st bin : {start_time} {start_time0} MJD: {mjds}\n')
        logfile.write(
            f'End of last bin  : {end_time} {end_time0} MJD: {mjde}\n')

        # calculate new rate by error propagation assuming there are enough counts in each bin to assume Gaussian statistics usinf infirmation on total counts background exposure and backrat
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
                yrate.append(nrate)
                yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
                if counts < mincounts:
                    uplimit.append(1)
                else:
                    uplimit.append(0)
                logfile.write(
                    f'bin number, Start,end time of bin, duration of bin, counts , rate, error,fracexp, fractime, fracarea, loop start, loop end {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate} {nrate_e**0.5/(nexp*delt[0])} {nexp} {ntime} {narea} {istart} {iend}\n')
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

        yrate.append(nrate)
        yrate_e.append(nrate_e**0.5/(nexp*delt[0]))
        if counts < mincounts:
            uplimit.append(1)
        else:
            uplimit.append(0)

        logfile.write(
            f'Start/end time of last bin {nbin} {start_bin} {end_bin} {end_bin-start_bin} {counts} {nrate/nexp} {nrate_e**0.5/nexp} {nexp} {istart} {iend}\n')
        # start of first bin at 0:
        xtime = xtime-xtime[0]+xtime_d[0]
        logfile.write(f'{xtime} {xtime_d}\n')
        mjd = np.array(xtime)/3600./24. + mjds
        mjd_d = np.array(xtime_d)/3600./24.

        logfile.write(f'Total counts: {tcounts} +/- {tcount_e}\n')
        logfile.write(f'Net counts {netcounts} +/- {netcounts_e}\n')
        logfile.write(f'binsize is {delt[0]}\n')
        logfile.write(f'Average rate: {trate} +/- {trate_ee}\n')

        logfile.write(f'Total exposure: {ttim}\n')
        logfile.write(f'Total fract. exposure: {texp}\n')

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

        ymin = min(yrate - yrate_e)
        ymax = max(yrate + yrate_e * np.logical_not(uplimit))
        xm = (xmax-xmin)*0.05
        ym = (ymax-ymin)*0.05
        pxmin.append(xmin - xm)
        pxmax.append(xmax + xm)
        pymin.append(ymin - ym)
        pymax.append(ymax + ym)

        # Plot limits
        logfile.write(
            f"Plot limits: {pxmin[i_ax]} {pxmax[i_ax]} {pymin[i_ax]} {pymax[i_ax]}\n")

        if xflag == 1:
            ax.errorbar(xtime, yrate, xerr=xtime_d, yerr=yrate_e,
                        uplims=uplimit, linestyle='None', color=color, fmt='o', zorder=1)
        else:
            ax.errorbar(mjd, yrate, xerr=mjd_d, yerr=yrate_e,
                        uplims=uplimit, linestyle='None', color=color, fmt='o', zorder=1)

    return pxmin, pxmax, pymin, pymax, time_rel


def format_axis_broken_new(fig, axs, pxmins, pxmaxs, pymin, pymax, ticknumber_x,
                           ticknumber_y, ncols, nrows, d, tilt, diag_color, big_ax):
    # proportion of vertical to horizontal extent of the slanted line
    d = np.tan(tilt)
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color=diag_color, mec=diag_color, mew=1, clip_on=False)

    # fig.subplots_adjust(wspace=0.05)

    for i_ax, ax in enumerate(axs):
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

        if i_ax == 0:
            ax.plot([1, 1], [0, 1], transform=ax.transAxes, **kwargs)
            ax.spines.right.set_visible(False)
            ax.tick_params(right=False, labelright=False)
        elif i_ax == len(axs) - 1:
            ax.plot([0, 0], [0, 1], transform=ax.transAxes, **kwargs)
            ax.spines.left.set_visible(False)
            ax.tick_params(left=False, labelleft=False)
        else:
            ax.plot([0, 0, 1, 1], [0, 1, 0, 1],
                    transform=ax.transAxes, **kwargs)
            ax.spines.right.set_visible(False)
            ax.spines.left.set_visible(False)
            ax.tick_params(left=False, labelleft=False,
                           right=False, labelright=False)

        tick_size_y = round_to_1((pymax - pymin) / ticknumber_y)
        yticks = []
        for j in range(-int(ticknumber_y), int(ticknumber_y)):
            if j * tick_size_y > pymin and j * tick_size_y < pymax:
                yticks.append(j * tick_size_y)
        ax.set_yticks(yticks)

        tick_size_x = np.round(
            (pxmaxs[i_ax] - pxmins[i_ax]) / ticknumber_x)
        xticks = []
        centre_x = np.round(
            (pxmaxs[i_ax] + pxmins[i_ax]) / 2.)
        for j in range(-int(ticknumber_x), int(ticknumber_x)):
            if j * tick_size_x + centre_x > pxmins[i_ax] and j * tick_size_x + centre_x < pxmaxs[i_ax]:
                xticks.append(j * tick_size_x + centre_x)
        if len(xticks) <= 1:
            xticks = [centre_x - tick_size_x / 2, centre_x + tick_size_x / 2]
        ax.set_xticks(xticks)

        ax.set_xbound([pxmins[i_ax], pxmaxs[i_ax]])
        ax.set_ybound([pymin, pymax])

    longest_y = ''
    for entry in axs[0].yaxis.get_ticklabels():
        if len(entry.get_text()) > len(longest_y):
            longest_y = entry.get_text()
    start_x = axs[0].xaxis.get_ticklabels()[0].get_text()
    end_x = axs[-1].xaxis.get_ticklabels()[-1].get_text()
    print(axs[-1].xaxis.get_ticklabels()[-1])

    big_ax.set_xbound([0, 1])
    big_ax.set_ybound([0, 1])
    big_ax.set_xticks([0, 1])
    big_ax.set_yticks([0, 1])
    big_ax.tick_params(left=False, bottom=False, right=False, top=False)
    big_ax.set_xticklabels([start_x, end_x], alpha=0.3)
    big_ax.set_yticklabels([longest_y, longest_y], alpha=0.3)

    print([start_x, end_x], [longest_y, longest_y])

    fig.set_tight_layout(True)
    fig.set_tight_layout(True)
    fig.subplots_adjust(wspace=0.05)

    # for i, ax in enumerate(axs):
    #     if obs_periods is not None:
    #         ax.set_xlim(obs_periods[i % ncols])
    #         ax.get_shared_x_axes().join(ax, last_row[i % ncols])
    # standardize_ticks(axs) #not sure if necessary
