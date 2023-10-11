#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 04:23:18 2022

@author: David Kaltenbrunner
"""
import fileinput
import getpass
import os
import shutil
import subprocess
import sys
import time
import warnings

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import scipy.stats
try:
    from xspec import AllData, AllModels, Fit, Model, Plot, PlotManager,\
        Spectrum, Xset
except ModuleNotFoundError:
    raise Exception('xspec needs to be installed before using this package. '
                    'xspec is fully integrated in the general build of '
                    'Heasoft, follow instruction on '
                    'http://heasarc.gsfc.nasa.gov/lheasoft/install.html')

from HiMaXBipy.io.package_data import get_path_of_data_dir, get_stan_dir
from HiMaXBipy.io.logging import setup_logfile, setup_logger, set_loglevel
from HiMaXBipy.lc_plotting.lc_plotting import plot_lc_UL, plot_lc_mincounts,\
    get_boundaries, get_boundaries_broken, format_axis, plot_lc_UL_broken_new,\
    plot_lc_mincounts_broken_new, format_axis_broken_new, plot_lc_UL_hr,\
    plot_lc_mincounts_hr
from HiMaXBipy.lc_plotting.lc_plotting_bayes import plot_lc_eROday_broken_bayes,\
    plot_lc_mincounts_broken_bayes
from HiMaXBipy.spectral_analysis.spectral_analysis import spec_model
from HiMaXBipy.spectral_analysis.spectral_analysis_bxa import fit_bxa, plot_bxa
from HiMaXBipy.spectral_analysis.standard_models_bxa import apl, apl_simple


class HiMaXBi:
    _skytile = ''
    _filelist = ''
    _esass = '/home/erosita/sw/eSASSusers_211214/bin/esass-init.sh'
    _LC_prebinning = '1.0'
    _LC_extracted = False
    _mjdref = 51543.875
    _ero_starttimes = np.array([58828, 59011, 59198, 59381, 59567])
    _energy_bins = [[0.2, 8.0]]
    _energy_bins_hr = [[0.2, 2.0], [2.0, 8.0]]
    _grouping = 1
    _ownership = 'x'
    _distance = 50.
    _Z = 0.49
    _debugging = False
    _NH = 6e20
    _NH_set = False
    _create_epochs = False

    def __init__(self, src_name, working_dir, data_dir, fix_path=True):
        '''
        Parameters
        ----------
        src_name : str
            Sets the name of the source used for naming files.
        working_dir : str
            Sets the working directory to which resulting files will be
            saved.
        data_dir : str
            Sets directory from where to use data files.
        fix_path : bool
            If set true, replaces '*/galaxy' by 'data40s/galaxy' in
            paths

        '''
        if type(src_name) == str:
            if src_name.find('-') != -1 or src_name.find(" ") != -1:
                warnings.warn(
                    'There can be problems when using "-" or " " as part of '
                    'the source name. All "-" and " " are replaced by "_".')
                src_name = src_name.replace('-', '_')
                src_name = src_name.replace(' ', '_')
            self._src_name = src_name
        else:
            raise Exception('src_name needs to be of type string.')

        working_dir = os.path.abspath(os.path.expanduser(working_dir))
        data_dir = os.path.abspath(os.path.expanduser(data_dir))

        if fix_path:
            working_dir = '/data40s/' + \
                working_dir[working_dir.find('galaxy'):]
            data_dir = '/data40s/' + data_dir[data_dir.find('galaxy'):]

        if os.path.exists(working_dir) and type(working_dir) == str:
            self._working_dir_full = os.path.abspath(working_dir)
        else:
            raise Exception('Not a valid path for a working directory.')

        if (os.path.exists(data_dir) and type(data_dir) == str):
            self._data_dir_full = os.path.abspath(data_dir)
        else:
            raise Exception('Not a valid path for a data directory.')

        for subdir in ['/working', '/results', '/logfiles']:
            if not os.path.exists(self._working_dir_full + subdir):
                os.mkdir(self._working_dir_full + subdir)
        for subdir in ['/results', '/logfiles']:
            for subsubdir in ['/lightcurves', '/spectra']:
                if not os.path.exists(self._working_dir_full + subdir
                                      + subsubdir):
                    os.mkdir(self._working_dir_full + subdir + subsubdir)

        self._sh_dir_ = get_path_of_data_dir()
        self._stan_dir_ = get_stan_dir()

        for (path, directories, filenames) in os.walk(self._sh_dir_):
            for filename in filenames:
                if os.path.exists(f'{self._working_dir_full}/working/'
                                  f'{filename}'):
                    os.remove(f'{self._working_dir_full}/working/{filename}')
                shutil.copy(self._sh_dir_ + '/' + filename,
                            self._working_dir_full + '/working')

        for (path, directories, filenames) in os.walk(self._stan_dir_):
            for filename in filenames:
                if os.path.exists(f'{self._working_dir_full}/working/'
                                  f'{filename}'):
                    os.remove(f'{self._working_dir_full}/working/{filename}')
                shutil.copy(self._stan_dir_ + '/' + filename,
                            self._working_dir_full + '/working')

        os.chdir(working_dir + '/working/')

        self._working_dir = os.path.relpath(working_dir)
        self._data_dir = os.path.relpath(data_dir)

        if self._working_dir.find("-") != -1 or self._data_dir.find("-") != -1:
            raise Exception(
                'Working and Data directories with a "-" in their full path ' +
                'cause problems during data analysis.')

        self._logger = setup_logger('HiMaXBipy', self._working_dir_full)

    def change_LogLevel(self, level):
        '''Change the level of all logs (files and stdout) to level.
        Parameters
        ----------
        level : Level as described by logging package

        '''
        set_loglevel(self._logger, level)

    def _replace_in_sh(self, path, replacements):
        '''
        Parameters
        ----------
        replacements : array-like shape (n, 2)
            array of n keywords to replace in sh file; 0th entry in each
            pair states the original keyword, 1st entry states the new
            keyword
        path : str
            path of sh file in which to replace entries

        Returns
        -------
        Name of the new sh file as a string.

        '''
        new_sh = path[:-3] + '_modified.sh'
        if os.path.exists(new_sh):
            os.remove(new_sh)
        shutil.copy(path, new_sh)
        for pair in replacements:
            for line in fileinput.input(new_sh, inplace=True):
                line = line.replace(str(pair[0]), str(pair[1]))
                sys.stdout.write(line)
        return new_sh

    def set_Ebins(self, bins):
        '''Set energy bins to be analysed in keV. The default is
        [0.2, 8.0].

        Parameters
        ----------
        bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin
            E_min and E_max must be given in keV. The default is
            [[0.2, 8.0]]
        '''
        self._logger.info('Ebins set.')
        if type(bins) != list and type(bins) != np.ndarray:
            raise Exception('bins must be array-like')
        else:
            for line in bins:
                if len(line) != 2:
                    raise Exception('Each line in bins needs 2 entries.')
                if ((type(line[0]) != float and type(line[0]) != int) or
                        (type(line[1]) != float and type(line[1]) != int)):
                    raise Exception(
                        'The entries of each line of bins need to be the '
                        'minimum and maximum energies given in keV of the '
                        'energy bins to analyse given as float or int.')
                elif (line[0] < 0.2 or line[0] > 8.0 or line[1] < 0.2 or
                      line[1] > 8.0 or line[0] >= line[1]):
                    raise Exception(
                        'The energies must be given in keV and must follow '
                        '0.2 <= E_min < E_max <= 8.0.')
        self._energy_bins = np.array(bins, dtype=np.float64).tolist()
        self._LC_extracted = False

    def set_Ebins_HR(self, bins):
        '''Set energy bins to be analysed in keV. The default is
        [0.2, 8.0].

        Parameters
        ----------
        bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin
            E_min and E_max must be given in keV. The default is
            [[0.2, 8.0]]
        '''
        self._logger.info('Ebins HR set.')
        if type(bins) != list and type(bins) != np.ndarray:
            raise Exception('bins must be array-like')
        else:
            if len(bins) != 2:
                raise Exception(
                    'Currently only a set of two energy bins is supported.')
            for line in bins:
                if len(line) != 2:
                    raise Exception('Each line in bins needs 2 entries.')
                if ((type(line[0]) != float and type(line[0]) != int) or
                        (type(line[1]) != float and type(line[1]) != int)):
                    raise Exception(
                        'The entries of each line of bins need to be the '
                        'minimum and maximum energies given in keV of the '
                        'energy bins to analyse given as float or int.')
                elif (line[0] < 0.2 or line[0] > 8.0 or line[1] < 0.2 or
                      line[1] > 8.0 or line[0] >= line[1]):
                    raise Exception(
                        'The energies must be given in keV and must follow '
                        '0.2 <= E_min < E_max <= 8.0.')
        self._energy_bins_hr = np.array(bins, dtype=np.float64).tolist()
        self._LC_extracted = False

    def set_distance(self, distance):
        '''Set distance to source in kpc. The default is 50.

        Parameters
        ----------
        distance : float or str
            Distance to observed object in kpc to calculate Flux.

        '''
        self._logger.info(f'Distance set to {distance} kpc.')
        if (type(distance) != float
            and type(distance) != str
                and type(distance) != int):
            raise Exception('distance must be a float, int or string.')
        try:
            if float(distance) <= 0:
                raise Exception('distance must be > 0.')
            self._distance = float(distance)
        except TypeError:
            raise Exception('distance must be convertible to float.')

    def set_metallicity(self, Z=-1., Z_file=''):
        '''Set metallicity at the source location. The default is 0.49
        for LMC. Either Z or Z_file has to be given, if both are given,
        Z will be used.

        Parameters
        ----------
        Z : float>0 or str
            Metallicity for local absorption. The default value is -1.
        file : str
            File that contains metallicity for local absorption. The
            default value is ''.

        '''
        self._logger.info(f'Metallicity set with Z {Z} and file \'{Z_file}\'.')
        if type(Z) != float and type(Z) != str:
            raise Exception('Z must be a float or string.')
        if type(Z_file) != str:
            raise Exception('Z_file must be a string.')
        try:
            if float(Z) == -1. and os.path.exists(Z_file):
                with open(Z_file) as file:
                    lines = file.readlines()
                Z = lines[0][:lines[0].find('\n')]
            if float(Z) <= 0:
                raise Exception('Z must be > 0 or Z_file must exist.')
        except TypeError:
            raise Exception('Z must be convertible to float.')
        self._Z = float(Z)

    def get_NH(self):
        '''Create DL.NH file with weighted NH from NH map by
        Dickey & Lockman  1990 at the source position and set NH to the
        value given with the task set_NH(NH_file='DL.NH').
        The task set_radec with the source position must be run before.
        '''
        if not '_RA' in dir(self):
            raise Exception('Source position needs to be set first.')
        self._logger.info('Creating DL.NH')
        sh_file = f'{self._working_dir_full}/working/get_nh.sh'
        process = subprocess.Popen(
            [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()  # Wait for process to complete.
        NH = ''
        string = '(cm**-2) '
        for line in process.stdout.readlines():
            line = str(line)[2:-3]
            if line.find('Weighted') > -1:
                NH = line[line.find(string)+len(string)+1:]
        if NH == '':
            raise Exception('Something went wrong when running get_nh.sh. '
                            'A likely error is not having installed the'
                            ' heasoft software or a wrong setup. Try if the '
                            'command "nh" is working correctly.')
        with open(f'{self._working_dir_full}/DL.NH', 'w') as file:
            file.writelines(NH)
        self.set_NH(NH_file=f'{self._working_dir_full}/DL.NH')

    def set_NH(self, NH=-1., NH_file=''):
        '''Set NH in the line of sight (LOS) towards the source. The
        default for LMC is 6e20. Either NH or NH_file has to be given,
        if both are given, NH will be used.

        Parameters
        ----------
        NH : float or str
            NH in LOS to source. The default value is -1.
        file : str
            File that contains NH in LOS to source. The
            default value is ''.

        '''
        self._logger.info(f'NH set with NH {NH} and file \'{NH_file}\'.')
        if type(NH) != float and type(NH) != str:
            raise Exception('NH must be a float or string.')
        if type(NH_file) != str:
            raise Exception('NH_file must be a string.')
        try:
            if float(NH) == -1 and os.path.exists(NH_file):
                with open(NH_file) as file:
                    lines = file.readlines()
                if lines[0].find('\n') > 0:
                    NH = lines[0][:lines[0].find('\n')]
                else:
                    NH = lines[0]
            if float(NH) <= 0:
                raise Exception('NH must be > 0 or NH_file must exist.')
        except TypeError:
            raise Exception('NH must be convertible to float (NH file might '
                            'be corrupted).')
        self._NH = float(NH)
        self._NH_set = True

    def set_mjd_ref(self, mjdref):
        '''Set Reference MJD date for eROSITA times. The default is
        51543.875.

        Parameters
        ----------
        mjdref : float or str
            Reference value to transform eROSITA times to MJD.

        '''
        self._logger.info(f'MJD set {mjdref}.')
        if type(mjdref) != float and type(mjdref) != str:
            raise Exception('mjdref must be a float or string.')
        try:
            self._mjdref = float(mjdref)
        except TypeError:
            raise Exception('mjdref must be convertible to float.')

    def set_esass(self, esass_location):
        '''Set path to eSASS in use.

        Parameters
        ----------
        esass_location : str
            full path of esass initialisation script.

        '''
        self._logger.info(f'eSASS set {esass_location}.')
        if type(esass_location) != str:
            raise Exception('esass_location must be string.')
        if not os.path.exists(esass_location):
            raise Exception(f'File {esass_location} does not exist.')
        self._esass = esass_location
        self._LC_extracted = False

    def set_skytile(self, skytile):
        '''Sets name of the skytile (e.g. 080156) in which the source
        lies.

        Parameters
        ----------
        skytile : str
            Name of the skytile in which the source lies.

        '''
        self._logger.info(f'Skytile set to {skytile}.')
        if type(skytile) == str:
            self._skytile = skytile
            self._LC_extracted = False
        else:
            raise Exception('Not a valid skytile name (must be string).')

    def set_radec(self, RA, Dec):
        '''Sets RA and Dec for the source.

        Parameters
        ----------
        RA : float
            Right ascension of the source in J2000.
        Dec : TYPE
            Declination of the source in J2000.

        '''
        if type(RA) != float or type(Dec) != float:
            raise Exception('Ra and Dec must be given as floats.')
        self._RA = RA
        self._Dec = Dec
        self._LC_extracted = False

    def set_filelist(self, filelist):
        '''Set the names of eventfiles to use for analysis.

        Parameters
        ----------
        filelist : str
            List of names of eventfiles to use separated by spaces and
            without foldernames.
        '''
        self._logger.info('Filelist set.')
        if type(filelist) != str:
            raise Exception('filelist must be string.')
        filelist = filelist.strip()
        if filelist[1] == 'b':
            self._ownership = 'b'
        elif filelist[1] == 'm':
            self._ownership = 'm'
        else:
            self._ownership = 'x'
            self._logger.warning('Unknown ownership.')
        temp1 = filelist
        temp3 = ''
        while temp1.find(' ') != -1:
            temp2 = temp1[:temp1.find(' ')]
            temp1 = temp1[temp1.find(' '):].strip()
            if not os.path.exists(self._data_dir + '/' + temp2):
                raise Exception(f'File {temp2} does not exist.')
            temp3 += self._data_dir + '/' + temp2 + ' '
        if not os.path.exists(self._data_dir + '/' + temp1):
            raise Exception(f'File {temp1} does not exist.')
        temp3 += self._data_dir + '/' + temp1
        self._filelist = temp3
        self._LC_extracted = False

    def set_LC_binning(self, lc_binning):
        '''Set the initial binning of lightcurves for extraction in
        seconds. The default is 1s.

        Parameters
        ----------
        binning : str or float
            Set the initial binning of the lightcurve in seconds.

        '''
        self._logger.info('Binning LC.')
        if not (type(lc_binning) == str or type(lc_binning) == float):
            raise Exception('lc_binning must be a string or float.')
        try:
            float(lc_binning)
        except ValueError:
            raise Exception('lc_binning must be a number.')
        if float(lc_binning) <= 0:
            raise Exception('lc_binning must be >0.')
        self._LC_prebinning = str(lc_binning)
        self._LC_extracted = False

    def set_grouping(self, grouping):
        '''Sets grouping of events per energy bin for spectral analysis.
        The default is 1.

        Parameters
        ----------
        grouping : str or int
            Set the grouping of extracted spectra. If c-statistic is
            used grouping=1 is recommended, for chi2 statistic
            grouping=20 is recommended.
        '''
        self._logger.info('Setting Grouping.')
        if not (type(grouping) == str or type(grouping) == int):
            raise Exception('binning must be a string or float.')
        try:
            int(grouping)
        except ValueError:
            raise Exception('grouping must be an integer.')
        if int(grouping) < 1:
            raise Exception('grouping must be an integer >= 1.')
        if float(grouping) != int(grouping):
            self._logger.warning(
                'grouping is treated as an integer. Float values will'
                ' be rounded down to the next integer.')
        self._grouping = int(grouping)

    def _extract_lc(self, logname='lc_extract_autosave'):
        '''
        Parameters
        ----------
        logname : str
            name of logfile to safe output of sh script to.

        Returns
        -------
        Creates sh files to extract light curve fits files and runs
        them.

        '''
        self._logger.info('Extracting LC.')
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if self._skytile == '' or self._filelist == '':
            raise Exception(
                'Set the region name and list of eventfiles first with the '
                'functions set_filelist and set_region.')
        if not os.path.exists(self._esass):
            raise Exception(f'File {self._esass} does not exist.')
        if not (os.path.exists(self._working_dir+'/src.reg')
                or not os.path.exists(self._working_dir+'/bkg.reg')):
            raise Exception(
                'Source and background extraction regions have to be defined '
                'before running the script.')
        for bin_e in self._energy_bins + self._energy_bins_hr:
            replacements = [['@source_name', self._src_name],
                            ['@main_name', self._working_dir],
                            ['@result_dir', '.'],
                            ['@region_code', self._skytile],
                            ['@sources_list', self._filelist],
                            ['@right_ascension', self._RA],
                            ['@declination', self._Dec],
                            ['@esass_location', self._esass],
                            ['@binning', self._LC_prebinning],
                            ['@emin', bin_e[0]],
                            ['@emax', bin_e[1]]]
            sh_file = self._working_dir_full + '/working/extract_lc.sh'
            sh_file = self._replace_in_sh(sh_file, replacements)
            process = subprocess.Popen(
                [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()  # Wait for process to complete.

            # iterate on the stdout line by line
            if not logname == '':
                filename = (f'{self._working_dir}/logfiles/lightcurves/'
                            f'{logname}_{bin_e[0]}keV_{bin_e[1]}'
                            'keV.log')
                logstate = setup_logfile(self._logger, filename)
                for line in process.stdout.readlines():
                    # to fix the weird b'something' format
                    self._logger.info(str(line)[2:-3] + '\n')
                self._logger.handlers = logstate
        self._LC_extracted = True
        # one month gap minimum to sort out possible short gaps due to
        # problems during observation
        self._find_obs_periods(60 * 60 * 24 * 30)
        self._eRASS_vs_epoch()

    def _find_obs_periods(self, gapsize):
        self._logger.info('Finding Observation Periods.')
        bin_e = self._energy_bins[0]
        with fits.open(f'./{self._src_name}_{self._skytile}_eROSITA_PATall_'
                       f'{self._LC_prebinning}s_{bin_e[0]}keV_{bin_e[1]}keV_'
                       f'020_LightCurve_00001.fits') as hdulist:
            time = hdulist[1].data.field('TIME').tolist()
            time.sort()
        self._obs_periods = []
        temp = []
        for i in range(len(time) - 1):
            if i == 0:
                # to definitely not lose any events
                temp = [time[i] - 2 * float(self._LC_prebinning)]
            if i == len(time) - 2:
                # to definitely not lose any events
                temp.append(time[i+1] + 2 * float(self._LC_prebinning))
                self._obs_periods.append(temp)
                continue
            if time[i + 1] - time[i] > gapsize:
                temp.append(time[i] + 2 * float(self._LC_prebinning))
                self._obs_periods.append(temp)
                temp = [time[i+1] - 2 * float(self._LC_prebinning)]
        self._obs_periods = np.array(
            self._obs_periods) / 3600. / 24. + self._mjdref
        # getting rid of spurious LC entries for obs_periods < 1 day
        index = []
        for i in range(len(self._obs_periods)):
            if self._obs_periods[i][1] - self._obs_periods[i][0] < 1:
                index.append(i)
        self._obs_periods = np.delete(self._obs_periods, index, axis=0)

    def _eRASS_vs_epoch(self):
        self._logger.info('Running eRASS vs. epoch.')
        self._period_names = []
        for i, period in enumerate(self._obs_periods):
            for date in self._ero_starttimes:
                if period[0] < date and period[1] > date:
                    self._create_epochs = True
        if self._create_epochs:
            for i in range(1, len(self._obs_periods) + 1):
                self._period_names.append(f'epoch{i}')
        else:
            for period in self._obs_periods:
                for j in range(len(self._ero_starttimes)):
                    if j == 0:
                        if period[1] < self._ero_starttimes[j]:
                            self._period_names.append(f'e{self._ownership}00')
                    if j == len(self._ero_starttimes) - 1:
                        if period[0] > self._ero_starttimes[j]:
                            self._period_names.append(
                                f'e{self._ownership}0{j + 1}')
                        continue
                    if (period[0] > self._ero_starttimes[j]
                            and period[1] < self._ero_starttimes[j + 1]):
                        self._period_names.append(
                            f'e{self._ownership}0{j + 1}')

    def plot_lc_full(self, fracexp='0.15', mincounts='10', mode='mincounts_ul',
                     show_eRASS=True, logname='lc_full_autosave.log',
                     time_axis='mjd', print_name=False, print_datetime=False,
                     label_style='serif', label_size=16, figsize=[16, 7],
                     colors=[], fileid='', toplab='', separate_TM=False,
                     vlines=[], ticknumber_y=5, ticknumber_x=8, E_bins=[],
                     lc_binning=-1, yscale='linear'):
        '''Function to create full lightcurve.

        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account
            for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted
            as an upper limit as well as minimum number of counts per
            bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or
            'mincounts_ul'. The default is 'ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines.
            The default is True.
        logname : str, optional
            Name of the logfile. The default is 'lc_full_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The
            default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is
            False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for
            matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 2.75].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib
            possible. For mode 'ul' and 'mincounts' the first entry is
            used, for mode 'mincounts_ul' the first entry sets color for
            'ul' part, and the second for 'mincounts' part. The default
            is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default
            is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with
            given color. The zorder entries need to be distinct negative
            integers < -2. The default is [].
        ticknumber_y : int, optional
            Sets the approximate number of tickmarks along the y axis.
            The default is 5.0.
        ticknumber_x : int, optional
            Sets the approximate number of tickmarks along the x axis.
            The default is 8.0.
        E_bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin E_min
            and E_max must be given in keV. The default is [[0.2, 8.0]]
        lc_binning : str or float, optional
            Sets initial lc binsize in seconds. The default is -1
            (meaning the current value is not changeds)
        yscale : str, optional
            Scale for yaxis of LC, either 'linear' or 'log'. The default
            is 'linear'.

        '''
        self._logger.info('Running plot_lc_full.')
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if (type(mincounts) != str and type(mincounts) != float
                and type(mincounts) != int):
            raise Exception('mincounts must be a string, float or int.')
        else:
            try:
                mincounts = float(mincounts)
            except ValueError:
                raise Exception('mincounts must be a number.')
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if mode != 'ul' and mode != 'mincounts' and mode != 'mincounts_ul':
                raise Exception(
                    'mode must be \'ul\', \'mincounts\' or \'mincounts_ul\'')
        if type(fracexp) != str and type(fracexp) != float:
            raise Exception('mincounts must be a string or float.')
        else:
            try:
                fracexp = float(fracexp)
            except ValueError:
                raise Exception('fracexp must be a number.')
        if type(show_eRASS) != bool:
            raise Exception('show_eRASS must be a bool.')
        if type(print_name) != bool:
            raise Exception('print_name must be a bool.')
        if type(print_datetime) != bool:
            raise Exception('print_datetime must be a bool.')
        if type(separate_TM) != bool:
            raise Exception('separate_TM must be a bool.')
        if type(time_axis) != str:
            raise Exception('time_axis must be a string.')
        else:
            if time_axis != 'mjd' and time_axis != 's':
                raise Exception('time_axis must be \'mjd\' or \'s\'')
        if type(label_style) != str:
            raise Exception('label_style must be a string.')
        if type(label_size) != float and type(label_size) != int:
            raise Exception('label_size must be a float or int.')
        if ((type(figsize) != list and type(figsize) != np.ndarray)
                or np.shape(figsize) != (2,)):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be an int.')
        if type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be an int.')
        if colors != []:
            if ((type(colors) != list and type(colors) != np.ndarray)
                    or (np.shape(colors) != (2,)
                        and np.shape(colors) != (1,))):
                raise Exception('colors must be (2,) or (1,) array-like.')
        if type(fileid) != str:
            raise Exception('fileid must be a string.')
        if type(toplab) != str:
            raise Exception('toplab must be a string.')
        if type(vlines) != list and type(vlines) != np.ndarray:
            raise Exception('vlines must be array-like')
        else:
            for line in vlines:
                if len(line) != 3:
                    raise Exception('Each line in vlines needs 3 entries.')
                if type(line[0]) != float and type(line[0]) != int:
                    raise Exception(
                        'The first entry in each line of vlines needs to be '
                        'the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception(
                        'The second entry in each line of vlines needs to be a'
                        ' matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
        if type(yscale) != str:
            raise Exception('mode must be a string.')
        else:
            if (yscale != 'linear'
                    and yscale != 'log'):
                raise Exception('yscale must be \'linear\' or \'log\'.')

        os.chdir(self._working_dir_full + '/working/')

        if lc_binning != -1:
            self.set_LC_binning(lc_binning=lc_binning)
        if np.array(E_bins).tolist() != []:
            self.set_Ebins(bins=E_bins)
        if not self._LC_extracted and not self._debugging:
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()
        logname_full = f'{self._working_dir}/logfiles/lightcurves/{logname}'
        logstate = setup_logfile(self._logger, logname_full)
        localtime = time.asctime(time.localtime(time.time()))

        self._logger.info(localtime + '\n')
        user = getpass.getuser()
        plt.rc('text', usetex=True)
        plt.rc('font', family=label_style, size=label_size)

        if separate_TM:
            TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            TM_list = [0]

        for bin_e in self._energy_bins:
            for TM in TM_list:
                if fileid == '':
                    pfile = (f'./{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                             f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}'
                             'keV_fullLC')
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                               f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                               'fullLC')
                else:
                    pfile = f'./{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_fullLC'
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_fullLC')
                if yscale != 'linear':
                    outfile += f'_{yscale}'
                replacements = [['@esass_location', self._esass],
                                ['@infile',
                                 f'./{self._src_name}_{self._skytile}_eROSITA_'
                                 f'PATall_1.0s_{bin_e[0]}keV_{bin_e[1]}keV_'
                                 f'{TM}20_LightCurve_00001.fits'],
                                ['@pfile', f'{pfile}.fits'],
                                ['@selection', f'FRACEXP>{fracexp}']]
                sh_file = self._working_dir_full + '/working/fselect_lc.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                fig1 = plt.figure(figsize=(figsize[0], figsize[1]))
                ax = fig1.add_subplot(111)
                ax.set_yscale(yscale)  # sets yscale to log if wanted

                self._logger.info(f'Now working on {pfile}.fits\n')
                hdulist = fits.open(f'{pfile}.fits')

                xflag = 0
                if time_axis == 'mjd':
                    xflag = 2
                elif time_axis == 's':
                    xflag = 1

                if colors == []:
                    colors = ['lightblue', 'black']

                pxmin, pxmax, pymin, pymax = 0, 0, 0, 0
                if mode == 'ul':
                    pxmin, pxmax, pymin, pymax = plot_lc_UL(
                        hdulist=hdulist, ax=ax, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1])
                elif mode == 'mincounts':
                    pxmin, pxmax, pymin, pymax = plot_lc_mincounts(
                        hdulist=hdulist, ax=ax, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1])
                elif mode == 'mincounts_ul':
                    pxmin1, pxmax1, pymin1, pymax1 = plot_lc_UL(
                        hdulist=hdulist, ax=ax, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[0])
                    pxmin2, pxmax2, pymin2, pymax2 = plot_lc_mincounts(
                        hdulist=hdulist, ax=ax, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1])
                    pxmin, pxmax, pymin, pymax = get_boundaries(
                        [[pxmin1, pxmax1, pymin1, pymax1],
                         [pxmin2, pxmax2, pymin2, pymax2]])

                format_axis(ax, pxmin, pxmax, pymin, pymax,
                            ticknumber_x, ticknumber_y)

                hdulist.close()

                # plot time in s from beginning (xflag=1) or in MJD
                if time_axis == 's':
                    ax.set_xlabel(r'Time (s)')  # , fontsize=12)
                elif time_axis == 'mjd':
                    ax.set_xlabel(r'MJD (days)')  # , fontsize=12)

                ax.set_ylabel(r'Count rate (cts/s)')  # , fontsize=12)

                if print_name:
                    # user name and time
                    ax.text(1.015, 0.0, user + ' - ' + localtime, rotation=90,
                            fontsize=8, verticalalignment='bottom',
                            horizontalalignment='left', transform=ax.transAxes)
                # eROSITA label
                    ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='left', transform=ax.transAxes)
                    ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='right',
                            transform=ax.transAxes)
                # label plot:
                    ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='center',
                            transform=ax.transAxes)

                for i in range(len(vlines)):
                    ax.vlines(vlines[i][0], pymin, pymax, colors=vlines[i]
                              [1], linestyle='dotted', zorder=vlines[i][2])
                if show_eRASS:
                    if time_axis == 'mjd':
                        ax.vlines(self._ero_starttimes, pymin, pymax,
                                  colors=['grey'], linestyle='dotted',
                                  zorder=-2)
                    elif time_axis == 's':
                        ax.vlines((np.array(self._ero_starttimes)
                                   - self._mjdref) * 3600 * 24, pymin, pymax,
                                  colors=['grey'], linestyle='dotted',
                                  zorder=-4)

                fig1.tight_layout()

                pltfile = outfile + ".pdf"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')
                pltfile = outfile + ".eps"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')
                pltfile = outfile + ".png"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')

        self._logger.handlers = logstate

    def plot_lc_broken(self, fracexp='0.15', mincounts='10',
                       mode='mincounts_ul', show_eRASS=True,
                       logname='lc_full_broken_autosave.log',
                       time_axis='mjd', print_name=False, print_datetime=False,
                       label_style='serif', label_size=16, figsize=[16, 7],
                       colors=[], fileid='', toplab='', separate_TM=False,
                       vlines=[], ticknumber_y=5, ticknumber_x=3, E_bins=[],
                       lc_binning=-1, d=12, tilt=45, diag_color="k",
                       short_time=True, fig_borders=[0.97, 0.1, 0.05, 0.98],
                       yscale='linear'):
        '''Function to create full lightcurve with gaps cut out.

        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account
            for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted
            as an upper limit as well as minimum number of counts per
            bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or
            'mincounts_ul'. The default is 'mincounts_ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines.
            The default is True.
        logname : str, optional
            Name of the logfile. The default is
            'lc_full_broken_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The
            default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is
            False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for
            matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 2.75].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib
            possible. For mode 'ul' and 'mincounts' the first entry is
            used, for mode 'mincounts_ul' the first entry sets color for
            'ul' part, and the second for 'mincounts' part. The default
            is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default
            is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with
            given color. The zorder entries need to be distinct negative
            integers < -2. The default is [].
        ticknumber_y : int, optional
            Sets the approximate number of tickmarks along the y axis.
            The default is 5.
        ticknumber_x : int, optional
            Sets the approximate number of tickmarks along the x axis in
            each section. The default is 3.
        E_bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin E_min
            and E_max must be given in keV. The default is [[0.2, 8.0]]
        lc_binning : str or float, optional
            Sets initial lc binsize in seconds. The default is -1
            (meaning the current value is not changeds)
        d : str, int or float, optional
            Size of gap markers in pt. The default is '12'.
        tilt : str, int or float, optional
            Tild of gap markers. The default is '45'.
        diag_color : str, optional
            Color of gap markers. The default is 'k'.
        short_time : bool, optional
            Shorten time stamps in x-axis by subtracting value of lowest
            enrty. The default is True.
        fig_borders : array-like (n,1), optional
            Sets the borders of the figure (top, bottom, left, right).
            The default is [0.97, 0.1, 0.05, 0.98].
        yscale : str, optional
            Scale for yaxis of LC, either 'linear' or 'log'. The default
            is 'linear'.
        '''
        self._logger.info('Running plot_lc_broken.')
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if (type(mincounts) != str and type(mincounts) != float
                and type(mincounts) != int):
            raise Exception('mincounts must be a string, float or int.')
        else:
            try:
                mincounts = float(mincounts)
            except ValueError:
                raise Exception('mincounts must be a number.')
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if (mode != 'ul' and mode != 'mincounts'
                    and mode != 'mincounts_ul'):
                raise Exception('mode must be \'ul\', \'mincounts\' or '
                                '\'mincounts_ul\'')
        if type(fracexp) != str and type(fracexp) != float:
            raise Exception('mincounts must be a string or float.')
        else:
            try:
                fracexp = float(fracexp)
            except ValueError:
                raise Exception('fracexp must be a number.')
        if type(show_eRASS) != bool:
            raise Exception('show_eRASS must be a bool.')
        if type(print_name) != bool:
            raise Exception('print_name must be a bool.')
        if type(print_datetime) != bool:
            raise Exception('print_datetime must be a bool.')
        if type(separate_TM) != bool:
            raise Exception('separate_TM must be a bool.')
        if type(time_axis) != str:
            raise Exception('time_axis must be a string.')
        else:
            if time_axis != 'mjd' and time_axis != 's':
                raise Exception('time_axis must be \'mjd\' or \'s\'')
        if type(label_style) != str:
            raise Exception('label_style must be a string.')
        if type(label_size) != float and type(label_size) != int:
            raise Exception('label_size must be a float or int.')
        if ((type(figsize) != list and type(figsize) != np.ndarray)
                or np.shape(figsize) != (2,)):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be an int.')
        if type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be an int.')
        if colors != []:
            if ((type(colors) != list and type(colors) != np.ndarray)
                    or (np.shape(colors) != (2,)
                        and np.shape(colors) != (1,))):
                raise Exception('colors must be (2,) or (1,) array-like.')
        if type(fileid) != str:
            raise Exception('fileid must be a string.')
        if type(toplab) != str:
            raise Exception('toplab must be a string.')
        if type(vlines) != list and type(vlines) != np.ndarray:
            raise Exception('vlines must be array-like')
        else:
            for line in vlines:
                if len(line) != 3:
                    raise Exception('Each line in vlines needs 3 entries.')
                if type(line[0]) != float and type(line[0]) != int:
                    raise Exception(
                        'The first entry in each line of vlines needs to be '
                        'the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception(
                        'The second entry in each line of vlines needs to be a'
                        ' matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
        if type(tilt) != float and type(tilt) != int and type(tilt) != str:
            raise Exception('tilt must be a float, str or int.')
        if type(d) != float and type(d) != str and type(d) != int:
            raise Exception('d must be a float or str.')
        if type(diag_color) != str:
            raise Exception('diag_color must be a str.')
        if type(fig_borders) != list and type(fig_borders) != np.ndarray:
            raise Exception('fig_borders must be array-like')
        else:
            if len(fig_borders) != 4:
                raise Exception('fig_borders needs exactly 4 entries')
            for entry in fig_borders:
                if type(entry) != float:
                    raise Exception(
                        'Entries in fig_borders need to be of type float.')
        if type(yscale) != str:
            raise Exception('mode must be a string.')
        else:
            if (yscale != 'linear'
                    and yscale != 'log'):
                raise Exception('yscale must be \'linear\' or \'log\'.')

        os.chdir(self._working_dir_full + '/working/')

        if lc_binning != -1:
            self.set_LC_binning(lc_binning=lc_binning)
        if np.array(E_bins).tolist() != []:
            self.set_Ebins(bins=E_bins)
        if not self._LC_extracted and not self._debugging:
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()
        logname_full = f'{self._working_dir}/logfiles/lightcurves/{logname}'
        logstate = setup_logfile(self._logger, logname_full)
        localtime = time.asctime(time.localtime(time.time()))

        self._logger.info(localtime + '\n')
        user = getpass.getuser()
        plt.rc('text', usetex=True)
        plt.rc('font', family=label_style, size=label_size)

        if separate_TM:
            TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            TM_list = [0]
        time_rel = 0
        pxmin = []
        pxmax = []
        pymin = 0
        pymax = 0
        xflag = 0

        for bin_e in self._energy_bins:
            for TM in TM_list:
                if fileid == '':
                    pfile = (f'./{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                             f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                             'brokenLC')
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                               f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                               'brokenLC')
                else:
                    pfile = f'./{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_brokenLC'
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_'
                               'brokenLC')
                if yscale != 'linear':
                    outfile += f'_{yscale}'
                replacements = [['@esass_location', self._esass],
                                ['@infile',
                                 f'./{self._src_name}_{self._skytile}_eROSITA_'
                                 f'PATall_1.0s_{bin_e[0]}keV_{bin_e[1]}keV_'
                                 f'{TM}20_LightCurve_00001.fits'],
                                ['@pfile', f'{pfile}.fits'],
                                ['@selection', f'FRACEXP>{fracexp}']]
                sh_file = self._working_dir_full + '/working/fselect_lc.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                fig1 = plt.figure(figsize=(figsize[0], figsize[1]))

                ncols, nrows = len(self._obs_periods), 1

                big_ax = fig1.add_subplot(111)
                big_ax.set_frame_on(False)
                big_ax.patch.set_facecolor("none")

                axs = []
                for _ in self._obs_periods:
                    ax = fig1.add_subplot(111)
                    axs.append(ax)

                self._logger.info(f'Now working on {pfile}.fits\n')
                hdulist = fits.open(f'{pfile}.fits')

                if time_axis == 'mjd':
                    xflag = 2
                elif time_axis == 's':
                    xflag = 1

                if colors == []:
                    colors = ['lightblue', 'black']

                if mode == 'ul':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_UL_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pymin = min(pymin)
                    pymax = max(pymax)
                elif mode == 'mincounts':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_mincounts_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pymin = min(pymin)
                    pymax = max(pymax)
                elif mode == 'mincounts_ul':
                    pxmin1, pxmax1, pymin1, pymax1, time_rel = plot_lc_UL_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[0], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pxmin2, pxmax2, pymin2, pymax2, _ = plot_lc_mincounts_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time, time_rel=time_rel)
                    pxmin, pxmax, pymin, pymax = get_boundaries_broken(
                        [[pxmin1, pxmax1, pymin1, pymax1],
                         [pxmin2, pxmax2, pymin2, pymax2]])

                hdulist.close()

                # plot time in s from beginning (xflag=1) or in MJD
                if short_time:
                    if time_axis == 's':
                        # , fontsize=12)
                        big_ax.set_xlabel(f'Time - {time_rel} (s)')
                    elif time_axis == 'mjd':
                        # , fontsize=12)
                        big_ax.set_xlabel(f'MJD  - {time_rel} (days)')
                else:
                    if time_axis == 's':
                        big_ax.set_xlabel(r'Time (s)')  # , fontsize=12)
                    elif time_axis == 'mjd':
                        big_ax.set_xlabel(r'MJD (days)')  # , fontsize=12)

                big_ax.set_ylabel(r'Count rate (cts/s)')  # , fontsize=12)

                if print_name:
                    # user name and time
                    big_ax.text(1.015, 0.0, f'{user} - {localtime}',
                                rotation=90, fontsize=8,
                                verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=big_ax.transAxes)
                # eROSITA label
                    big_ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=big_ax.transAxes)
                    big_ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='right',
                                transform=big_ax.transAxes)
                # label plot:
                    big_ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='center',
                                transform=big_ax.transAxes)

                for ax in axs:
                    ax.set_yscale(yscale)
                    for i in range(len(vlines)):
                        ax.vlines(vlines[i][0] - time_rel, pymin, pymax,
                                  colors=vlines[i][1], linestyle='dotted',
                                  zorder=vlines[i][2])
                    if show_eRASS:
                        if time_axis == 'mjd':
                            ax.vlines(self._ero_starttimes - time_rel, pymin, pymax,
                                      colors='grey', linestyle='dotted',
                                      zorder=-2)
                        elif time_axis == 's':
                            ax.vlines((np.array(self._ero_starttimes)
                                       - self._mjdref) * 3600 * 24 - time_rel,
                                      pymin, pymax, colors='grey',
                                      linestyle='dotted', zorder=-4)

                format_axis_broken_new(fig1, axs, pxmin, pxmax, pymin, pymax,
                                       ticknumber_x, ticknumber_y, ncols,
                                       nrows, d, tilt, diag_color, big_ax,
                                       yscale)

                fig1.set_tight_layout(True)
                fig1.set_tight_layout(False)
                wspace = 8.0 / figsize[0] * 0.05
                fig1.subplots_adjust(
                    wspace=wspace, top=fig_borders[0], bottom=fig_borders[1],
                    left=fig_borders[2], right=fig_borders[3])

                width_ratios = []
                height_ratios = [1]
                for i_ax in range(len(axs)):
                    width_ratios.append(pxmax[i_ax] - pxmin[i_ax])

                gs = gridspec.GridSpec(ncols=ncols,
                                       nrows=nrows,
                                       height_ratios=height_ratios,
                                       width_ratios=width_ratios)

                for i_ax, ax in enumerate(axs):
                    ax.set_position(gs[i_ax].get_position(fig1))

                self._width_ratios = width_ratios
                self._fig = fig1
                self._axes = axs
                self._big_ax = big_ax

                pltfile = outfile + ".pdf"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')
                pltfile = outfile + ".eps"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')
                pltfile = outfile + ".png"
                plt.savefig(pltfile)
                self._logger.info(f'{pltfile} created\n')

        self._logger.handlers = logstate

    def plot_lc_parts(self, fracexp='0.15', mincounts='10',
                      mode='mincounts_ul', show_eRASS=True,
                      logname='lc_parts_autosave.log', time_axis='mjd',
                      print_name=False, print_datetime=False,
                      label_style='serif', label_size=16, figsize=[8, 5.5],
                      colors=[], fileid='', toplab='', separate_TM=False,
                      vlines=[], ticknumber_y=5, ticknumber_x=8, eRASSi=[],
                      E_bins=[], lc_binning=-1):
        '''Function to create lightcurve of eRASSi/epochs.

        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account
            for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted
            as an upper limit as well as minimum number of counts per
            bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or
            'mincounts_ul'. The default is 'mincounts_ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines.
            The default is True.
        logname : str, optional
            Name of the logfile. The default is 'lc_parts_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The
            default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is
            False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for
            matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 5.5].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib
            possible. For mode 'ul' and 'mincounts' the first entry is
            used, for mode 'mincounts_ul' the first entry sets color for
            'ul' part, and the second for 'mincounts' part. The default
            is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default
            is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with
            given color. The zorder entries need to be distinct negative
            integers < -2. The default is [].
        ticknumber_y : int, optional
            Sets the approximate number of tickmarks along the y axis.
            The default is 5.0.
        ticknumber_x : int, optional
            Sets the approximate number of tickmarks along the x axis.
            The default is 8.0.
        eRASSi : array-like of ints
            List of from which eRASS eventfiles were used. The default
            is [].
        E_bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin E_min
            and E_max must be given in keV. The default is []
            (meaning the current value is not changed)
        lc_binning : str or float, optional
            Sets initial lc binsize in seconds. The default is -1
            (meaning the current value is not changeds)

        '''
        self._logger.info('Running plot_lc_parts.')
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if (type(mincounts) != str and type(mincounts) != float
                and type(mincounts) != int):
            raise Exception('mincounts must be a string, float or int.')
        else:
            try:
                mincounts = float(mincounts)
            except ValueError:
                raise Exception('mincounts must be a number.')
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if mode != 'ul' and mode != 'mincounts' and mode != 'mincounts_ul':
                raise Exception(
                    'mode must be \'ul\', \'mincounts\' or \'mincounts_ul\'')
        if type(fracexp) != str and type(fracexp) != float:
            raise Exception('mincounts must be a string or float.')
        else:
            try:
                fracexp = float(fracexp)
            except ValueError:
                raise Exception('fracexp must be a number.')
        if type(show_eRASS) != bool:
            raise Exception('show_eRASS must be a bool.')
        if type(print_name) != bool:
            raise Exception('print_name must be a bool.')
        if type(print_datetime) != bool:
            raise Exception('print_datetime must be a bool.')
        if type(separate_TM) != bool:
            raise Exception('separate_TM must be a bool.')
        if type(time_axis) != str:
            raise Exception('time_axis must be a string.')
        else:
            if time_axis != 'mjd' and time_axis != 's':
                raise Exception('time_axis must be \'mjd\' or \'s\'')
        if type(label_style) != str:
            raise Exception('label_style must be a string.')
        if type(label_size) != float and type(label_size) != int:
            raise Exception('label_size must be a float or int.')
        if ((type(figsize) != list and type(figsize) != np.ndarray)
                or np.shape(figsize) != (2,)):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be an int.')
        if type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be an int.')
        if colors != []:
            if ((type(colors) != list and type(colors) != np.ndarray)
                    or (np.shape(colors) != (2,)
                        and np.shape(colors) != (1,))):
                raise Exception('colors must be (2,) or (1,) array-like.')
        if type(fileid) != str:
            raise Exception('fileid must be a string.')
        if type(toplab) != str:
            raise Exception('toplab must be a string.')
        if type(vlines) != list and type(vlines) != np.ndarray:
            raise Exception('vlines must be array-like')
        else:
            for line in vlines:
                if len(line) != 3:
                    raise Exception('Each line in vlines needs 3 entries.')
                if type(line[0]) != float and type(line[0]) != int:
                    raise Exception(
                        'The first entry in each line of vlines needs to be '
                        'the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception(
                        'The second entry in each line of vlines needs to be a'
                        ' matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
        if type(eRASSi) != list and type(eRASSi) != np.ndarray:
            raise Exception('eRASSi must be array-like.')
        else:
            for i in range(len(eRASSi)):
                if type(eRASSi[i]) != int:
                    raise Exception('Entries in eRASSi must be integer.')
        if np.array(eRASSi) != np.sort(eRASSi):
            raise Exception('eRASSi must be sorted ascending.')
        if eRASSi != []:
            self._logger.warning('Use of eRASSi is currently not '
                                 'supported.')
        os.chdir(self._working_dir_full + '/working/')

        if lc_binning != -1:
            self.set_LC_binning(lc_binning=lc_binning)
        if np.array(E_bins).tolist() != []:
            self.set_Ebins(bins=E_bins)
        if not self._LC_extracted and not self._debugging:
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()

        xflag = 0
        pxmin, pxmax, pymin, pymax = 0, 0, 0, 0

        logname_full = f'{self._working_dir}/logfiles/lightcurves/{logname}'
        logstate = setup_logfile(self._logger, logname_full)

        for i, period in enumerate(self._obs_periods):
            naming = self._period_names[i]

            localtime = time.asctime(time.localtime(time.time()))

            self._logger.info(localtime + '\n')
            user = getpass.getuser()
            plt.rc('text', usetex=True)
            plt.rc('font', family=label_style, size=label_size)

            if separate_TM:
                TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
            else:
                TM_list = [0]

            for bin_e in self._energy_bins:
                for TM in TM_list:
                    if fileid == '':
                        pfile = (f'./{self._src_name}_{self._skytile}_LC_TM'
                                 f'{TM}20_fracexp{fracexp}_{bin_e[0]}keV_'
                                 f'{bin_e[1]}keV_LC_{naming}')
                        outfile = (
                            f'{self._working_dir}/results/lightcurves/'
                            f'{self._src_name}_{self._skytile}_LC_TM{TM}20_fra'
                            f'cexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_LC_'
                            f'{naming}')
                    else:
                        pfile = (f'./{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_'
                                 f'{naming}')
                        outfile = (f'{self._working_dir}/results/lightcurves/'
                                   f'{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_'
                                   f'{naming}')
                    selection = (f'FRACEXP>{fracexp} && TIME < '
                                 f'{(period[1] - self._mjdref) * 24 * 3600} &&'
                                 ' TIME > '
                                 f'{(period[0] - self._mjdref) * 24 * 3600}')
                    replacements = [['@esass_location', self._esass],
                                    ['@infile', f'./{self._src_name}_'
                                    f'{self._skytile}_eROSITA_PATall_1.0s_'
                                     f'{bin_e[0]}keV_{bin_e[1]}keV_{TM}20_'
                                     'LightCurve_00001.fits'],
                                    ['@pfile', f'{pfile}.fits'],
                                    ['@selection', selection]]
                    sh_file = self._working_dir_full + '/working/fselect_lc.sh'
                    sh_file = self._replace_in_sh(sh_file, replacements)
                    process = subprocess.Popen(
                        [sh_file], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
                    process.wait()  # Wait for process to complete.

                    fig1 = plt.figure(figsize=(figsize[0], figsize[1]))
                    ax = fig1.add_subplot(111)

                    self._logger.info(f'Now working on {pfile}.fits\n')
                    hdulist = fits.open(f'{pfile}.fits')

                    if time_axis == 'mjd':
                        xflag = 2
                    elif time_axis == 's':
                        xflag = 1

                    if colors == []:
                        colors = ['lightblue', 'black']

                    if mode == 'ul':
                        pxmin, pxmax, pymin, pymax = plot_lc_UL(
                            hdulist=hdulist, ax=ax, log=self._logger,
                            mjdref=self._mjdref, xflag=xflag,
                            mincounts=mincounts, color=colors[0])
                    elif mode == 'mincounts':
                        pxmin, pxmax, pymin, pymax = plot_lc_mincounts(
                            hdulist=hdulist, ax=ax, log=self._logger,
                            mjdref=self._mjdref, xflag=xflag,
                            mincounts=mincounts, color=colors[0])
                    elif mode == 'mincounts_ul':
                        pxmin1, pxmax1, pymin1, pymax1 = plot_lc_UL(
                            hdulist=hdulist, ax=ax, log=self._logger,
                            mjdref=self._mjdref, xflag=xflag,
                            mincounts=mincounts, color=colors[0])
                        pxmin2, pxmax2, pymin2, pymax2 = plot_lc_mincounts(
                            hdulist=hdulist, ax=ax, log=self._logger,
                            mjdref=self._mjdref, xflag=xflag,
                            mincounts=mincounts, color=colors[1])
                        pxmin, pxmax, pymin, pymax = get_boundaries(
                            [[pxmin1, pxmax1, pymin1, pymax1],
                             [pxmin2, pxmax2, pymin2, pymax2]])

                    # plot time in s from beginning (xflag=1) or in MJD
                    if time_axis == 's':
                        ax.set_xlabel(r'Time (s)')  # , fontsize=12)
                    elif time_axis == 'mjd':
                        ax.set_xlabel(r'MJD (days)')  # , fontsize=12)

                    ax.set_ylabel(r'Count rate (cts/s)')  # , fontsize=12)

                    if print_name:
                        # user name and time
                        ax.text(1.015, 0.0, user + ' - ' + localtime,
                                rotation=90, fontsize=8,
                                verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=ax.transAxes)
                    # eROSITA label
                        ax.text(0.0, 1.015, 'eROSITA', rotation=0,
                                fontsize=10, verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=ax.transAxes)
                        ax.text(1.0, 1.015, 'MPE', rotation=0,
                                fontsize=10, verticalalignment='bottom',
                                horizontalalignment='right',
                                transform=ax.transAxes)
                    # label plot:
                        ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='center',
                                transform=ax.transAxes)

                    for i in range(len(vlines)):
                        ax.vlines(vlines[i][0], pymin, pymax, colors=vlines[i]
                                  [1], linestyle='dotted', zorder=vlines[i][2])
                    if show_eRASS:
                        if time_axis == 'mjd':
                            ax.vlines(self._ero_starttimes, pymin, pymax,
                                      colors=['grey'], linestyle='dotted',
                                      zorder=-2)
                        elif time_axis == 's':
                            ax.vlines((np.array(self._ero_starttimes)
                                       - self._mjdref) * 3600 * 24, pymin, pymax,
                                      colors=['grey'], linestyle='dotted',
                                      zorder=-4)

                    format_axis(ax, pxmin, pxmax, pymin, pymax,
                                ticknumber_x, ticknumber_y)

                    fig1.tight_layout()

                    pltfile = outfile + ".pdf"
                    plt.savefig(pltfile)
                    self._logger.info(f'{pltfile} created\n')
                    pltfile = outfile + ".eps"
                    plt.savefig(pltfile)
                    self._logger.info(f'{pltfile} created\n')
                    pltfile = outfile + ".png"
                    plt.savefig(pltfile)
                    self._logger.info(f'{pltfile} created\n')

        self._logger.handlers = logstate

    def plot_lc_HR(self, fracexp='0.15', mincounts='10',
                   mode='mincounts_ul', show_eRASS=True,
                   logname='lc_full_hr_autosave.log',
                   time_axis='mjd', print_name=False, print_datetime=False,
                   label_style='serif', label_size=16, figsize=[16, 21],
                   colors=[], fileid='', toplab='', separate_TM=False,
                   vlines=[], ticknumber_y=5, ticknumber_x=3, E_bins=[],
                   lc_binning=-1, d=12, tilt=45, diag_color="k",
                   short_time=True, fig_borders=[0.97, 0.1, 0.05, 0.98]):
        '''Function to create full lightcurve with gaps cut out.

        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account
            for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted
            as an upper limit as well as minimum number of counts per
            bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or
            'mincounts_ul'. The default is 'mincounts_ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines.
            The default is True.
        logname : str, optional
            Name of the logfile. The default is
            'lc_full_broken_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The
            default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is
            False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for
            matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 2.75].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib
            possible. For mode 'ul' and 'mincounts' the first entry is
            used, for mode 'mincounts_ul' the first entry sets color for
            'ul' part, and the second for 'mincounts' part. The default
            is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default
            is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with
            given color. The zorder entries need to be distinct negative
            integers < -2. The default is [].
        ticknumber_y : int, optional
            Sets the approximate number of tickmarks along the y axis.
            The default is 5.
        ticknumber_x : int, optional
            Sets the approximate number of tickmarks along the x axis in
            each section. The default is 3.
        E_bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin E_min
            and E_max must be given in keV. The default is [[0.2, 2.0],
            [2.0, 8.0]]
        lc_binning : str or float, optional
            Sets initial lc binsize in seconds. The default is -1
            (meaning the current value is not changeds)
        d : str, int or float, optional
            Size of gap markers in pt. The default is '12'.
        tilt : str, int or float, optional
            Tild of gap markers. The default is '45'.
        diag_color : str, optional
            Color of gap markers. The default is 'k'.
        short_time : bool, optional
            Shorten time stamps in x-axis by subtracting value of lowest
            enrty. The default is True.
        fig_borders : array-like (n,1), optional
            Sets the borders of the figure (top, bottom, left, right).
            The default is [0.97, 0.1, 0.05, 0.98].
        '''
        self._logger.info('Running plot_lc_HR.')
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if (type(mincounts) != str and type(mincounts) != float
                and type(mincounts) != int):
            raise Exception('mincounts must be a string, float or int.')
        else:
            try:
                mincounts = float(mincounts)
            except ValueError:
                raise Exception('mincounts must be a number.')
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if (mode != 'ul' and mode != 'mincounts'
                    and mode != 'mincounts_ul'):
                raise Exception('mode must be \'ul\', \'mincounts\' or '
                                '\'mincounts_ul\'')
        if type(fracexp) != str and type(fracexp) != float:
            raise Exception('mincounts must be a string or float.')
        else:
            try:
                fracexp = float(fracexp)
            except ValueError:
                raise Exception('fracexp must be a number.')
        if type(show_eRASS) != bool:
            raise Exception('show_eRASS must be a bool.')
        if type(print_name) != bool:
            raise Exception('print_name must be a bool.')
        if type(print_datetime) != bool:
            raise Exception('print_datetime must be a bool.')
        if type(separate_TM) != bool:
            raise Exception('separate_TM must be a bool.')
        if type(time_axis) != str:
            raise Exception('time_axis must be a string.')
        else:
            if time_axis != 'mjd' and time_axis != 's':
                raise Exception('time_axis must be \'mjd\' or \'s\'')
        if type(label_style) != str:
            raise Exception('label_style must be a string.')
        if type(label_size) != float and type(label_size) != int:
            raise Exception('label_size must be a float or int.')
        if ((type(figsize) != list and type(figsize) != np.ndarray)
                or np.shape(figsize) != (2,)):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be an int.')
        if type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be an int.')
        if colors != []:
            if ((type(colors) != list and type(colors) != np.ndarray)
                    or (np.shape(colors) != (2,)
                        and np.shape(colors) != (1,))):
                raise Exception('colors must be (2,) or (1,) array-like.')
        if type(fileid) != str:
            raise Exception('fileid must be a string.')
        if type(toplab) != str:
            raise Exception('toplab must be a string.')
        if type(vlines) != list and type(vlines) != np.ndarray:
            raise Exception('vlines must be array-like')
        else:
            for line in vlines:
                if len(line) != 3:
                    raise Exception('Each line in vlines needs 3 entries.')
                if type(line[0]) != float and type(line[0]) != int:
                    raise Exception(
                        'The first entry in each line of vlines needs to be '
                        'the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception(
                        'The second entry in each line of vlines needs to be a'
                        ' matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
        if type(tilt) != float and type(tilt) != int and type(tilt) != str:
            raise Exception('tilt must be a float, str or int.')
        if type(d) != float and type(d) != str and type(d) != int:
            raise Exception('d must be a float or str.')
        if type(diag_color) != str:
            raise Exception('diag_color must be a str.')
        if type(fig_borders) != list and type(fig_borders) != np.ndarray:
            raise Exception('fig_borders must be array-like')
        else:
            if len(fig_borders) != 4:
                raise Exception('fig_borders needs exactly 4 entries')
            for entry in fig_borders:
                if type(entry) != float:
                    raise Exception(
                        'Entries in fig_borders need to be of type float.')

        os.chdir(self._working_dir_full + '/working/')

        if lc_binning != -1:
            self.set_LC_binning(lc_binning=lc_binning)
        if np.array(E_bins).tolist() != []:
            self.set_Ebins_HR(bins=E_bins)
        if not self._LC_extracted and not self._debugging:
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()
        logname_full = f'{self._working_dir}/logfiles/lightcurves/{logname}'
        logstate = setup_logfile(self._logger, logname_full)
        localtime = time.asctime(time.localtime(time.time()))

        self._logger.info(localtime + '\n')
        user = getpass.getuser()
        plt.rc('text', usetex=True)
        plt.rc('font', family=label_style, size=label_size)

        if separate_TM:
            TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            TM_list = [0]
        time_rel = 0
        pxmin = []
        pxmax = []
        pymin = 0
        pymax = 0
        xflag = 0
        axs_full = []
        width_ratios = []
        height_ratios = []

        for TM in TM_list:
            pfiles = []
            if fileid == '':
                outfile = (f'{self._working_dir}/results/lightcurves/'
                           f'{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                           f'fracexp{fracexp}_hrLC')
            else:
                outfile = (f'{self._working_dir}/results/lightcurves/'
                           f'{fileid}_hrLC')
            fig1 = plt.figure(figsize=(figsize[0], figsize[1]))
            for i, bin_e in enumerate(self._energy_bins_hr):
                if fileid == '':
                    pfile = (f'./{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                             f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                             'hrLC')
                else:
                    pfile = f'./{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_hrLC'
                replacements = [['@esass_location', self._esass],
                                ['@infile',
                                 f'./{self._src_name}_{self._skytile}_eROSITA_'
                                 f'PATall_1.0s_{bin_e[0]}keV_{bin_e[1]}keV_'
                                 f'{TM}20_LightCurve_00001.fits'],
                                ['@pfile', f'{pfile}.fits'],
                                ['@selection', f'FRACEXP>{fracexp}']]
                sh_file = self._working_dir_full + '/working/fselect_lc.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                ncols, nrows = len(self._obs_periods), 3

                big_ax = fig1.add_subplot(310+1+i)
                big_ax.set_frame_on(False)
                big_ax.patch.set_facecolor("none")

                axs = []
                for _ in self._obs_periods:
                    ax = fig1.add_subplot(111)
                    axs.append(ax)
                    axs_full.append(ax)

                self._logger.info(f'Now working on {pfile}.fits\n')
                hdulist = fits.open(f'{pfile}.fits')
                pfiles.append(pfile)

                if time_axis == 'mjd':
                    xflag = 2
                elif time_axis == 's':
                    xflag = 1

                if colors == []:
                    colors = ['lightblue', 'black']

                if mode == 'ul':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_UL_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pymin = min(pymin)
                    pymax = max(pymax)
                elif mode == 'mincounts':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_mincounts_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pymin = min(pymin)
                    pymax = max(pymax)
                elif mode == 'mincounts_ul':
                    pxmin1, pxmax1, pymin1, pymax1, time_rel = plot_lc_UL_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[0], obs_periods=self._obs_periods,
                        short_time=short_time)
                    pxmin2, pxmax2, pymin2, pymax2, _ = plot_lc_mincounts_broken_new(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time, time_rel=time_rel)
                    pxmin, pxmax, pymin, pymax = get_boundaries_broken(
                        [[pxmin1, pxmax1, pymin1, pymax1],
                         [pxmin2, pxmax2, pymin2, pymax2]])

                hdulist.close()

                big_ax.set_ylabel(r'Count rate (cts/s)')  # , fontsize=12)

                for ax in axs:
                    for k in range(len(vlines)):
                        ax.vlines(vlines[k][0] - time_rel, pymin, pymax,
                                  colors=vlines[k][1], linestyle='dotted',
                                  zorder=vlines[k][2])
                    if show_eRASS:
                        if time_axis == 'mjd':
                            ax.vlines(self._ero_starttimes - time_rel, pymin, pymax,
                                      colors='grey', linestyle='dotted',
                                      zorder=-2)
                        elif time_axis == 's':
                            ax.vlines((np.array(self._ero_starttimes)
                                       - self._mjdref) * 3600 * 24 - time_rel,
                                      pymin, pymax, colors='grey',
                                      linestyle='dotted', zorder=-4)

                yscale = 'linear' #if needed this needs to be programmed
                format_axis_broken_new(fig1, axs, pxmin, pxmax, pymin, pymax,
                                       ticknumber_x, ticknumber_y, ncols,
                                       nrows, d, tilt, diag_color, big_ax,
                                       yscale)

                if i == 0:
                    width_ratios = []
                    height_ratios = [1/3, 1/3, 1/3]
                    for i_ax in range(len(axs)):
                        width_ratios.append(pxmax[i_ax] - pxmin[i_ax])

            ncols, nrows = len(self._obs_periods), 3

            big_ax = fig1.add_subplot(313)
            big_ax.set_frame_on(False)
            big_ax.patch.set_facecolor("none")

            axs = []
            for _ in self._obs_periods:
                ax = fig1.add_subplot(111)
                axs.append(ax)
                axs_full.append(ax)

            self._logger.info(f'Now working on HR\n')
            hdulist1 = fits.open(f'{pfiles[0]}.fits')
            hdulist2 = fits.open(f'{pfiles[1]}.fits')

            if time_axis == 'mjd':
                xflag = 2
            elif time_axis == 's':
                xflag = 1

            if colors == []:
                colors = ['lightblue', 'black']

            if mode == 'ul':
                pxmin, pxmax, pymin, pymax, time_rel = plot_lc_UL_hr(
                    hdulist_1=hdulist1, hdulist_2=hdulist2, axs=axs, log=self._logger,
                    mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                    color=colors[1], obs_periods=self._obs_periods,
                    short_time=short_time)
                pymin = min(pymin)
                pymax = max(pymax)
            elif mode == 'mincounts':
                pxmin, pxmax, pymin, pymax, time_rel = plot_lc_mincounts_hr(
                    hdulist_1=hdulist1, hdulist_2=hdulist2, axs=axs,
                    log=self._logger, mjdref=self._mjdref, xflag=xflag,
                    mincounts=mincounts, color=colors[1],
                    obs_periods=self._obs_periods, short_time=short_time)
                pymin = min(pymin)
                pymax = max(pymax)
            elif mode == 'mincounts_ul':
                pxmin1, pxmax1, pymin1, pymax1, time_rel = plot_lc_UL_hr(
                    hdulist_1=hdulist1, hdulist_2=hdulist2, axs=axs,
                    log=self._logger, mjdref=self._mjdref, xflag=xflag,
                    mincounts=mincounts, color=colors[0],
                    obs_periods=self._obs_periods, short_time=short_time)
                # pxmin2, pxmax2, pymin2, pymax2, _ = plot_lc_mincounts_hr(
                #     hdulist_1=hdulist1, hdulist_2=hdulist2, axs=axs,
                #     log=self._logger, mjdref=self._mjdref, xflag=xflag,
                #     mincounts=mincounts, color=colors[1],
                #     obs_periods=self._obs_periods, short_time=short_time,
                #     time_rel=time_rel)
                # pxmin, pxmax, pymin, pymax = get_boundaries_broken(
                #     [[pxmin1, pxmax1, pymin1, pymax1],
                #         [pxmin2, pxmax2, pymin2, pymax2]])
                pymin = min(pymin1)
                pymax = max(pymax1)

            hdulist1.close()
            hdulist2.close()

            # plot time in s from beginning (xflag=1) or in MJD
            if short_time:
                if time_axis == 's':
                    # , fontsize=12)
                    big_ax.set_xlabel(f'Time - {time_rel} (s)')
                elif time_axis == 'mjd':
                    # , fontsize=12)
                    big_ax.set_xlabel(f'MJD  - {time_rel} (days)')
            else:
                if time_axis == 's':
                    big_ax.set_xlabel(r'Time (s)')  # , fontsize=12)
                elif time_axis == 'mjd':
                    big_ax.set_xlabel(r'MJD (days)')  # , fontsize=12)

            big_ax.set_ylabel(r'Hardness Ratio')  # , fontsize=12)

            if print_name:
                # user name and time
                big_ax.text(1.015, 0.0, f'{user} - {localtime}',
                            rotation=90, fontsize=8,
                            verticalalignment='bottom',
                            horizontalalignment='left',
                            transform=big_ax.transAxes)
            # eROSITA label
                big_ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='left',
                            transform=big_ax.transAxes)
                big_ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='right',
                            transform=big_ax.transAxes)
            # label plot:
                big_ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                            verticalalignment='bottom',
                            horizontalalignment='center',
                            transform=big_ax.transAxes)

            for ax in axs:
                for k in range(len(vlines)):
                    ax.vlines(vlines[k][0] - time_rel, pymin, pymax,
                              colors=vlines[k][1], linestyle='dotted',
                              zorder=vlines[k][2])
                if show_eRASS:
                    if time_axis == 'mjd':
                        ax.vlines(self._ero_starttimes - time_rel, pymin, pymax,
                                  colors='grey', linestyle='dotted',
                                  zorder=-2)
                    elif time_axis == 's':
                        ax.vlines((np.array(self._ero_starttimes)
                                   - self._mjdref) * 3600 * 24 - time_rel,
                                  pymin, pymax, colors='grey',
                                  linestyle='dotted', zorder=-4)

            yscale = 'linear' #this needs to be changed entirely anyways
            format_axis_broken_new(fig1, axs, pxmin, pxmax, pymin, pymax,
                                   ticknumber_x, ticknumber_y, ncols,
                                   nrows, d, tilt, diag_color, big_ax, yscale)

            gs = gridspec.GridSpec(ncols=ncols,
                                   nrows=nrows,
                                   height_ratios=height_ratios,
                                   width_ratios=width_ratios)

            fig1.set_tight_layout(True)
            fig1.set_tight_layout(False)
            wspace = 8.0 / figsize[0] * 0.05
            fig1.subplots_adjust(
                wspace=wspace, top=fig_borders[0], bottom=fig_borders[1],
                left=fig_borders[2], right=fig_borders[3])

            for i_ax, ax in enumerate(axs_full):
                ax.set_position(gs[i_ax].get_position(fig1))

            self.gs = gs

            self._width_ratios = width_ratios
            self._fig = fig1
            self._axes = axs_full
            self._big_ax = big_ax

            pltfile = outfile + ".pdf"
            plt.savefig(pltfile)
            self._logger.info(f'{pltfile} created\n')
            pltfile = outfile + ".eps"
            plt.savefig(pltfile)
            self._logger.info(f'{pltfile} created\n')
            pltfile = outfile + ".png"
            plt.savefig(pltfile)
            self._logger.info(f'{pltfile} created\n')

        self._logger.handlers = logstate

    def plot_lc_bayes_broken(self, fracexp='0.15', mincounts='10',
                             mode='mincounts_scan', show_eRASS=True,
                             logfile='', stan_model='',
                             alpha_bg=0.5, quantiles=[],
                             time_axis='mjd', print_name=False,
                             print_datetime=False, label_style='serif',
                             label_size=16, figsize=[16, 7],
                             colors=[], fileid='', toplab='',
                             separate_TM=False, vlines=[], ticknumber_y=5,
                             ticknumber_x=3, E_bins=[], lc_binning=-1, d=12,
                             tilt=45, diag_color="k", short_time=True,
                             fig_borders=[0.97, 0.1, 0.05, 0.98],
                             bbp0=0.01, bbmode='both', yscale='linear'):
        '''Function to create full lightcurve with bayesian estimates
        for source and background countrates with time gaps cut out.

        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account
            for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted
            as an upper limit as well as minimum number of counts per
            bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'scan', 'mincounts',
            'mincounts_scan', 'mincounts_bb' (mincounts and bayesian
            blocks overlay) or 'scan_bb' (scan with bayesian blocks
            overlay). The default is 'mincounts_scan'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines.
            The default is True.
        logfile : str, optional
            Name of the logfile. The default is ''.
        fexp_cut : float, optional
            Minimum value of fractional exposure for time bins to be
            used. The default is 0.15.
        stan_model : str, optional
            Name of the .stan file to use as a model when fitting. If
            the default is used a standard Poisson model is used. The
            default is ''.
        alpha_bg : float, optional
            alpha value for plotting the estimated background countrate.
            The default is 0.5.
        quantiles : array-like (3,) or (0,) float or int, optional
            Quantiles to plot as lower boundary, expected value and
            upper boundary. If the default is used, 1 sigma percentiles
            and the median are used. The default is [].
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The
            default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is
            False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for
            matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 2.75].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib
            possible. For mode 'ul' and 'mincounts' the first entry is
            used, for mode 'mincounts_ul' the first entry sets color for
            'ul' part, and the second for 'mincounts' part. The default
            is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default
            is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with
            given color. The zorder entries need to be distinct negative
            integers < -2. The default is [].
        ticknumber_y : int, optional
            Sets the approximate number of tickmarks along the y axis.
            The default is 5.
        ticknumber_x : int, optional
            Sets the approximate number of tickmarks along the x axis in
            each section. The default is 3.
        E_bins : array-like (n,2), optional
            Sets energy bins that should be analysed. For each bin E_min
            and E_max must be given in keV. The default is [[0.2, 8.0]]
        lc_binning : str or float, optional
            Sets initial lc binsize in seconds. The default is -1
            (meaning the current value is not changeds)
        d : str, int or float, optional
            Size of gap markers in pt. The default is '12'.
        tilt : str, int or float, optional
            Tild of gap markers. The default is '45'.
        diag_color : str, optional
            Color of gap markers. The default is 'k'.
        short_time : bool, optional
            Shorten time stamps in x-axis by subtracting value of lowest
            enrty. The default is True.
        fig_borders : array-like (n,1), optional
            Sets the borders of the figure (top, bottom, left, right).
            The default is [0.97, 0.1, 0.05, 0.98].
        bbp0 : float, optional
            Value used for p0 when running astropy's bayesian blocks.
            The default is '0.01'.
        bbmode : str, optional
            Mode which bayesian block contours to draw. Either 'sc'
            (source counts only), 'both' (one overlay for source and
            background counts each) or 'sum' (overlay for sc+bg counts).
            The default is 'both'.
        yscale : str, optional
            Scale for yaxis of LC, either 'linear' or 'log'. The default
            is 'linear'.
        '''
        self._logger.info(f'Running plot_lc_HR in mode {mode}.')
        if type(logfile) != str:
            raise Exception('logdir must be a string.')
        if (type(mincounts) != str and type(mincounts) != float
                and type(mincounts) != int):
            raise Exception('mincounts must be a string, float or int.')
        else:
            try:
                mincounts = float(mincounts)
            except ValueError:
                raise Exception('mincounts must be a number.')
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if (mode != 'scan'
                and mode != 'mincounts'
                and mode != 'mincounts_scan'
                and mode != 'mincounts_bb'
                    and mode != 'scan_bb'):
                raise Exception('mode must be \'scan\', \'mincounts\', '
                                '\'mincounts_scan\', \'mincounts_bb\' or'
                                '\'scan_bb\'.')
        if type(fracexp) != str and type(fracexp) != float:
            raise Exception('mincounts must be a string or float.')
        else:
            try:
                fracexp = float(fracexp)
            except ValueError:
                raise Exception('fracexp must be a number.')
        if type(show_eRASS) != bool:
            raise Exception('show_eRASS must be a bool.')
        if type(print_name) != bool:
            raise Exception('print_name must be a bool.')
        if type(print_datetime) != bool:
            raise Exception('print_datetime must be a bool.')
        if type(separate_TM) != bool:
            raise Exception('separate_TM must be a bool.')
        if type(stan_model) != str:
            raise Exception('stan_model must be a string.')
        if type(time_axis) != str:
            raise Exception('time_axis must be a string.')
        else:
            if time_axis != 'mjd' and time_axis != 's':
                raise Exception('time_axis must be \'mjd\' or \'s\'')
        if type(label_style) != str:
            raise Exception('label_style must be a string.')
        if type(label_size) != float and type(label_size) != int:
            raise Exception('label_size must be a float or int.')
        if type(alpha_bg) != float:
            raise Exception('alpha_bg must be a float.')
        if ((type(figsize) != list and type(figsize) != np.ndarray)
                or np.shape(figsize) != (2,)):
            raise Exception('figsize must be (2,) array-like.')
        if ((type(quantiles) != list and type(quantiles) != np.ndarray)
                or (np.shape(quantiles) != (3,) and
                    np.shape(quantiles) != (0,))):
            raise Exception('quantiles must be (3,) or (0,) array-like.')
        if type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be an int.')
        if type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be an int.')
        if colors != []:
            if ((type(colors) != list and type(colors) != np.ndarray)
                    or (np.shape(colors) != (2,)
                        and np.shape(colors) != (1,))):
                raise Exception('colors must be (2,) or (1,) array-like.')
        if type(fileid) != str:
            raise Exception('fileid must be a string.')
        if type(toplab) != str:
            raise Exception('toplab must be a string.')
        if type(vlines) != list and type(vlines) != np.ndarray:
            raise Exception('vlines must be array-like')
        else:
            for line in vlines:
                if len(line) != 3:
                    raise Exception('Each line in vlines needs 3 entries.')
                if type(line[0]) != float and type(line[0]) != int:
                    raise Exception(
                        'The first entry in each line of vlines needs to be '
                        'the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception(
                        'The second entry in each line of vlines needs to be a'
                        ' matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception(
                        'The third entry in each line of vlines needs to be a '
                        'negative integer < -2.')
        if type(tilt) != float and type(tilt) != int and type(tilt) != str:
            raise Exception('tilt must be a float, str or int.')
        if type(d) != float and type(d) != str and type(d) != int:
            raise Exception('d must be a float or str.')
        if type(diag_color) != str:
            raise Exception('diag_color must be a str.')
        if type(fig_borders) != list and type(fig_borders) != np.ndarray:
            raise Exception('fig_borders must be array-like')
        else:
            if len(fig_borders) != 4:
                raise Exception('fig_borders needs exactly 4 entries')
            for entry in fig_borders:
                if type(entry) != float:
                    raise Exception(
                        'Entries in fig_borders need to be of type float.')
        if type(bbp0) != float:
            raise Exception('tilt must be a float.')
        if type(bbmode) != str:
            raise Exception('bbmode must be a str.')
        else:
            if (bbmode != 'sc'
                and bbmode != 'both'
                    and bbmode != 'sum'):
                raise Exception('bbmode must be \'sc\', \'both\' or \'sum\'')
        if type(yscale) != str:
            raise Exception('mode must be a string.')
        else:
            if (yscale != 'linear'
                    and yscale != 'log'):
                raise Exception('yscale must be \'linear\' or \'log\'.')

        os.chdir(self._working_dir_full + '/working/')

        if stan_model == '':
            stan_model = self._working_dir_full + '/working/lc_model_simu.stan'
        if np.shape(quantiles) == (0,):
            quantiles = scipy.stats.norm().cdf([-1, 0, 1]) * 100

        if lc_binning != -1:
            self.set_LC_binning(lc_binning=lc_binning)
        if np.array(E_bins).tolist() != []:
            self.set_Ebins(bins=E_bins)
        if not self._LC_extracted and not self._debugging:
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()
        localtime = time.asctime(time.localtime(time.time()))

        user = getpass.getuser()
        plt.rc('text', usetex=True)
        plt.rc('font', family=label_style, size=label_size)

        if separate_TM:
            TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            TM_list = [0]
        time_rel = 0
        pxmin = []
        pxmax = []
        pymin = 0
        pymax = 0
        xflag = 0

        for bin_e in self._energy_bins:
            if logfile == '':
                logfile = f'LC_{bin_e[0]}keV_{bin_e[1]}keV_fexp{fracexp}.log'
            logfile = self._working_dir_full + '/logfiles/lightcurves/' + logfile
            if os.path.exists(logfile):
                os.remove(logfile)
            logstate = setup_logfile(self._logger, logfile)
            for TM in TM_list:
                if fileid == '':
                    pfile = (f'./{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                             f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                             f'brokenLC_bayes_{mode}')
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{self._src_name}_{self._skytile}_LC_TM{TM}20_'
                               f'fracexp{fracexp}_{bin_e[0]}keV_{bin_e[1]}keV_'
                               f'brokenLC_bayes_{mode}')
                else:
                    pfile = (f'./{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_'
                             f'brokenLC_bayes_{mode}')
                    outfile = (f'{self._working_dir}/results/lightcurves/'
                               f'{fileid}_{bin_e[0]}keV_{bin_e[1]}keV_'
                               f'brokenLC_bayes_{mode}')
                if yscale != 'linear':
                    outfile += f'_{yscale}'
                replacements = [['@esass_location', self._esass],
                                ['@infile',
                                f'./{self._src_name}_{self._skytile}_eROSITA_'
                                 f'PATall_1.0s_{bin_e[0]}keV_{bin_e[1]}keV_'
                                 f'{TM}20_LightCurve_00001.fits'],
                                ['@pfile', f'{pfile}.fits'],
                                ['@selection', f'FRACEXP>{fracexp}']]
                sh_file = self._working_dir_full + '/working/fselect_lc.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                fig1 = plt.figure(figsize=(figsize[0], figsize[1]))

                ncols, nrows = len(self._obs_periods), 1

                big_ax = fig1.add_subplot(111)
                big_ax.set_frame_on(False)
                big_ax.patch.set_facecolor("none")

                axs = []
                for _ in self._obs_periods:
                    ax = fig1.add_subplot(111)
                    axs.append(ax)

                hdulist = fits.open(f'{pfile}.fits')

                if time_axis == 'mjd':
                    xflag = 2
                elif time_axis == 's':
                    xflag = 1

                if colors == []:
                    colors = ['lightblue', 'black']

                if mode == 'scan':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_eROday_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, color=colors[1],
                        obs_periods=self._obs_periods, short_time=short_time,
                        stan_model=stan_model, quantiles=quantiles,
                        time_rel=time_rel, fexp_cut=float(fracexp),
                        alpha_bg=alpha_bg, yscale=yscale)
                elif mode == 'mincounts':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_mincounts_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time, stan_model=stan_model,
                        quantiles=quantiles, time_rel=time_rel,
                        fexp_cut=float(fracexp), alpha_bg=alpha_bg,
                        yscale=yscale)
                elif mode == 'mincounts_scan':
                    pxmin1, pxmax1, pymin1, pymax1, time_rel = plot_lc_eROday_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, color=colors[0],
                        obs_periods=self._obs_periods, short_time=short_time,
                        stan_model=stan_model, quantiles=quantiles,
                        time_rel=time_rel, fexp_cut=float(fracexp),
                        alpha_bg=alpha_bg, yscale=yscale)
                    pxmin2, pxmax2, pymin2, pymax2, time_rel = plot_lc_mincounts_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time, stan_model=stan_model,
                        quantiles=quantiles, time_rel=time_rel,
                        fexp_cut=float(fracexp), alpha_bg=alpha_bg,
                        yscale=yscale)
                    pymin = min([pymin1, pymin2])
                    pymax = max([pymax1, pymax2])
                    pxmin = np.min([pxmin1, pxmin2], axis=0)
                    pxmax = np.max([pxmax1, pxmax2], axis=0)
                elif mode == 'mincounts_bb':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_mincounts_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, mincounts=mincounts,
                        color=colors[1], obs_periods=self._obs_periods,
                        short_time=short_time, stan_model=stan_model,
                        quantiles=quantiles, time_rel=time_rel,
                        fexp_cut=float(fracexp), alpha_bg=alpha_bg,
                        bblocks=True, bbmode=bbmode, bbp0=bbp0, yscale=yscale)
                elif mode == 'scan_bb':
                    pxmin, pxmax, pymin, pymax, time_rel = plot_lc_eROday_broken_bayes(
                        hdulist=hdulist, axs=axs, log=self._logger,
                        mjdref=self._mjdref, xflag=xflag, color=colors[1],
                        obs_periods=self._obs_periods, short_time=short_time,
                        stan_model=stan_model, quantiles=quantiles,
                        time_rel=time_rel, fexp_cut=float(fracexp),
                        alpha_bg=alpha_bg, bblocks=True, bbmode=bbmode,
                        bbp0=bbp0, yscale=yscale)

                hdulist.close()

                # plot time in s from beginning (xflag=1) or in MJD
                if short_time:
                    if time_axis == 's':
                        # , fontsize=12)
                        big_ax.set_xlabel(f'Time - {time_rel} (s)')
                    elif time_axis == 'mjd':
                        # , fontsize=12)
                        big_ax.set_xlabel(f'MJD  - {time_rel} (days)')
                else:
                    if time_axis == 's':
                        big_ax.set_xlabel(r'Time (s)')  # , fontsize=12)
                    elif time_axis == 'mjd':
                        big_ax.set_xlabel(r'MJD (days)')  # , fontsize=12)

                big_ax.set_ylabel(r'Count rate (cts/s)')  # , fontsize=12)

                if print_name:
                    # user name and time
                    big_ax.text(1.015, 0.0, f'{user} - {localtime}',
                                rotation=90, fontsize=8,
                                verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=big_ax.transAxes)
                # eROSITA label
                    big_ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='left',
                                transform=big_ax.transAxes)
                    big_ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='right',
                                transform=big_ax.transAxes)
                # label plot:
                    big_ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                                verticalalignment='bottom',
                                horizontalalignment='center',
                                transform=big_ax.transAxes)

                for ax in axs:
                    ax.set_yscale(yscale)  # set yscale to log if wanted
                    for i in range(len(vlines)):
                        ax.vlines(vlines[i][0] - time_rel, pymin, pymax,
                                  colors=vlines[i][1], linestyle='dotted',
                                  zorder=vlines[i][2])
                    if show_eRASS:
                        if time_axis == 'mjd':
                            ax.vlines(self._ero_starttimes - time_rel, pymin,
                                      pymax, colors='grey', linestyle='dotted',
                                      zorder=-2)
                        elif time_axis == 's':
                            ax.vlines((np.array(self._ero_starttimes)
                                       - self._mjdref) * 3600 * 24 - time_rel,
                                      pymin, pymax, colors='grey',
                                      linestyle='dotted', zorder=-4)

                format_axis_broken_new(fig1, axs, pxmin, pxmax, pymin, pymax,
                                       ticknumber_x, ticknumber_y, ncols,
                                       nrows, d, tilt, diag_color, big_ax,
                                       yscale)

                fig1.set_tight_layout(True)
                fig1.set_tight_layout(False)
                wspace = 8.0 / figsize[0] * 0.05
                fig1.subplots_adjust(
                    wspace=wspace, top=fig_borders[0], bottom=fig_borders[1],
                    left=fig_borders[2], right=fig_borders[3])

                width_ratios = []
                height_ratios = [1]
                for i_ax in range(len(axs)):
                    width_ratios.append(pxmax[i_ax] - pxmin[i_ax])

                gs = gridspec.GridSpec(ncols=ncols,
                                       nrows=nrows,
                                       height_ratios=height_ratios,
                                       width_ratios=width_ratios)

                for i_ax, ax in enumerate(axs):
                    ax.set_position(gs[i_ax].get_position(fig1))

                self._width_ratios = width_ratios
                self._fig = fig1
                self._axes = axs
                self._big_ax = big_ax

                pltfile = outfile + ".pdf"
                plt.savefig(pltfile)
                pltfile = outfile + ".eps"
                plt.savefig(pltfile)
                pltfile = outfile + ".png"
                plt.savefig(pltfile)

            self._logger.handlers = logstate

    def _extract_spectrum(self, logname, mode, filelist, epoch):
        self._logger.debug(f'Extracting spectrum for {filelist}.')
        replacements = [['@source_name', self._src_name],
                        ['@main_name', self._working_dir],
                        ['@result_dir', '.'],
                        ['@region_code', self._skytile],
                        ['@sources_list', filelist],
                        ['@right_ascension', self._RA],
                        ['@declination', self._Dec],
                        ['@esass_location', self._esass],
                        ['@group', self._grouping],
                        ['@mode', mode],
                        ['@epoch', epoch]]
        sh_file = self._working_dir_full + '/working/extract_spectrum.sh'
        sh_file = self._replace_in_sh(sh_file, replacements)
        process = subprocess.Popen(
            [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()  # Wait for process to complete.

        # iterate on the stdout line by line
        logstate = setup_logfile(self._logger, logname)
        for line in process.stdout.readlines():
            # to fix weird b'something' format
            self._logger.debug(str(line)[2:-3] + '\n')
        self._logger.handlers = logstate

    def plot_spectra(self, mode='all', log_prefix='spectrum',
                     log_suffix='.log', Z=-1, distance=-1,
                     skip_varabs=False, NH=-1., rebin=True,
                     rebin_params=[3, 10], rescale=False, rescale_F=[1e-6, 1.],
                     rescale_chi=[-5., 5.], model_file='', save_settings=False,
                     grouping=-1, fit_statistic='cstat', colors=[], markers=[],
                     title='', varabs_starting_pars=[1., 6 * 10**-2],
                     separate=False, plot_command=["ldata", "delchi"],
                     return_array=False, conf_contours=False, abund='wilm',
                     skip_eRASS=[]):
        '''Fit and plot spectrum in prefered mode.

        Parameters
        ----------
        mode : str, optional
            Sets how data should be arranged for analysis.
            'simultaneous' fits all eRASS spectra simultaneously
            (powerlaw index and absorption parameters are tied,
            normalisation  is free), 'merged' creates a single spectrum
            of all existing data, 'individual' fits each eRASS/epoch
            individually with no parameters tied, 'all' does each of the
            possibilities mentioned. The default is 'all'.
        log_prefix : str, optional
            Sets a prefix for logfiles. The default is 'spectrum'.
        log_suffix : str, optional
            Sets a suffix for logfiles (including data type extension).
            The default is 'autosafe.log'.
        Z : float or str, optional
            Sets the metalicity to use for the varabs model
            (e.g. '0.49' for LMC sources). The default is -1 which uses
            the previously set metalicity.
        distance : float or str, optional
            Sets the distance to the source in kpc to calculate
            luminosities. The default is -1 which uses the previously
            set value (default saved is 50).
        skip_varabs : bool, optional
            Skips model fitting of models using varabs
            (can be used for low count data to reduce computation time
            since the varabs fit usually does not converge then). The
            default is False.
        NH : float or str, optional
            Sets the MW absorption modeled by Tbabs in units of cm^-2.
            The default is -1 to use
            preset NH.
        rebin : bool, optional
            Rebin energy bins during plotting (does not influence fit).
            The default is True.
        rebin_params :  int array-like (2,), optional
            Sets rebin parameters for plotting. The default is [3, 10].
        rescale : bool, optional
            Rescale y axis during plotting (does not influence fit). The
            default is False.
        rescale_F :  float array-like (2,), optional
            Sets rescale parameters for top plot during plotting. The
            default is [1e-6, 1.].
        rescale_chi :  float array-like (2,), optional
            Sets rescale parameters for bottom plot during plotting. The
            default is [-5., 5.].
        model_file : str, optional
            Loads a *.xcm setting file with model and plotting settings.
            When used no additional settings are done, which means the
            file must always include both model and plotting settings.
            The default is ''.
        save_settings : bool, optional
            If True a setting file (*.xcm) with all plot and model
            settings used is saved for each fit done with the
            corresponding names. The default is False.
        grouping : int or str, optional
            Sets the minimum number of events per energy bin during
            spectrum extraction. The default is -1 which uses the
            previously set value (default saved is 1).
        fit_statistic : str, optional
            Sets the fit-statistic used. Possible entries are 'cstat'
            and 'chi'. If set to 'chi' a grouping parameter <20 will be
            overwritten to 20. Future releases will allow 'bxa' and
            '3ml' as input. The default is 'cstat'.
        colors : array-like (n,), optional
            If set, overwrites default colors during plotting. n has to
            be equal to the number of datasets used (number of
            eRASS/epochs). The 0th component will be used for merged and
            individual spectra. Entries have to be integers. The default
            is [].
        markers : array-like (n,), optional
            If set, overwrites default markers during plotting. n has to
            be equal to the number of datasets used (number of
            eRASS/epochs). The 0th component will be used for merged and
            individual spectra. Entries have to be integers. The default
            is [].
        title : str, optional
            Set title of spectrum plots. The default is ''.
        varabs_starting_pars : array-like (2,), optional
            Sets starting value of powerlaw index and absorption for
            models including varabs. The default is [1, absorption]
            (absorption parameter set before or 6 * 10 ** -2 if left to
            default).
        separate : bool, optional
            If True creates a plot for each epoch/eRASS separately
            instead of a combined one.
        plot_command : array-like (2,), optional
            Sets the two panels to be plot (must be available in xspec).
            The default is ["ldata", "delchi"]. Example alternative
            ["ldata", "rat"]
        return_array : bool, optional
            If True returns a numpy array of model and data points. Not
            available yet, will be included in future releases. The
            default is False.
        conf_contours : bool, optional
            If True plots confidence contours of model parameters. Not
            available yet, will be included in future releases. The
            default is False.
        abund : str, optional
            Sets the abundance table used. Default is 'wilm'
        latest_eRASS : int or str, optional
            Sets the number of latest eRASS in use. The default is 5.

        '''
        self._logger.info('Running plot_spectra.')
        # Checking if input is sensible
        if type(mode) != str:
            raise Exception('mode must be a string.')
        else:
            if (mode != 'all' and mode != 'individual' and mode != 'merged'
                    and mode != 'simultaneous'):
                raise Exception(
                    'mode must be \'individual\', \'merged\', \'simultaneous'
                    '\' or \'all\'.')
        if type(log_prefix) != str:
            raise Exception('log_prefix must be a string.')
        if type(log_suffix) != str:
            raise Exception('log_suffix must be a string.')
        if type(skip_varabs) != bool:
            raise Exception('skip_varabs must be a bool.')
        if type(rebin) != bool:
            raise Exception('rebin must be a bool.')
        if type(rescale) != bool:
            raise Exception('rescale must be a bool.')
        if type(save_settings) != bool:
            raise Exception('save_settings must be a bool.')
        if type(return_array) != bool:
            raise Exception('return_array must be a bool.')
        if type(conf_contours) != bool:
            raise Exception('conf_contours must be a bool.')
        if type(NH) != str and type(NH) != float:
            raise Exception('NH must be a string or float.')
        else:
            try:
                NH = float(NH)
                if NH <= 0 and not NH == -1.:
                    raise Exception('NH must be > 0.')
            except ValueError:
                raise Exception('NH must be a number.')
        if NH != -1.:
            self.set_NH(NH=NH)
        if type(rebin_params) != list and type(rebin_params) != np.ndarray:
            raise Exception('rebin_params must be array-like')
        else:
            for entry in rebin_params:
                if type(entry) != int:
                    raise Exception(
                        'The entries of rebin_params need to be given as int.')
                if entry <= 0:
                    raise Exception(
                        'The entries of rebin_params need to be >0.')
            if len(rebin_params) != 2:
                raise Exception('rebin_params needs exactly 2 entries.')
        if type(rescale_F) != list and type(rescale_F) != np.ndarray:
            raise Exception('rescale_F must be array-like')
        else:
            for entry in rescale_F:
                if type(entry) != float:
                    raise Exception(
                        'The entries of rescale_F need to be given as float.')
                if entry <= 0:
                    raise Exception('The entries of rescale_F need to be >0.')
            if len(rescale_F) != 2:
                raise Exception('rescale_F needs exactly 2 entries.')
        if type(rescale_chi) != list and type(rescale_chi) != np.ndarray:
            raise Exception('rescale_chi must be array-like')
        else:
            for entry in rescale_chi:
                if type(entry) != float:
                    raise Exception(
                        'The entries of rescale_chi need to be given as float.')
            if len(rescale_chi) != 2:
                raise Exception('rescale_chi needs exactly 2 entries.')
        if type(model_file) != str:
            raise Exception('model_file must be string.')
        if not os.path.exists(self._working_dir + '/' + model_file):
            raise Exception(f'File {model_file} does not exist.')
        if Z != -1:
            self.set_metallicity(Z=Z)
        if distance != -1:
            self.set_distance(distance)
        if grouping != -1:
            self.set_grouping(grouping)
        if type(fit_statistic) != str:
            raise Exception('fit_statistic must be a string.')
        else:
            if fit_statistic != 'chi' and fit_statistic != 'cstat':
                raise Exception(
                    'mode must be \'chi\' or \'cstat\'. (\'3ml\' and \'bxa\' '
                    'not supported yet.)')
            if fit_statistic == 'chi':
                if self._grouping < 20:
                    self._grouping = 20
                    self._logger.warning('Grouping set to 20 due to '
                                         'chi fit-statistic usage.')
        if type(colors) != list and type(colors) != np.ndarray:
            raise Exception('colors must be array-like')
        else:
            for entry in colors:
                if type(entry) != int:
                    raise Exception(
                        'The entries of colors need to be given as int.')
            if (colors != []
                    and len(colors) < len(self._filelist.split(sep=' '))):
                raise Exception(
                    'If colors are set, they need to be set for every eRASS '
                    'used.')
        if type(markers) != list and type(markers) != np.ndarray:
            raise Exception('markers must be array-like')
        else:
            for entry in markers:
                if type(entry) != int:
                    raise Exception(
                        'The entries of markers need to be given as int.')
            if (markers != []
                    and len(markers) < len(self._filelist.split(sep=' '))):
                raise Exception(
                    'If markers are set, they need to be set for every eRASS '
                    'used.')
        if type(title) != str:
            raise Exception('title must be a string.')
        if (type(varabs_starting_pars) != list
                and type(varabs_starting_pars) != np.ndarray):
            raise Exception('varabs_starting_pars must be array-like')
        else:
            for entry in varabs_starting_pars:
                if type(entry) != float:
                    raise Exception(
                        'The entries of varabs_starting_pars need to be given '
                        'as float.')
                if entry <= 0:
                    raise Exception(
                        'The entries of varabs_starting_pars need to be >0.')
            if len(varabs_starting_pars) != 2:
                raise Exception(
                    'varabs_starting_pars needs exactly 2 entries.')
        if type(plot_command) != list and type(plot_command) != np.ndarray:
            raise Exception('plot_command must be array-like')
        else:
            for entry in plot_command:
                if type(entry) != str:
                    raise Exception(
                        'The entries of plot_command need to be given as str.')
        if type(abund) != str:
            raise Exception('abund must be a string.')

        if type(skip_eRASS) != list and type(skip_eRASS) != np.ndarray:
            raise Exception('skip_eRASS must be array-like')
        else:
            for entry in skip_eRASS:
                if type(entry) != int:
                    raise Exception(
                        'The entries of skip_eRASS need to be given as int.')
                if entry <= 0:
                    raise Exception(
                        'The entries of skip_eRASS need to be >0.')
        skip_eRASS = np.array(skip_eRASS)

        # For RMF and ARF files to work as intended
        os.chdir(self._working_dir_full + '/working/')
        # Prerequisites
        if self._skytile == '' or self._filelist == '':
            raise Exception(
                'Set the region name and list of eventfiles first with the '
                'functions set_filelist and set_region.')
        if not self._LC_extracted and not self._debugging:  # for debugging
            self._extract_lc()
        if self._debugging:
            self._LC_extracted = True
            self._find_obs_periods(60 * 60 * 24 * 30)
            self._eRASS_vs_epoch()
        if not self._NH_set:
            self._logger.warning('Source specific NH not set yet. '
                                 'Consider rerunning the script with a value '
                                 'given.')

        table_name = self._working_dir + '/results/spectra/' + \
            log_prefix  # not sure if this works the inteded way
        log_prefix = self._working_dir + '/logfiles/spectra/' + log_prefix
        if mode == 'all':
            self._plot_spectra_simultaneous(table_name, log_prefix,
                                            skip_varabs, self._NH/(1e22),
                                            separate, rebin, rebin_params,
                                            rescale_F, rescale_chi, abund,
                                            skip_eRASS, varabs_starting_pars,
                                            plot_command, model_file, title,
                                            save_settings, log_suffix, colors,
                                            markers, fit_statistic)
            self._plot_spectra_merged(log_prefix, skip_eRASS, table_name,
                                      skip_varabs, self._NH/(1e22), rebin,
                                      rebin_params, rescale_F, rescale_chi,
                                      abund, varabs_starting_pars,
                                      plot_command, model_file, title,
                                      save_settings, log_suffix, colors,
                                      markers, fit_statistic)
            self._plot_spectra_individual(table_name, log_prefix, skip_eRASS,
                                          skip_varabs, self._NH/(1e22), rebin,
                                          rebin_params, rescale_F, rescale_chi,
                                          abund, varabs_starting_pars,
                                          plot_command, model_file, title,
                                          save_settings, log_suffix, colors,
                                          markers, fit_statistic)
        elif mode == 'individual':
            self._plot_spectra_individual(table_name, log_prefix, skip_eRASS,
                                          skip_varabs, self._NH/(1e22), rebin,
                                          rebin_params, rescale_F, rescale_chi,
                                          abund, varabs_starting_pars,
                                          plot_command, model_file, title,
                                          save_settings, log_suffix, colors,
                                          markers, fit_statistic)
        elif mode == 'simultaneous':
            self._plot_spectra_simultaneous(table_name, log_prefix,
                                            skip_varabs, self._NH/(1e22),
                                            separate, rebin, rebin_params,
                                            rescale_F, rescale_chi, abund,
                                            skip_eRASS, varabs_starting_pars,
                                            plot_command, model_file, title,
                                            save_settings, log_suffix, colors,
                                            markers, fit_statistic)
        elif mode == 'merged':
            self._plot_spectra_merged(log_prefix, skip_eRASS, table_name,
                                      skip_varabs, self._NH/(1e22), rebin,
                                      rebin_params, rescale_F, rescale_chi,
                                      abund, varabs_starting_pars,
                                      plot_command, model_file, title,
                                      save_settings, log_suffix, colors,
                                      markers, fit_statistic)

    def _plot_spectra_simultaneous(self, table_name, log_prefix, skip_varabs,
                                   absorption, separate, rebin, rebin_params,
                                   rescale_F, rescale_chi, abund, skip_eRASS,
                                   varabs_starting_pars, plot_command,
                                   model_file, title, save_settings,
                                   log_suffix, colors, markers, fit_statistic):
        self._logger.info('Running _plot_spectra_simultaneous (log done by '
                          'xspec).')
        if self._create_epochs:
            period = 'epoch'
        else:
            period = f'e{self._ownership}'
        file_list_xspec = ''
        for epoch_counter in range(len(self._period_names)):
            if (period != 'epoch'
                    and (skip_eRASS == int(self._period_names[epoch_counter]
                                           [-1])).any()):
                continue
            start = ((self._obs_periods[epoch_counter][0] - self._mjdref) * 24.
                     * 3600.)
            stop = ((self._obs_periods[epoch_counter][1] - self._mjdref) * 24.
                    * 3600.)
            replacements = [['@esass_location', self._esass],
                            ['@infiles', self._filelist],
                            ['@outfile',
                            f'./{self._period_names[epoch_counter]}_'
                             'simultaneous.fits'],
                            ['@start', f'{start}'], ['@stop', f'{stop}']]
            sh_file = self._working_dir_full + '/working/trim_eventfile.sh'
            sh_file = self._replace_in_sh(sh_file, replacements)
            process = subprocess.Popen(
                [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()  # Wait for process to complete.

            self._extract_spectrum(f'{log_prefix}_simultaneous_'
                                   f'{self._period_names[epoch_counter]}'
                                   f'{log_suffix}', 'simultaneous',
                                   f'./{self._period_names[epoch_counter]}_'
                                   'simultaneous.fits',
                                   f'{self._period_names[epoch_counter]}')

            if epoch_counter == 0:
                file_list_xspec += (f'./{self._src_name}_{self._skytile}_'
                                    f'{self._period_names[epoch_counter]}_'
                                    'eROSITA_simultaneous_PATall_TMon020_'
                                    f'SourceSpec_00001_g{self._grouping}.fits')
            else:
                file_list_xspec += (f' {epoch_counter+1}:{epoch_counter+1} '
                                    f'./{self._src_name}_{self._skytile}_'
                                    f'{self._period_names[epoch_counter]}_'
                                    'eROSITA_simultaneous_PATall_TMon020_'
                                    f'SourceSpec_00001_g{self._grouping}.fits')

        self._standard_spec_an(separate, skip_eRASS, table_name, 'simultaneous',
                               file_list_xspec, skip_varabs, absorption, rebin,
                               rebin_params, rescale_F, rescale_chi, abund,
                               period, varabs_starting_pars, plot_command,
                               model_file, title, save_settings, log_suffix,
                               colors, markers, fit_statistic)

        if self._create_epochs:
            file_list_xspec = ''
            for epoch_counter in range(len(self._ero_starttimes)):
                if (skip_eRASS
                        == int(self._period_names[epoch_counter][-1])).any():
                    continue
                if epoch_counter == 0:
                    epoch = '0_1'
                    start = (58500 - self._mjdref) * 24. * 3600.
                    stop = ((self._ero_starttimes[epoch_counter + 1]
                             - self._mjdref) * 3600. * 24.)
                elif epoch_counter == len(self._ero_starttimes) - 1:
                    epoch = f'{epoch_counter + 1}'
                    start = ((self._ero_starttimes[epoch_counter]
                              - self._mjdref) * 3600. * 24.)
                    stop = ((self._ero_starttimes[epoch_counter] + 200
                             - self._mjdref) * 3600. * 24.)
                else:
                    epoch = f'{epoch_counter + 1}'
                    start = ((self._ero_starttimes[epoch_counter]
                              - self._mjdref) * 3600. * 24.)
                    stop = ((self._ero_starttimes[epoch_counter + 1]
                             - self._mjdref) * 3600. * 24.)

                replacements = [['@esass_location', self._esass],
                                ['@infiles', self._filelist],
                                ['@outfile',
                                f'./e{self._ownership}{epoch}_simultaneous'
                                 '.fits'],
                                ['@start', f'{start}'],
                                ['@stop', f'{stop}']]
                sh_file = self._working_dir_full + '/working/trim_eventfile.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                self._extract_spectrum(f'{log_prefix}_simultaneous_e'
                                       f'{self._ownership}{epoch}{log_suffix}',
                                       'simultaneous', f'./e{self._ownership}'
                                       f'{epoch}_simultaneous.fits',
                                       f'e{self._ownership}{epoch}')

                if epoch_counter == 0:
                    file_list_xspec += (f'./{self._src_name}_{self._skytile}_e'
                                        f'{self._ownership}{epoch}_eROSITA_'
                                        f'simultaneous_PATall_TMon020_Source'
                                        f'Spec_00001_g{self._grouping}.fits')
                else:
                    file_list_xspec += (f' {epoch_counter+1}:{epoch_counter+1}'
                                        f' ./{self._src_name}_{self._skytile}_'
                                        f'e{self._ownership}{epoch}_eROSITA_'
                                        f'simultaneous_PATall_TMon020_Source'
                                        f'Spec_00001_g{self._grouping}.fits')

            self._standard_spec_an(separate, skip_eRASS, table_name,
                                   'simultaneous', file_list_xspec,
                                   skip_varabs, absorption, rebin,
                                   rebin_params, rescale_F, rescale_chi, abund,
                                   f'e{self._ownership}', varabs_starting_pars,
                                   plot_command, model_file, title,
                                   save_settings, log_suffix, colors, markers,
                                   fit_statistic)

    def _plot_spectra_merged(self, log_prefix, skip_eRASS, table_name,
                             skip_varabs, absorption, rebin, rebin_params,
                             rescale_F, rescale_chi, abund,
                             varabs_starting_pars, plot_command, model_file,
                             title, save_settings, log_suffix, colors,
                             markers, fit_statistic):
        self._logger.info('Running _plot_spectra_merged (log done by xspec).')
        if self._create_epochs:
            period = 'epoch'
        else:
            period = f'e{self._ownership}'
        suffix = ''
        for entry in self._period_names:
            suffix += entry[len(period):]
        self._extract_spectrum(f'{log_prefix}_merged_{period}{suffix}'
                               f'{log_suffix}', 'merged', self._filelist,
                               f'{period}{suffix}')

        file_list_xspec = (f'./{self._src_name}_{self._skytile}_{period}'
                           f'{suffix}_eROSITA_merged_PATall_TMon020_SourceSpec'
                           f'_00001_g{self._grouping}.fits')

        self._standard_spec_an(False, skip_eRASS, table_name,
                               f'merged_{period}{suffix}', file_list_xspec,
                               skip_varabs, absorption, rebin, rebin_params,
                               rescale_F, rescale_chi, abund, period,
                               varabs_starting_pars, plot_command, model_file,
                               title, save_settings, log_suffix, colors,
                               markers, fit_statistic)

    def _plot_spectra_individual(self, table_name, log_prefix, skip_eRASS,
                                 skip_varabs, absorption, rebin, rebin_params,
                                 rescale_F, rescale_chi, abund,
                                 varabs_starting_pars, plot_command,
                                 model_file, title, save_settings, log_suffix,
                                 colors, markers, fit_statistic):
        self._logger.info('Running _plot_spectra_individual (log done by '
                          'xspec).')
        if self._create_epochs:
            period = 'epoch'
        else:
            period = f'e{self._ownership}'
        for epoch_counter in range(len(self._period_names)):
            if (period != 'epoch'
                    and (skip_eRASS == int(self._period_names[epoch_counter]
                                           [-1])).any()):
                continue
            start = ((self._obs_periods[epoch_counter][0] - self._mjdref) * 24.
                     * 3600.)
            stop = ((self._obs_periods[epoch_counter][1] - self._mjdref) * 24.
                    * 3600.)
            replacements = [['@esass_location', self._esass],
                            ['@infiles', self._filelist],
                            ['@outfile',
                            f'./{self._period_names[epoch_counter]}_'
                             'individual.fits'],
                            ['@start',
                            f'{start}'],
                            ['@stop', f'{stop}']]
            sh_file = self._working_dir_full + '/working/trim_eventfile.sh'
            sh_file = self._replace_in_sh(sh_file, replacements)
            process = subprocess.Popen(
                [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()  # Wait for process to complete.

            self._extract_spectrum(f'{log_prefix}_individual_'
                                   f'{self._period_names[epoch_counter]}'
                                   f'{log_suffix}', 'individual',
                                   f'./{self._period_names[epoch_counter]}_'
                                   'individual.fits',
                                   f'{self._period_names[epoch_counter]}')

            file_list_xspec = (f'./{self._src_name}_{self._skytile}_'
                               f'{self._period_names[epoch_counter]}_eROSITA_'
                               'individual_PATall_TMon020_SourceSpec_00001_'
                               f'g{self._grouping}.fits')

            temp_name = f'{self._period_names[epoch_counter][len(period):]}'
            self._standard_spec_an(False, skip_eRASS, table_name, 'individual '
                                   f'{temp_name}', file_list_xspec,
                                   skip_varabs, absorption, rebin,
                                   rebin_params, rescale_F, rescale_chi, abund,
                                   period, varabs_starting_pars, plot_command,
                                   model_file, title, save_settings,
                                   log_suffix, colors, markers, fit_statistic)

        if self._create_epochs:
            for epoch_counter in range(len(self._ero_starttimes)):
                if (skip_eRASS
                        == int(self._period_names[epoch_counter][-1])).any():
                    continue
                if epoch_counter == 0:
                    epoch = '0_1'
                    start = (58500 - self._mjdref) / 24. / 3600.
                    stop = (self._ero_starttimes[epoch_counter + 1]
                            - self._mjdref) * 3600. * 24.
                elif epoch_counter == len(self._ero_starttimes) - 1:
                    epoch = f'{epoch_counter + 1}'
                    start = (self._ero_starttimes[epoch_counter]
                             - self._mjdref) * 3600. * 24.
                    stop = (self._ero_starttimes[epoch_counter] + 200
                            - self._mjdref) * 3600. * 24.
                else:
                    epoch = f'{epoch_counter + 1}'
                    start = (self._ero_starttimes[epoch_counter]
                             - self._mjdref) * 3600. * 24.
                    stop = (self._ero_starttimes[epoch_counter + 1]
                            - self._mjdref) * 3600. * 24.

                replacements = [['@esass_location', self._esass],
                                ['@infiles', self._filelist],
                                ['@outfile', f'./e{self._ownership}{epoch}_'
                                'individual.fits'],
                                ['@start', f'{start}'], ['@stop', f'{stop}']]
                sh_file = self._working_dir_full + '/working/trim_eventfile.sh'
                sh_file = self._replace_in_sh(sh_file, replacements)
                process = subprocess.Popen(
                    [sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()  # Wait for process to complete.

                self._extract_spectrum(f'{log_prefix}_individual_e'
                                       f'{self._ownership}{epoch}{log_suffix}',
                                       'individual', f'./e{self._ownership}'
                                       f'{epoch}_individual.fits',
                                       f'e{self._ownership}{epoch}')

                file_list_xspec = (f'./{self._src_name}_{self._skytile}_e'
                                   f'{self._ownership}{epoch}_eROSITA_'
                                   'individual_PATall_TMon020_SourceSpec_00001'
                                   f'_g{self._grouping}.fits')

                temp_name = self._period_names[epoch_counter][len(period):]
                self._standard_spec_an(False, skip_eRASS, table_name,
                                       f'individual {temp_name}',
                                       file_list_xspec, skip_varabs,
                                       absorption, rebin, rebin_params,
                                       rescale_F, rescale_chi, abund,
                                       f'e{self._ownership}',
                                       varabs_starting_pars, plot_command,
                                       model_file, title, save_settings,
                                       log_suffix, colors, markers,
                                       fit_statistic)

    def plot_spectra_bayesian(self):
        fit_bxa(Xset, Fit, PlotManager, AllData, AllModels, Spectrum, Model,
                abund, distance, E_ranges, func, galnh, log, prompting, quantiles,
                src_files, statistic, suffix, resume, working_dir, Z)
        plot_bxa(Plot, rebinning, src_files, ax_spec, ax_res, colors,
                 src_markers, bkg_markers, bkg_linestyle, epoch_type)
        return

    def _standard_spec_an(self, separate, skip_eRASS, table_name, mode,
                          file_list, skip_varabs, absorption, rebin,
                          rebin_params, rescale_F, rescale_chi, abund, period,
                          varabs_starting_pars, plot_command, model_file,
                          title, save_settings, log_suffix, colors, markers,
                          fit_statistic):
        bands = {}
        if separate:
            list_visibles = range(1, len(file_list.split(sep=' ')) + 1)
        else:
            list_visibles = [-1]
        for t, energy_bin in enumerate(self._energy_bins):
            for visible in list_visibles:
                if mode.find(' ') == -1:
                    bands[f'table_{t}'] = open(f'{table_name}_{mode}_'
                                               f'{energy_bin[0]}keV_'
                                               f'{energy_bin[1]}keV.tex', 'w')
                else:
                    bands[f'table_{t}'] = open(f'{table_name}_'
                                               f'{mode.split()[0]}_'
                                               f'{energy_bin[0]}keV_'
                                               f'{energy_bin[1]}keV_{period}'
                                               f'{mode.split()[1]}.tex', 'w')

                bands[f'table_{t}'].write('\\begin{{tabular}}{{cccccc}}\n')
                bands[f'table_{t}'].write('\\hline\\hline\n')
                bands[f'table_{t}'].write('Data & Part & Power-law & N$_'
                                          '{{\\rm H, varab}}$ & \\mbox'
                                          '{{F$_{{\\rm x}}$}} & \\mbox'
                                          '{{L$_{{\\rm x}}$}} \\\\ \n')
                bands[f'table_{t}'].write('-- & -- & index & $\\times 10^'
                                          '{{20}}$ cm$^{{-2}}$ & $\\times$erg '
                                          'cm$^{{-2}}$s$^{{-1}}$ & $\\times'
                                          '$erg s$^{{-1}}$ \\\\ \n')
                bands[f'table_{t}'].write('\\hline\n')
                bands[f'table_{t}'].write('& & & & & \\\\ \n')

                AllData(file_list)
                AllData.ignore('bad')
                AllData.ignore('*:**-0.2 8.0-**')

                if model_file != '':
                    Xset.restore(self._working_dir + model_file)

                if mode.find(' ') != -1:
                    epoch = f'{period}{mode.split()[1]} {mode.split()[0]}'
                elif mode.find('_') != -1:
                    epoch = f'merged {mode.split(sep="_")[1]}'
                else:
                    epoch = f'{mode}'

                spec_model(Xset, AllModels, AllData, Model, Fit, Plot,
                           self._working_dir, bands[f'table_{t}'], self._Z,
                           self._distance, skip_varabs, epoch, absorption,
                           separate, visible, rebin,
                           [rebin_params[0], rebin_params[1], rescale_F[0],
                            rescale_F[1], rescale_chi[0], rescale_chi[1]],
                           abund, energy_bin[0], energy_bin[1],
                           varabs_starting_pars, plot_command, title,
                           colors, markers, fit_statistic)

                if skip_varabs:
                    parts = ['1', '3_1', '3_2']
                else:
                    parts = ['1', '2', '3_1', '3_2', '3_3', '3_4']

                for part in parts:
                    for extension in ['.ps', '.qdp', '.pco']:
                        if separate:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{visible}{extension}')
                        elif mode.find(' ') != -1:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode.split()[0]}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{mode.split()[1]}{extension}')
                        elif mode.find('_') != -1:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode.split()[0]}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{mode.split(sep="_")[1]}'
                                      f'{extension}')
                        else:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV'
                                      f'{extension}')
                    for extension in ['.log']:
                        if separate:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{visible}{log_suffix}')
                        elif mode.find(' ') != -1:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode.split()[0]}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{mode.split()[1]}{log_suffix}')
                        elif mode.find('_') != -1:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode.split()[0]}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                      f'{period}{mode.split(sep="_")[1]}'
                                      f'{log_suffix}')
                        else:
                            os.rename(f'./xspec_part{part}{extension}',
                                      f'{self._working_dir}/results/spectra/'
                                      f'xspec_part{part}_{mode}_'
                                      f'{energy_bin[0]}keV_{energy_bin[1]}keV'
                                      f'{log_suffix}')

                if mode.find(' ') == -1:
                    xcm_name = (f'{table_name}_{mode}_{energy_bin[0]}keV_'
                                f'{energy_bin[1]}keV.xcm')
                else:
                    xcm_name = (f'{table_name}_{mode.split()[0]}_'
                                f'{energy_bin[0]}keV_{energy_bin[1]}keV_'
                                f'{period}{mode.split()[1]}.xcm')
                if save_settings:
                    Xset.save(xcm_name)

        for part_table in bands:
            bands[part_table].write('\\hline\n')
            bands[part_table].write("\\end{{tabular}}")
            bands[part_table].close()

    def wipe_working_dir(self):
        '''Task to clean up working directory to free up space.
        '''
        self._logger.info('Cleaning up working directory.')
        for (path, dirs, filenames) in os.walk(f'{self._working_dir_full}'
                                               '/working'):
            for filename in filenames:
                os.remove(os.path.join(path, filename))
        self._LC_extracted = False
        self._debugging = False

    def run_standard(self):
        '''
        Returns
        -------
        Runs entire analysis chain with standard settings.

        '''
        # self.plot_lc_full()
        # self.plot_lc_parts()
        self.debug = False
        self.plot_lc_broken()
        self.plot_spectra()
        self.plot_lc_bayes_broken()
