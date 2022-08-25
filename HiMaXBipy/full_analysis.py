from HiMaXBipy.io.package_data import get_path_of_data_dir
from HiMaXBipy.lc_plotting.lc_plotting import plot_lc_UL, plot_lc_mincounts, get_boundaries, format_axis
import os
import shutil
import sys, fileinput
import subprocess
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import numpy as np
import getpass

class HiMaXBi:
    def __init__(self, src_name, working_dir, data_dir):
        '''
        Parameters
        ----------
        src_name : str
            Sets the name of the source used for naming files.
        working_dir : TYPE
            Sets the working directory to which resulting files will be saved.
        data_dir : TYPE
            Sets directory from where to use data files.

        Returns
        -------
        None.

        '''
        if type(src_name) == str:
            self.__src_name = src_name
        else:
            raise Exception('src_name needs to be of type string.')
            
        if os.path.exists(working_dir) and type(working_dir) == str:
            self.__working_dir = os.path.abspath(working_dir)
        else:
            raise Exception('Not a valid path for a working directory.')
            
        if os.path.exists(data_dir) and type(data_dir) == str:
            self.__data_dir = os.path.abspath(data_dir)
        else:
            raise Exception('Not a valid path for a data directory.')
        self.data_files = ''

        if not os.path.exists(self.__working_dir + '/working'):
            os.mkdir(self.__working_dir + '/working')
        if not os.path.exists(self.__working_dir + '/results'):
            os.mkdir(self.__working_dir + '/results')
        if not os.path.exists(self.__working_dir + '/logiles'):
            os.mkdir(self.__working_dir + '/logfiles')
        
        self.__sh_dir__ = get_path_of_data_dir()
        
        for (path, directories, filenames) in os.walk(self.__sh_dir__):
            for filename in filenames:
                shutil.copy(self.__sh_dir__ + filename, self.__working_dir + '/working')
                
        self.__region = ''
        self.__filenames = ''
        self.__esass = '/home/erosita/sw/eSASSusers_201009/bin/esass-init.sh'
        self.__LC_prebinning = '1.0'
        self.__LC_extracted = False
        self.__mjdref = 51543.875
        self.__ero_starttimes = [58828, 59011, 59198, 59381, 59567]
    
    def __replace_in_ssh(self, path, replacements):
        '''
        Parameters
        ----------
        replacements : array-like shape (n, 2)
            array of n keywords to replace in ssh file; 0th entry in each pair states the original keyword, 1st entry states the new keyword
        path : string
            path of ssh file in which to replace entries

        Returns
        -------
        None.

        '''
        for pair in replacements:
            for line in fileinput.input(path, inplace=True):
                line = line.replace(pair[0], pair[1])
                sys.stdout.write(line)
    
    def set_mjd_ref(self, mjdref):
        '''
        Parameters
        ----------
        mjdref : float or string
            Reference value to transform eROSITA times to MJD.
        Returns
        -------
        Sets reference value to transform eROSITA times to MJD.

        '''
        if type(mjdref) != float and type(mjdref) != str:
            raise Exception('mjdref must be a float or string.')
        try:
            self.__mjdref = float(mjdref)
        except TypeError:
            raise Exception('mjdref must be convertible to float.')
    
    def set_esass(self, esass_location):
        '''
        Parameters
        ----------
        esass_location : str
            full path of esass initialisation script.

        Returns
        -------
        Sets path of esass in use.

        '''
        if type(esass_location) != str:
            raise Exception('esass_location must be string.')
        if not os.path.exists(esass_location):
            raise Exception(f'File {esass_location} does not exist.')
        self.__esass = esass_location
        self.__LC_extracted = False
    
    def set_skytile(self, skytile):
        '''
        Parameters
        ----------
        skytile : string
            Name of the skytile in which the source lies.

        Returns
        -------
        Sets name of the skytile (e.g. 080156) in which the source lies.

        '''
        if type(skytile) == str:
            self.__skytile = skytile
            self.__LC_extracted = False
        else:
            raise Exception('Not a valid skytile name (must be string).')
    
    def set_radec(self, RA, Dec):
        '''
        Parameters
        ----------
        RA : float
            Right ascension of the source in J2000.
        Dec : TYPE
            Declination of the source in J2000.

        Returns
        -------
        Sets RA and Dec for the source.

        '''
        if type(RA) != float or type(Dec) != float:
            raise Exception('Ra and Dec must be given as floats.')
        self.__RA = RA
        self.__Dec = Dec
        self.__LC_extracted = False
    
    def set_filelist(self, filelist):
        '''
        Parameters
        ----------
        filelist : string
            List of names of eventfiles to use separated by spaces and without foldernames.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if type(filelist) != str:
            raise Exception('filelist must be string.')
        temp1 = filelist
        temp3 = ''
        while temp1.find(' ') != -1:
            temp2 = temp1[:temp1.find(' ')]
            temp1 = temp1[temp1.find(' '):].strip()
            if not os.path.exists(self.__data_dir + '/' + temp2):
                raise Exception(f'File {temp2} does not exist.')
            temp3 += self.__data_dir + '/' + temp2
        if not os.path.exists(self.__data_dir + '/' + temp1):
            raise Exception(f'File {temp1} does not exist.')
        temp3 += self.__data_dir + '/' + temp1
        self.__filelist = temp3
        self.__LC_extracted = False
    
    def set_LC_binning(self, binning):
        '''
        Parameters
        ----------
        binning : string or float
            Set the initial binning of the lightcurve in seconds.

        Returns
        -------
        None.

        '''
        if not (type(binning) == str or type(binning) == float):
            raise Exception('binning must be a string or float.')
        self.__binning = str(binning)
        self.__LC_extracted = False
    
    def extract_lc(self, logname = 'lc_extract_autosave.log'):
        '''
        Parameters
        ----------
        logname : string
            name of logfile to safe output of sh script to.
            
        Returns
        -------
        Creates sh files to extract light curve fits files and runs them.

        '''
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if self.__region == '' or self.__filelist == '':
            raise Exception('Set the region name and list of eventfiles first with the functions set_filelist and set_region.')
        replacements = [['@source_name', self.__src_name],
                        ['@main_name', self.__working_dir],
                        ['@result_dir', self.__working_dir + '/working'],
                        ['@region_code', self.__skytile],
                        ['@sources_list', self.__filelist],
                        ['@right_ascension', self.__RA],
                        ['@declination', self.__Dec],
                        ['@esass_location', self.__esass],
                        ['@binning', self.__binning]]
        sh_file = self.__working_dir + '/working/extract_lc.sh'
        self.__replace_in_ssh(sh_file, replacements)
        process = subprocess.Popen([sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait() # Wait for process to complete.

        # iterate on the stdout line by line
        if not logname == '':
            with open(self.__working_dir + '/logfiles/' + logname, 'w') as logfile:
                for line in process.stdout.readlines():
                    logfile.writelines(line)
        self.__LC_extracted = True
        self.__find_gap_centres(60 * 60 * 24 * 30) #one month gap minimum to sort out possible short gaps due to problems during observation
        self.__eRASS_vs_epoch()
    
    def __find_gap_centres(self, gapsize):
        with fits.open(f'{self.__working_dir}/working/{self.__src_name}_{self.__skytile}_eROSITA_PATall_1.0s020_LightCurve_00001.fits') as hdulist:
            time = hdulist[1].data.field('TIME').tolist().sort()
        self.__gap_centres = []
        for i in range(len(time) - 1):
            if time[i + 1] - time[i] > gapsize:
                self.__gap_centres.append((time[i + 1] - time[i]) / 2.)
    
    def __eRASS_vs_epoch(self):
        eRASSi = self.__filelist.split(sep = ' ')
        self.__create_epochs = False
        times_list = []
        for entry in eRASSi:
            with fits.open(entry) as hdulist:
                times_list.append(hdulist[1].data.field('TIME').tolist())
        for i in range(len(eRASSi)):
            for entry in self.__gap_centres:
                if entry > min(times_list[i]) and entry < max(times_list[i]):
                    self.__create_epochs = True
    
    def plot_lc_full(self, fracexp = '0.15', mincounts = '10', mode = 'ul', show_eRASS = True, logname = 'lc_full_autosave.log', time_axis = 'mjd', print_name = False, print_datetime = False, label_style = 'serif', label_size = 12, figsize = [8, 5.5], colors = [], fileid = '', toplab = '', separate_TM = False, vlines = [], ticknumber_y = 5.0, ticknumber_x = 8.0):
        '''
        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted as an upper limit as well as minimum number of counts per bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or 'mincounts_ul'. The default is 'ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines. The default is True.
        logname : str, optional
            Name of the logfile. The default is 'lc_full_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 5.5].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib possible. For mode 'ul' and 'mincounts' the first entry is used, for mode 'mincounts_ul' the first entry sets color for 'ul' part, and the second for 'mincounts' part. The default is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with given color. The zorder entries need to be distinct negative integers < -2. The default is [].
        ticknumber_y : float or int, optional
            Sets the approximate number of tickmarks along the y axis. The default is 5.0.
        ticknumber_x : float or int, optional
            Sets the approximate number of tickmarks along the x axis. The default is 8.0.

        Returns
        -------
        Creates full LC.

        '''
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if type(mincounts) != str and type(mincounts) != float and type(mincounts) != int:
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
                raise Exception('mode must be \'ul\', \'mincounts\' or \'mincounts_ul\'')
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
        if (type(figsize) != list and type(figsize) != np.ndarray) or np.shape(figsize) != (2,):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != float and type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be a float or int.')
        if type(ticknumber_y) != float and type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be a float or int.')
        if (type(colors) != list and type(colors) != np.ndarray) or (np.shape(colors) != (2,) and np.shape(colors) != (1,)):
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
                    raise Exception('The first entry in each line of vlines needs to be the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception('The second entry in each line of vlines needs to be a matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception('The third entry in each line of vlines needs to be a negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception('The third entry in each line of vlines needs to be a negative integer < -2.')
        
        if not self.__LC_extracted:
            self.extract_lc()
        logfile = open(self.__working_dir + '/logfiles/' + logname, 'w')
        localtime = time.asctime(time.localtime(time.time()))
        
        logfile.writelines(localtime + '\n')
        user = getpass.getuser()
        plt.rc('text', usetex=True)
        plt.rc('font', family=label_style, size = label_size)
        
        if separate_TM:
            TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            TM_list = [0]
        
        for TM in TM_list:
            if fileid == '':
                pfile = f'{self.__working_dir}/working/{self.__src_name}_{self.__skytile}_LC_TM{TM}20_fracexp{fracexp}_fullLC'
                outfile = f'{self.__working_dir}/results/{self.__src_name}_{self.__skytile}_LC_TM{TM}20_fracexp{fracexp}_fullLC'
            else:
                pfile = f'{self.__working_dir}/working/{fileid}'
                outfile = f'{self.__working_dir}/results/{fileid}'
            replacements = [['@infile', f'{self.__working_dir}/working/{self.__src_name}_{self.__skytile}_eROSITA_PATall_1.0s{TM}20_LightCurve_00001.fits'],
                            ['@pfile', f'{pfile}.fits'],
                            ['@selection', f'FRACEXP>{fracexp}']]
            sh_file = self.__working_dir + '/working/fselect_lc.sh'
            self.__replace_in_ssh(sh_file, replacements)
            
            fig1 = plt.figure(figsize=(figsize[0],figsize[1]))
            ax = fig1.add_subplot(111)
            
            logfile.writelines(f'Now working on {pfile}.fits')
            hdulist = fits.open(f'{pfile}.fits')
            
            if time_axis == 'mjd':
                xflag = 2
            elif time_axis == 's':
                xflag = 1
            
            if colors == []:
                colors = ['lightblue', 'black']
            
            if mode == 'ul':
                pxmin, pxmax, pymin, pymax = plot_lc_UL(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[0])
            elif mode == 'mincounts':
                pxmin, pxmax, pymin, pymax = plot_lc_mincounts(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[0])
            elif mode == 'mincounts_ul':
                pxmin1, pxmax1, pymin1, pymax1 = plot_lc_UL(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, colors = colors[0])
                pxmin2, pxmax2, pymin2, pymax2 = plot_lc_mincounts(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[1])
                pxmin, pxmax, pymin, pymax = get_boundaries([[pxmin1, pxmax1, pymin1, pymax1], [pxmin2, pxmax2, pymin2, pymax2]])
            
            format_axis(ax, pxmin, pxmax, pymin, pymax, ticknumber_x, ticknumber_y)
            
            # plot time in s from beginning (xflag=1) or in MJD
            if time_axis == 's':
               ax.set_xlabel(r'Time (s)')#, fontsize=12)
            elif time_axis == 'mjd':
               ax.set_xlabel(r'MJD (days)')#, fontsize=12)
        
            ax.set_ylabel(r'Count rate (cts/s)')#, fontsize=12)
            
            if print_name:
            # user name and time
                ax.text(1.015, 0.0, user +' - '+ localtime, rotation=90, fontsize=8, 
                         verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)
            # eROSITA label
                ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10, 
                         verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)
                ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10, 
                         verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes)
            # label plot:
                ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                         verticalalignment='bottom', horizontalalignment='center', transform=ax.transAxes)
            
            for i in range(len(vlines)):
                ax.vlines(vlines[i][0], -5, 5, colors = vlines[i][1], linestyle = 'dotted', zorder=vlines[i][2])
            if show_eRASS:
                if time_axis == 'mjd':
                    ax.vlines(self.__ero_starttimes, -5, 5, colors = 'grey', linestyle = 'dotted', zorder=-2)
                elif time_axis == 's':
                    ax.vlines((np.array(self.__ero_starttimes) - self.__mjdref) * 3600 * 24, -5, 5, colors = 'grey', linestyle = 'dotted', zorder=-4)
   
            fig1.tight_layout()
   
            pltfile = outfile + ".pdf"
            plt.savefig(pltfile)
            logfile.writelines(f'{pltfile} created')
            pltfile = outfile + ".eps"
            plt.savefig(pltfile)
            logfile.writelines(f'{pltfile} created')
            pltfile = outfile + ".png"
            plt.savefig(pltfile)
            logfile.writelines(f'{pltfile} created')
            
        logfile.close()

    def plot_lc_parts(self, fracexp = '0.15', mincounts = '10', mode = 'mincounts_ul', show_eRASS = True, logname = 'lc_parts_autosave.log', time_axis = 'mjd', print_name = False, print_datetime = False, label_style = 'serif', label_size = 12, figsize = [8, 5.5], colors = [], fileid = '', toplab = '', separate_TM = False, vlines = [], ticknumber_y = 5.0, ticknumber_x = 8.0, eRASSi = []):
        '''
        Parameters
        ----------
        fracexp : str or float, optional
            Fractional exposure lower limit for times taken into account for LC (noise reduction). The default is '0.15'.
        mincounts : str, float or int, optional
            Minimum number of counts for counts per bin to not be noted as an upper limit as well as minimum number of counts per bin for mode mincounts/mincounts_ul. The default is '10'.
        mode : str, optional
            Type of LC to be produced. Either 'ul', 'mincounts' or 'mincounts_ul'. The default is 'mincounts_ul'.
        show_eRASS : bool, optional
            True to show start/end dates of eRASSi as vertical lines. The default is True.
        logname : str, optional
            Name of the logfile. The default is 'lc_parts_autosave.log'.
        time_axis : str, optional
            Defines the unit of time axis. Either 'mjd' or 's'. The default is 'mjd'.
        print_name : bool, optional
            Print name of person who runs the skript. The default is False.
        print_datetime : bool, optional
            Print date-time when skript was run. The default is False.
        label_style : str, optional
            Sets fontstyle of plots. Any possible style available for matplotlib.pyplot.rc. The default is 'serif'.
        label_size : float or int, optional
            Sets fontsize. The default is 12.
        figsize : array-like (2,), optional
            Sets width and height of figure. The default is [8, 5.5].
        colors : array of str (1,) or (2,), optional
            Sets colors of plots. Any color available to matplotlib possible. For mode 'ul' and 'mincounts' the first entry is used, for mode 'mincounts_ul' the first entry sets color for 'ul' part, and the second for 'mincounts' part. The default is [].
        fileid : str, optional
            Name of outputfile without filespecific ending. The default is ''.
        toplab : str, optional
            Sets label of the plot. The default is ''.
        separate_TM : bool, optional
            Create LC for each TM. The default is False.
        vlines : array of mjd-color-zorder combinations (n, 3), optional
            Adds additional vertical lines in the plot at given MJD with given color. The zorder entries need to be distinct negative integers < -2. The default is [].
        ticknumber_y : float or int, optional
            Sets the approximate number of tickmarks along the y axis. The default is 5.0.
        ticknumber_x : float or int, optional
            Sets the approximate number of tickmarks along the x axis. The default is 8.0.
        eRASSi : array-like of ints
            List of from which eRASS eventfiles were used. The default is [].

        Returns
        -------
        Creates LC of eRASSi/epochs.

        '''
        if type(logname) != str:
            raise Exception('logname must be a string.')
        if type(mincounts) != str and type(mincounts) != float and type(mincounts) != int:
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
                raise Exception('mode must be \'ul\', \'mincounts\' or \'mincounts_ul\'')
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
        if (type(figsize) != list and type(figsize) != np.ndarray) or np.shape(figsize) != (2,):
            raise Exception('figsize must be (2,) array-like.')
        if type(ticknumber_x) != float and type(ticknumber_x) != int:
            raise Exception('ticknumber_x must be a float or int.')
        if type(ticknumber_y) != float and type(ticknumber_y) != int:
            raise Exception('ticknumber_y must be a float or int.')
        if (type(colors) != list and type(colors) != np.ndarray) or (np.shape(colors) != (2,) and np.shape(colors) != (1,)):
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
                    raise Exception('The first entry in each line of vlines needs to be the MJD given as float or int.')
                if type(line[1]) != str:
                    raise Exception('The second entry in each line of vlines needs to be a matplotlib color given as a string.')
                if type(line[2]) != int:
                    raise Exception('The third entry in each line of vlines needs to be a negative integer < -2.')
                elif line[2] >= -1:
                    raise Exception('The third entry in each line of vlines needs to be a negative integer < -2.')
        if type(eRASSi) != list and type(eRASSi) != np.ndarray:
            raise Exception('eRASSi must be array-like.')
        else:
            for i in range(len(eRASSi)):
                if type(eRASSi[i]) != int:
                    raise Exception('Entries in eRASSi must be integer.')
        if np.array(eRASSi) != np.sort(eRASSi):
            raise Exception('eRASSi must be sorted ascending.')
        if eRASSi != []:
            print('Use of eRASSi is currently not supported.')
        
        if not self.__LC_extracted:
            self.extract_lc()
        
        logfile = open(self.__working_dir + '/logfiles/' + logname, 'w')
        
        for i in range(1, len(self.__gap_centres) + 2):
            if self.__create_epochs:
                naming = f'epoch{i}'
            else:
                naming = f'em0{i}'
            
            localtime = time.asctime(time.localtime(time.time()))
            
            logfile.writelines(localtime + '\n')
            user = getpass.getuser()
            plt.rc('text', usetex=True)
            plt.rc('font', family=label_style, size = label_size)
            
            if separate_TM:
                TM_list = [0, 1, 2, 3, 4, 5, 6, 7]
            else:
                TM_list = [0]
            
            for TM in TM_list:
                if fileid == '':
                    pfile = f'{self.__working_dir}/working/{self.__src_name}_{self.__skytile}_LC_TM{TM}20_fracexp{fracexp}_LC_{naming}'
                    outfile = f'{self.__working_dir}/results/{self.__src_name}_{self.__skytile}_LC_TM{TM}20_fracexp{fracexp}_LC_{naming}'
                else:
                    pfile = f'{self.__working_dir}/working/{fileid}'
                    outfile = f'{self.__working_dir}/results/{fileid}'
                if i == 1:
                    selection = f'FRACEXP>{fracexp} && TIME < {self.__gap_centres[0]}'
                elif i == len(self.__gap_centres) + 1:
                    selection = f'FRACEXP>{fracexp} && TIME > {self.__gap_centres[-1]}'
                else:
                    selection = f'FRACEXP>{fracexp} && TIME < {self.__gap_centres[i-1]} && TIME > {self.__gap_centres[i-2]}'
                replacements = [['@infile', f'{self.__working_dir}/working/{self.__src_name}_{self.__skytile}_eROSITA_PATall_1.0s{TM}20_LightCurve_00001.fits'],
                                ['@pfile', f'{pfile}.fits'],
                                ['@selection', selection]]
                sh_file = self.__working_dir + '/working/fselect_lc.sh'
                self.__replace_in_ssh(sh_file, replacements)
                
                fig1 = plt.figure(figsize=(figsize[0],figsize[1]))
                ax = fig1.add_subplot(111)
                
                logfile.writelines(f'Now working on {pfile}.fits')
                hdulist = fits.open(f'{pfile}.fits')
                
                if time_axis == 'mjd':
                    xflag = 2
                elif time_axis == 's':
                    xflag = 1
                
                if colors == []:
                    colors = ['lightblue', 'black']
                
                if mode == 'ul':
                    pxmin, pxmax, pymin, pymax = plot_lc_UL(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[0])
                elif mode == 'mincounts':
                    pxmin, pxmax, pymin, pymax = plot_lc_mincounts(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[0])
                elif mode == 'mincounts_ul':
                    pxmin1, pxmax1, pymin1, pymax1 = plot_lc_UL(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, colors = colors[0])
                    pxmin2, pxmax2, pymin2, pymax2 = plot_lc_mincounts(hdulist = hdulist, ax = ax, logfile = logfile, mjdref = self.__mjdref, xflag = xflag, mincounts = mincounts, color = colors[1])
                    pxmin, pxmax, pymin, pymax = get_boundaries([[pxmin1, pxmax1, pymin1, pymax1], [pxmin2, pxmax2, pymin2, pymax2]])
                
                format_axis(ax, pxmin, pxmax, pymin, pymax, ticknumber_x, ticknumber_y)
                
                # plot time in s from beginning (xflag=1) or in MJD
                if time_axis == 's':
                   ax.set_xlabel(r'Time (s)')#, fontsize=12)
                elif time_axis == 'mjd':
                   ax.set_xlabel(r'MJD (days)')#, fontsize=12)
            
                ax.set_ylabel(r'Count rate (cts/s)')#, fontsize=12)
                
                if print_name:
                # user name and time
                    ax.text(1.015, 0.0, user +' - '+ localtime, rotation=90, fontsize=8, 
                             verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)
                # eROSITA label
                    ax.text(0.0, 1.015, 'eROSITA', rotation=0, fontsize=10, 
                             verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)
                    ax.text(1.0, 1.015, 'MPE', rotation=0, fontsize=10, 
                             verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes)
                # label plot:
                    ax.text(0.5, 1.015, toplab, rotation=0, fontsize=10,
                             verticalalignment='bottom', horizontalalignment='center', transform=ax.transAxes)
                
                for i in range(len(vlines)):
                    ax.vlines(vlines[i][0], -5, 5, colors = vlines[i][1], linestyle = 'dotted', zorder=vlines[i][2])
                if show_eRASS:
                    if time_axis == 'mjd':
                        ax.vlines(self.__ero_starttimes, -5, 5, colors = 'grey', linestyle = 'dotted', zorder=-2)
                    elif time_axis == 's':
                        ax.vlines((np.array(self.__ero_starttimes) - self.__mjdref) * 3600 * 24, -5, 5, colors = 'grey', linestyle = 'dotted', zorder=-4)
       
                fig1.tight_layout()
       
                pltfile = outfile + ".pdf"
                plt.savefig(pltfile)
                logfile.writelines(f'{pltfile} created')
                pltfile = outfile + ".eps"
                plt.savefig(pltfile)
                logfile.writelines(f'{pltfile} created')
                pltfile = outfile + ".png"
                plt.savefig(pltfile)
                logfile.writelines(f'{pltfile} created')
            
        logfile.close()
    
    def run_standard(self):
        '''
        Returns
        -------
        Runs entire analysis chain with standard settings.

        '''
        self.plot_lc_full()
        self.plot_lc_parts()
        #XXX