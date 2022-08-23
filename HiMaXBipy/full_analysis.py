from HiMaXBipy.io.package_data import get_path_of_data_dir
import os
import shutil
import sys, fileinput
import subprocess

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
    
    def set_region(self, region):
        '''
        Parameters
        ----------
        region : string
            Name of the region in which the source lies.

        Returns
        -------
        Sets name of the region (e.g. 080156) in which the source lies.

        '''
        if type(region) == str:
            self.__region = region
        else:
            raise Exception('Not a valid region name.')
    
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

    def extract_lc(self, logname):
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
        replacements = [['source_name', self.__src_name],
                        ['main_name', self.__working_dir],
                        ['result_dir', self.__working_dir + '/working'],
                        ['region_code', self.__region],
                        ['sources_list', self.__filelist],
                        ['right_ascension', self.__RA],
                        ['declination', self.__Dec]]
        sh_file = self.__working_dir + '/working/extract_lc.sh'
        self.__replace_in_ssh(sh_file, replacements)
        process = subprocess.Popen([sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait() # Wait for process to complete.

        # iterate on the stdout line by line
        with open(self.__working_dir + '/logfiles/' + logname, 'w') as logfile:
            for line in process.stdout.readlines():
                logfile.writelines(line)
        