from HiMaXBipy.io.package_data import get_path_of_data_dir
import os
import shutil

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
            self.src_name = src_name
        else:
            raise Exception('src_name needs to be of type string.')
            
        if os.path.exists(working_dir) and type(working_dir) == str:
            self.working_dir = os.path.abspath(working_dir)
        else:
            raise Exception('Not a valid path for a working directory.')
            
        if os.path.exists(data_dir) and type(data_dir) == str:
            self.data_dir = os.path.abspath(data_dir)
        else:
            raise Exception('Not a valid path for a data directory.')
        self.data_files = ''
        
        self.__sh_dir__ = get_path_of_data_dir()
        
        for (path, directories, filenames) in os.walk(self.__sh_dir__):
            for filename in filenames:
                shutil.copy(self.__sh_dir__ + filename, self.working_dir)
        