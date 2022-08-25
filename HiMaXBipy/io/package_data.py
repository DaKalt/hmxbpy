import pkg_resources
from math import log10, floor

def round_to_1(x):
    if x == 0:
        return x
    else:
        return round(x, -int(floor(log10(abs(x)))))

def get_path_of_data_file(data_file):
     file_path = pkg_resources.resource_filename("HiMaXBipy", 'sh_files/%s' % data_file)

     return file_path

def get_path_of_data_dir():
     file_path = pkg_resources.resource_filename("HiMaXBipy", 'sh_files')

     return file_path
