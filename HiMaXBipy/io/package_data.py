#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Kaltenbrunner
"""
import pkg_resources
from math import log10, floor
import os
import shutil


def round_to_1(x):
    if x == 0:
        return x
    else:
        return round(x, -int(floor(log10(abs(x)))))


def get_path_of_data_file(data_file):
    file_path = pkg_resources.resource_filename("HiMaXBipy",
                                                'sh_files/%s' % data_file)

    return file_path


def get_path_of_data_dir():
    file_path = pkg_resources.resource_filename("HiMaXBipy", 'sh_files')

    return file_path

def get_stan_dir():
    file_path = pkg_resources.resource_filename("HiMaXBipy", 'stan_files')

    return file_path


def install_tex_sty():
    file_path = pkg_resources.resource_filename("HiMaXBipy",
                                                'tex_style/type1ec.sty')
    if not os.path.exists(os.path.expanduser('~/texmf')):
        os.mkdir(os.path.expanduser('~/texmf'))
        for subfolder in ['/bibtex', '/doc', '/tex']:
            os.mkdir(os.path.expanduser('~/texmf') + subfolder)
        for subfolder in ['/bib', '/bst']:
            os.mkdir(os.path.expanduser('~/texmf/bibtex') + subfolder)
        for subfolder in ['/context', '/generic', '/latex',
                          '/plain', '/xelatex', '/xetex']:
            os.mkdir(os.path.expanduser('~/texmf/tex') + subfolder)
    else:
        if not os.path.exists(os.path.expanduser('~/texmf/bibtex')):
            os.mkdir(os.path.expanduser('~/texmf/bibtex'))
            for subfolder in ['/bib', '/bst']:
                os.mkdir(os.path.expanduser('~/texmf/bibtex') + subfolder)
        else:
            for subfolder in ['/bib', '/bst']:
                if (not os.path.exists(os.path.expanduser('~/texmf/bibtex')
                                       + subfolder)):
                    os.mkdir(os.path.expanduser('~/texmf/bibtex') + subfolder)
        if not os.path.exists(os.path.expanduser('~/texmf/doc')):
            os.mkdir(os.path.expanduser('~/texmf/doc'))

        if not os.path.exists(os.path.expanduser('~/texmf/tex')):
            os.mkdir(os.path.expanduser('~/texmf/tex'))
            for subfolder in ['/context', '/generic', '/latex',
                              '/plain', '/xelatex', '/xetex']:
                os.mkdir(os.path.expanduser('~/texmf/tex') + subfolder)
        else:
            for subfolder in ['/context', '/generic', '/latex',
                              '/plain', '/xelatex', '/xetex']:
                if (not os.path.exists(os.path.expanduser('~/texmf/tex')
                                       + subfolder)):
                    os.mkdir(os.path.expanduser('~/texmf/tex') + subfolder)
    if not os.path.exists(os.path.expanduser('~/texmf/tex/latex/type1ec.sty')):
        shutil.copy(file_path, os.path.expanduser('~/texmf/tex/latex'))
        print('Folder structure successfully built, file copied.')
    else:
        print('File already exists.')
