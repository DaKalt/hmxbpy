# HiMaXbipy
Tool to analyse eROSITA data of HMXB (lightcurve, spectrum, RGB image)

Requires:
-python3
-HEASOFT (possibility to run xspec in python)

Require latex style 'type1ec.sty', which is part of the package cm-super to run the necessary version of matplotlib.
If the style is missing (matplotlib will cause an error, stating that the style file is needed), it can be installed
in the home directory running the following commands in python:

from HiMaXBipy.io.package_data import install_tex_sty
install_tex_sty()

(Latex needs to look for styles in ~/texmf, but it usually does.)
