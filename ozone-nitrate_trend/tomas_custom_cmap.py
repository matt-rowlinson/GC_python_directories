#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as colors
import os

def get_colormap(fileName):
    """ Retrieve a colormap from a file and return as a cmap object"""
    # Set location of file
    folder = 'custom_cmaps/'
    #    cmap =  gmtColormap(fileName, GMTPath=folder)
    cmap = gmtColormap(folder+fileName)
    return cmap


def gmtColormap_openfile(cptf, name=None):
    """
    Read a GMT color map from an OPEN cpt file
    Parameters
    -------
    cptf : open file or url handle
    path to .cpt file
    name : str, optional
    name for color map
    if not provided, the file name will be used
    Returns
    -------
    (cmap)
    Notes
    -------
    - Credit: j08lue (https://github.com/j08lue/pycpt)
    """
    import matplotlib.colors as mcolors
    # generate cmap name
    if name is None:
        name = '_'.join(os.path.basename(cptf.name).split('.')[:-1])
    # process file
    x = []
    r = []
    g = []
    b = []
    lastls = None
    for l in cptf.readlines():
        ls = l.split()
        # skip empty lines
        if not ls:
            continue
        # parse header info
        if ls[0] in ["#", b"#"]:
            if ls[-1] in ["HSV", b"HSV"]:
                colorModel = "HSV"
            else:
                colorModel = "RGB"
            continue
        # skip BFN info
        if ls[0] in ["B", b"B", "F", b"F", "N", b"N"]:
            continue
        # parse color vectors
        x.append(float(ls[0]))
        r.append(float(ls[1]))
        g.append(float(ls[2]))
        b.append(float(ls[3]))
        # save last row
        lastls = ls
    x.append(float(lastls[4]))
    r.append(float(lastls[5]))
    g.append(float(lastls[6]))
    b.append(float(lastls[7]))
    x = np.array(x)
    r = np.array(r)
    g = np.array(g)
    b = np.array(b)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            # convert HSV to RGB
            rr, gg, bb = colorsys.hsv_to_rgb(r[i]/360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    elif colorModel == "RGB":
        r /= 255.
        g /= 255.
        b /= 255.
    red = []
    blue = []
    green = []
    xNorm = (x - x[0])/(x[-1] - x[0])
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])
    # return colormap
    cdict = dict(red=red, green=green, blue=blue)
    return mcolors.LinearSegmentedColormap(name=name, segmentdata=cdict)

def gmtColormap(cptfile, name=None):
    """
    Read a GMT color map from a cpt file
    Parameters
    -------
    cptfile : str or open file-like object
    path to .cpt file
    name : str, optional
    name for color map
    if not provided, the file name will be used
    Returns
    -------
    (cmap)
    Notes
    -------
    - Credit: j08lue (https://github.com/j08lue/pycpt)
    """
    with open(cptfile, 'r') as cptf:
        return gmtColormap_openfile(cptf, name=name)

def main():
    fileName = 'TMS_custom_colormap_CMRmap.cpt'
    cmap=get_colormap(fileName)
    return cmap

if __name__ == "__main__":
    main()
