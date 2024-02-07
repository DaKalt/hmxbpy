import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pyregion

def plot_rgb(figsize, label_style, label_size, src_color, bkg_color,
             src_region, bkg_region, fig_borders, red_lims, green_lims,
             blue_lims, non_linearity, size_factor, r_file, g_file, b_file):
    plt.rc('text', usetex=True)
    plt.rc('font', family=label_style, size=label_size)

    # image
    r_min = red_lims[0]
    r_max = red_lims[1]
    g_min = green_lims[0]
    g_max = green_lims[1]
    b_min = blue_lims[0]
    b_max = blue_lims[1]

    r_img = fits.open(r_file)[0].data
    g_img = fits.open(g_file)[0].data
    b_img = fits.open(b_file)[0].data
    header_r = fits.open(r_file)[0].header

    img = np.zeros((r_img.shape[0], r_img.shape[1], 3), dtype=float)
    wcs = WCS(header_r)

    # regions
    src_region = pyregion.open(src_region).as_imagecoord(header_r)
    patch_list_src, artist_list_src = src_region.get_mpl_patches_texts()

    bkg_region = pyregion.open(bkg_region).as_imagecoord(header_r)
    patch_list_bkg, artist_list_bkg = bkg_region.get_mpl_patches_texts()

    maxs_x = []
    maxs_y = []
    mins_x = []
    mins_y = []

    for p in patch_list_src + patch_list_bkg:
        center = p.get_center()
        width = p.get_width()
        height = p.get_height()
        maxs_x.append(center[0] + width / 2.)
        maxs_y.append(center[1] + height / 2.)
        mins_x.append(center[0] - width / 2.)
        mins_y.append(center[1] - height / 2.)

    max_x = max(maxs_x)
    max_y = max(maxs_y)
    min_x = min(mins_x)
    min_y = min(mins_y)

    center_x = (max_x + min_x) / 2.
    center_y = (max_y + min_y) / 2.
    width = max_x - min_x
    height = max_y - min_y
    size = size_factor * max(width, height)

    r_img[r_img<r_min] = r_min
    r_img[r_img>r_max] = r_max
    g_img[g_img<g_min] = g_min
    g_img[g_img>g_max] = g_max
    b_img[b_img<b_min] = b_min
    b_img[b_img>b_max] = b_max

    r_img = (r_img - r_min) / (r_max - r_min)
    g_img = (g_img - g_min) / (g_max - g_min)
    b_img = (b_img - b_min) / (b_max - b_min)

    img[:,:,0] = np.arcsinh(non_linearity * r_img) / np.arcsinh(non_linearity)
    img[:,:,1] = np.arcsinh(non_linearity * g_img) / np.arcsinh(non_linearity)
    img[:,:,2] = np.arcsinh(non_linearity * b_img) / np.arcsinh(non_linearity)

    # plotting
    fig = plt.figure(figsize = (figsize,0.95*figsize))
    ax = fig.add_subplot(111, projection=wcs)

    ax.imshow(img, aspect = 'equal')
    for p in patch_list_src:
        p.set_edgecolor(src_color)
        ax.add_patch(p)
    for p in patch_list_bkg:
        p.set_edgecolor(bkg_color)
        ax.add_patch(p)
    ax.set_xlim(center_x - size / 2., center_x + size / 2.)
    ax.set_ylim(center_y - size / 2., center_y + size / 2.)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    fig.canvas.draw()
    fig.subplots_adjust(top = fig_borders[0], bottom = fig_borders[1],
                        left = fig_borders[2], right = fig_borders[3])

    return fig
