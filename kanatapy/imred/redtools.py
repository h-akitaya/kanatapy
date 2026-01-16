# Image reduction tools.
# since 2024/08/06

__version__ = "0.1.1"
__author__ = "Akitaya, Hiroshi"
__email__ = "akitaya@perc.it-chiba.ac.jp"

import sys
import os
from datetime import datetime

import numpy as np
from astropy.io import fits
from scipy.constants import hp


def mark_ignore_flag(fn: str, remove_flag: bool = False):
    """
    Mark or unmark an ignore flag on a image. (REDIGNR -> Ture or False)
    :param fn: file name
    :bool unmark: unmark a flag. (REDIGNR -> False)
    :return:
    """
    hdul = fits.open(fn, mode='update')
    if remove_flag is False:
        flag = True
    else:
        flag = False
    hdul[0].header['REDIGNR'] = flag
    hdul.flush()
    hdul.close()
    return flag

def imcombine(fn_list: list, fn_out: str = "", overwrite: bool = True,
              mode: str = 'median', verbose=False) -> fits.HDUList:
    """ Combine fits images.
    median mode.
    :param list fn_list: file name list to be combined
    :param str fn_out: output file name
    :param bool overwrite: overwrite mode
    :param bool verbose: verbose mode
    :param str mode: combine mode (median, average)
    :return: HDUList.
    """

    # Read image data.
    imgs = []
    for fn in fn_list:
        with fits.open(fn) as hdul:
            imgs.append(hdul[0].data)
            if verbose:
                print(fn)
    if len(imgs) < 1:
        sys.stderr.write('No images found.\n')
        sys.exit(1)
    imgs_ndarray = np.stack(imgs)
    del imgs

    # Combine images.
    if mode == 'median':
        imgs_comb = np.median(imgs_ndarray, axis=0)
    else:
        imgs_comb = np.average(imgs_ndarray, axis=0)

    del imgs_ndarray

    # Write history into the fits header.
    hdul[0].data = imgs_comb
    hdul[0].header['history'] = f'imcombine done.'
    hdul[0].header['history'] = f'Combined: {str(fn_list)}'
    hdul[0].header['history'] = f'Combined file name: {fn_out}'

    # Write to a fits file.
    if fn_out != "":
        hdul.writeto(fn_out, overwrite=overwrite)
        print(f'Write to {fn_out}.')

    return hdul


def subtract_bias(fn: str, fn_bias: str, overwrite: bool = False,
                  subext: str = '_bs') -> str:
    """ Bias subtraction.
    :param str fn: objective file name
    :param str fn_bias: bias file name
    :param bool overwrite: overwrite mode
    :param str subext: sub-extension of the output file name
    :return: output file name.
    """
    try:
        hdul_img = fits.open(fn)
        hdul_bias = fits.open(fn)
    except IOError:
        sys.stderr.write('File open error.')
        return ""

    fn_splt = os.path.splitext(fn)
    fn_new = fn_splt[0] + subext + fn_splt[1]
    hdul_img[0].data -= hdul_bias[0].data

    # Write history.
    hdul_img[0].header['history'] = f'bias_subtraction'
    hdul_img[0].header['history'] = f'Bias subtraction: {fn} - {fn_bias} = {fn_new}'
    hdul_img.writeto(fn_new, overwrite=overwrite)

    hdul_img.close()
    hdul_bias.close()

    return fn_new


def group_by_hwpangle(fns: list):
    """ Group images by HWPANGLE
    :param list fns: file name list to be grouped
    :return: dictionary of the file list."""
    hwpangle_dict = {}
    for fn in fns:
        with fits.open(fn) as hdul:
            if not hdul[0].header.has_key('HWPANGLE'):
                continue
            hwpangle = hdul[0].header['HWPANGLE']
            if hwpangle not in hwpangle_dict:
                hwpangle_dict[hwpangle] = []
            hwpangle_dict[hwpangle].append(fn)
    return hwpangle_dict

def create_flat_image(flat_fns: list, bias_fn: str, flat_fn: str,
                      overwrite: bool = True, mode: str = 'median') -> fits.HDUList:
    """ Create flat image.
    :param list flat_fns: list of source images
    :param str bias_fn: bias file name
    :param str flat_fn: output fiat file name
    :param bool overwrite: overwrite mode
    :param str mode: combine mode (median, average)
    :return: HDUList of the result
    """
    # Open files.
    img_flat_comb_hdul = imcombine(flat_fns, fn_out="", mode=mode)
    try:
        hdul_bias = fits.open(bias_fn)
    except IOError:
        sys.stderr.write(f'Bias file {bias_fn} open error.\n')
        sys.exit(1)

    # Calculations.
    # Bias subtraction.
    img_flat_bs = img_flat_comb_hdul[0].data - hdul_bias[0].data

    # Normalization.
    flat_median_val = np.median(img_flat_bs)
    img_flat_nrm = img_flat_bs / flat_median_val
    img_flat_comb_hdul[0].data = img_flat_nrm

    # Write history.
    img_flat_comb_hdul[0].header['history'] = f'imcombine'
    img_flat_comb_hdul[0].header['history'] = f'Combined: {str(flat_fns)}'
    img_flat_comb_hdul[0].header['history'] = f'Bias subtracted: {bias_fn}'
    img_flat_comb_hdul[0].header['history'] = f'Normalized: {flat_median_val}'
    img_flat_comb_hdul[0].header['history'] = f'File name: {flat_fn}'
    img_flat_comb_hdul.writeto(flat_fn, overwrite=overwrite)
    print(f'Write to {flat_fn}')

    hdul_bias.close()
    img_flat_comb_hdul.close()
    return hdul_bias


def select_images_by_header(fns: list, cond_dict: dict, logic: str = 'or') -> list:
    """ Select filenames by fits header values
    :param list fns: file name list to be searched
    :param dict cond_dict: dictionary of a pair of header_key and value
    :param str logic: logic
    :return: file name list
    """
    fns_result = []
    for fn in fns:
        with fits.open(fn) as hdul:
            if logic == 'or':
                for key in cond_dict:
                    if hdul[0].header[key] == cond_dict[key]:
                        fns_result.append(fn)
            else:
                flag = True
                for key in cond_dict:
                    if hdul[0].header[key] != cond_dict[key]:
                        flag = False
                if flag:
                    fns_result.append(fn)

    return fns_result


def show_all_headers(fn: str) -> None:
    """ Show all fits header information
    :param str fn: input file name
    :return: none
    """
    with fits.open(fn) as hdul:
        for hdr, value in hdul[0].header.items():
            print(hdr, value)


def easy_imshow(fn: str, std_factor: float = 3, figsize: tuple = (5, 5)):
    """Easy image viewer.
    :param fn: file name
    :param std_factor: factor (stddev x factor) for image count range
    :param figsize: figure size in tuple
    """
    import matplotlib.pyplot as plt

    with fits.open(fn) as hdul:
        tmp_img = hdul[0].data
        med = np.median(tmp_img)
        std = np.std(tmp_img)
        plt.figure(figsize=figsize)
        plt.imshow(tmp_img, vmin=med - std_factor * std, vmax=med + std_factor * std)


def imstatistics(fn: str, show: bool = False):
    """IRAF imstatistics like tool.
    :param str fn: file name
    :param bool show: show results
    :return: statistic values; (average, median, stddev, max, min)
    """
    with fits.open(fn) as hdul:
        hdul_data = hdul[0].data
        (ave, midpt, stdev, maxval, minval) = get_imstatistics_ndarray(hdul_data)
        if show:
            print("average: {}".format(ave))
            print("median: {}".format(midpt))
            print("stdev: {}".format(stdev))
            print("max: {}".format(maxval))
            print("min: {}".format(minval))
        return ave, midpt, stdev, maxval, minval


def get_imstatistics_ndarray(hdul_data):
    ave = np.average(hdul_data)
    midpt = np.median(hdul_data)
    stdev = np.std(hdul_data)
    maxval = np.max(hdul_data)
    minval = np.min(hdul_data)
    return ave, midpt, stdev, maxval, minval


def get_filename_with_subsuffix(fn: str, sub_suffix: str):
    """
    Image name with sub-suffix.
    (e.g.) (XXXX).fits -> (XXXX)_(sub_suffix).fits
    :param fn: Input filename. (xxxxxxxxx.yyy)
    :param sub_suffix: Sub-suffix. (_zzz)
    :return: New filename. (xxxxxxxx_zzz.yyy)
    """
    fn_part = os.path.splitext(fn)
    return fn_part[0] + sub_suffix + fn_part[1]


def get_now_datetime_str():
    return datetime.now().strftime('%Y-%m-%dT%H:%M:%S')