# HONIR reduction tools.
# since 2024/8/6 H. Akitaya

__version__ = "0.1.1"
__author__ = "Akitaya, Hiroshi"
__email__ = "akitaya@perc.it-chiba.ac.jp"

import glob
import sys
import inspect

from astropy.io import fits

from kanatapy.imred import redtools as rt
from kanatapy.ccd.hpkossub import HPKOsSub


def port_replacement_fix(fn: str, fn_out: str, overwrite: bool = True):
    """Wrapper function for kanatapy.ccd.hn_port_rep_fix.port_replacement_fix().
    :param str fn: input file name
    :param str fn_out: output file name
    :param bool overwrite: overwrite mode
    """
    from kanatapy.ccd.hn_port_rep_fix import port_replacement_fix
    port_replacement_fix(fn, fn_out, overwrite=overwrite)


def ossub(fn, compat_hntrimccd=True, overwrite=True):
    """HONIR CCD image overscan region removal.
    :param str fn: input file name
    :param bool compat_hntrimccd: compatibility to the iraf hntrimccd
    :param bool overwrite: overwrite mode
    :return: output file name
    """
    hnccd = HPKOsSub(fn, howpol=False, compat_hntrimccd=compat_hntrimccd, overwrite=overwrite)
    fn_out = hnccd.ossub_all()
    print(fn_out)
    return fn_out


def flat_fielding(fn_in: str, fn_flat: str, fn_out: str = '', overwrite: bool = False):
    """
    Bias subtraction.
    :param fn_in: Input filename.
    :param fn_flat: Flat filename.
    :param fn_out: Output filename.
    :param overwrite: Overwrite mode.
    :return:
    """
    hdul_in = fits.open(fn_in)
    hdul_flat = fits.open(fn_flat)

    # Check ccd-mode. (for CCD)
    if hdul_in[0].header['HN-ARM'] == 'opt':
        ccd_mode_in = hdul_in[0].header['CCD-MODE']
        ccd_mode_flat = hdul_flat[0].header['CCD-MODE']
        if ccd_mode_in != ccd_mode_flat:
            sys.stderr.write(f'CCD-MODE differs: {ccd_mode_in} {ccd_mode_flat}\n')
            sys.exit(1)

    if fn_out == '':
        fn_out = rt.get_filename_with_subsuffix(fn_in, '_fl')

    hdul_in[0].data /= hdul_flat[0].data

    hdul_in[0].header['HISTORY'] = (f'{__name__}: {__version__}, '
                                    f'{inspect.currentframe().f_code.co_name}')
    hdul_in[0].header['HISTORY'] = f'Flat fielded.'
    hdul_in[0].header['HISTORY'] = f'{rt.get_now_datetime_str()}'
    hdul_in[0].header['HISTORY'] = f'{fn_in} - {fn_flat} -> {fn_out}.'

    try:
        hdul_in.writeto(fn_out, overwrite=overwrite)
    except IOError:
        sys.stderr.write(f'Image write error: {fn_out}')
    return fn_out


def bias_subtract(fn_in: str, fn_bias: str, fn_out: str = '', overwrite: bool = False):
    """
    Bias subtraction.
    :param fn_in: Input filename.
    :param fn_bias: Bias filename.
    :param fn_out: Output filename.
    :param overwrite: Overwrite mode.
    :return:
    """
    hdul_in = fits.open(fn_in)
    hdul_bias = fits.open(fn_bias)

    # Check ccd-mode.
    ccd_mode_in = hdul_in[0].header['CCD-MODE']
    ccd_mode_bias = hdul_bias[0].header['CCD-MODE']
    if ccd_mode_in != ccd_mode_bias:
        sys.stderr.write(f'CCD-MODE differs: {ccd_mode_in} {ccd_mode_bias}\n')
        sys.exit(1)

    if fn_out == '':
        fn_out = rt.get_filename_with_subsuffix(fn_in, '_bs')

    hdul_in[0].data -= hdul_bias[0].data
    hdul_in[0].header['HISTORY'] = (f'{__name__}: {__version__}, '
                                    f'{inspect.currentframe().f_code.co_name}')
    hdul_in[0].header['HISTORY'] = f'Bias subtracted.'
    hdul_in[0].header['HISTORY'] = f'{rt.get_now_datetime_str()}'
    hdul_in[0].header['HISTORY'] = f'{fn_in} - {fn_bias} -> {fn_out}.'

    try:
        hdul_in.writeto(fn_out, overwrite=overwrite)
    except IOError:
        sys.stderr.write(f'Image write error: {fn_out}')
    return fn_out


def select_bias_images(fn_in, ccdmode: str = None):
    """Select bias images.
    :param fn_in: file name pattern or list.
    :param ccdmode: CCD readout mode (default: ignore mode.
        full: Fill readout, partial: Partial readout.
    :return:
    """

    if type(fn_in) is not list:
        fn_list = glob.glob(fn_in)
        fn_list.sort()
    else:
        fn_list = fn_in

    filter = {'DATA-TYP': 'BIAS',
              'OBJECT': 'BIAS'}
    if ccdmode is not None:
        if ccdmode.lower() == 'full':
            filter['CCD-MODE'] = 'Full'
        elif ccdmode.lower() == 'partial':
            filter['CCD-MODE'] = 'Partial'
    fns_bias = rt.select_images_by_header(fn_list, filter, logic='and')
    return fns_bias
