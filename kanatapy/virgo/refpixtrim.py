#!/usr/bin/env python3
"""
Kanata VIRGO tool for removing reference pixel regions.
"""

#    VIRGO REDUCTION SUB ROUTINE
#       refpixtrim.py: : Reference pixel region subtraction and merge 4 port images
#
#       Original code:
#           hpkossub.py and hntrimvirgo.cl by H. Akitaya
#
#       Python version:
#         Ver. 1.0  2023-03-30: H. Akitaya
#         Ver. 1.1  2024-08-07 H. Akitaya; comment modified.


import os
from datetime import datetime
import argparse

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy import stats

import astropy.io.fits as fits


class RefPixTrim(object):
    """
    HONIR Virgo reference pixel regions subtraction class.
    """
    __N_PORTS = 4  # Number of readout ports.

    def __init__(self, _fn, verbose=True, sub_extension='_bt', overwrite=False,
                 debug=False) -> None:
        """
        Initialization.
        :rtype: None
        :param str _fn: input filename
        :param bool verbose: verbose mode
        :param str sub_extension: sub-extension of output file. (a.fits -> a_subext.fits)
        :param bool overwrite: overwrite mode
        :param bool debug: debug mode
        """
        self.fn = _fn  # Input file name.
        self.hdulist = None  # HDU list for the input imred.
        self.hdu = None  # The first HDU of the input imred.
        self.verbose = verbose  # Verbose mode (bool).
        self.y1 = None  # Y-start index of the effective area.
        self.y2 = None  # Y-end index of the effective area.
        self.sub_extension = sub_extension  # Sub-extention (xxxxSUBEXT.ext).
        self.overwrite = overwrite  # Overwrite mode (bool).
        self.stat = {'stddevs': [], 'means': []}  # Statistics tables.
        self.fn_out = self.get_subext_fn(_fn)  # Output file name.
        self.debug = debug  # Debug mode.

    def set_fn(self, _fn: str) -> None:
        """
        Set file name.
        :param _fn: Input file name. (e.g.) HN0012345opt00.fits
        """
        self.fn = _fn  # Input file name.
        self.fn_out = self.get_subext_fn(_fn)  # Output file name.

    def read_image(self) -> None:
        """ Read fits imred.
        """
        try:
            self.hdulist = fits.open(self.fn)
        except OSError:
            raise OSError('File {} read error.\n'.format(self.fn))

        if self.verbose:
            print('Fits file {} opened.'.format(self.fn))
        self.hdu = self.hdulist[0]

    def check_processed(self) -> bool:
        """ Check the imred was processed or not (check fits header REFPIXTR).
        """
        if 'REFPIXTR' in self.hdu.header:
            if self.hdu.header['REFPIXTR'] is True:
                return True
        return False

    def read_yrange(self) -> None:
        """ Read size of y-range from fits headers.
        """
        self.y1 = self.hdu.header['EFPYMIN1'] - 1
        self.y2 = self.y1 + self.hdu.header['EFPYRNG1'] - 1

    def get_port_area_data(self, port_n: int) -> np.ndarray:
        """
        Get a port data array of number port_n.
        """
        if port_n < 1 or 4 < port_n:
            raise ValueError('port_ns should be 1-4. {} given.')
        # Define fits header keys of coordinates of the port range.
        kyd_img_x1 = 'PORT{:1d}X1'.format(port_n)
        kyd_img_x2 = 'PORT{:1d}X2'.format(port_n)

        # Index numbers of the imred area,
        x1 = self.hdu.header[kyd_img_x1] - 1
        x2 = self.hdu.header[kyd_img_x2] - 1

        # Cut out the imred area as ndarray.
        img = self.hdu.data[self.y1:self.y2 + 1, x1:x2 + 1]

        # Return imred area ndarray of the selected port.
        return img

    def get_effective_area_data(self, port_n: int) -> np.ndarray:
        """ Get data of an effective area of port # port_n.
        """
        if port_n < 1 or 4 < port_n:
            raise ValueError('port_ns should be 1-4. {} given.')

        # Define fits header keys of coordinates of the effective range
        # in the port no. port_n.
        kyd_eff_x1 = 'EFPXMIN{:1d}'.format(port_n)
        kyd_eff_rng = 'EFPXRNG{:1d}'.format(port_n)

        # Index numbers of the effective imred area,
        x1 = self.hdu.header[kyd_eff_x1] - 1
        x_rng = self.hdu.header[kyd_eff_rng]

        # Cut out the imred area as ndarray.
        img_eff = self.hdu.data[self.y1:self.y2 + 1, x1:x1 + x_rng]

        # Return imred area ndarray of the selected port.
        return img_eff

    def get_psos_data(self, img: np.ndarray, port_n: int) -> (np.ndarray, np.ndarray, np.ndarray):
        """ Get pre/overscan regions data array.
        """
        if port_n < 1 or 4 < port_n:
            raise ValueError('port_ns should be 1-4. {} given.')

        # Odd port number.
        if port_n == 1 or port_n == 3:
            plabel = '13'
        # Even port number.
        elif port_n == 2 or port_n == 4:
            plabel = '24'
        else:
            plabel = '13'

        # Define fits header keys of pre/overscan regions.
        kyd_ps_x1 = 'PSRE{:2s}X1'.format(plabel)
        kyd_ps_x2 = 'PSRE{:2s}X2'.format(plabel)
        kyd_os_x1 = 'OSRE{:2s}X1'.format(plabel)
        kyd_os_x2 = 'OSRE{:2s}X2'.format(plabel)

        # Index numbers of pre/0verscan regions.
        ps_x1 = self.hdu.header[kyd_ps_x1] - 1
        ps_x2 = self.hdu.header[kyd_ps_x2] - 1
        os_x1 = self.hdu.header[kyd_os_x1] - 1
        os_x2 = self.hdu.header[kyd_os_x2] - 1

        # Cut out pre/overscan regions array
        img_ps = img[0:self.y2 - self.y1 + 1, ps_x1:(ps_x2 + 1)]
        img_os = img[0:self.y2 - self.y1 + 1, os_x1:(os_x2 + 1)]

        # Concatenate pre/overscan regions
        img_psos = np.concatenate((img_ps, img_os), axis=1)

        # Return concatenated array, prescan region array, and ovserscan array.
        return img_psos, img_ps, img_os

    @staticmethod
    def calc_stddev_and_mean(img: np.ndarray,
                             low: float = 2.5, high: float = 2.5) -> (float, float):
        """ Calculate a standard deviation and mean in an imred array 'img' with sigma-clipping method.
        """
        # Get sigma-clipped array.
        clp, _, _ = stats.sigmaclip(np.ravel(img), low, high)
        # Calculate a standard deviation.
        stddev = np.std(clp)
        mean = np.mean(clp)
        return stddev, mean

    def subtract_overscan_region(self, img_eff: np.ndarray, img_psos: np.ndarray,
                                 median=True) -> np.ndarray:
        """ Subtract overscan region.
        If median is False, fit overscan region values in y-direction with
        a polynomial function. Otherwise, median value is subtracted.
        """
        # Median subtraction. (fit == False)
        if median is True:
            midpt = np.median(img_psos, axis=None)
            img_sub = img_eff - midpt
            return img_sub
        else:
            return img_eff

        """ # function fitting used for HPK CCD tool.
        # Function fitting subtraction. (fit == True)
        else:
            psosr_1d = np.median(img_psos, axis=1)
            #  print(img_eff.shape, len(psosr_1d))  # Debug.
            x = np.arange(0, self.y2 - self.y1 + 1, 1)
            # Fitting Legendre 2nd order function.
            ospsr_fit = Legendre.fit(x, psosr_1d, 2)
            # A y-direction discrete values from the fitting function.
            fit_array = ospsr_fit.linspace(self.y2 - self.y1 + 1)

            # Debug: plot fitting
            if self.debug:
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(x, psosr_1d, marker='o', label='original data')
                ax.plot(x, fit_array[1], label='fitting')
                plt.legend()
                fig.savefig('img.png'.format())

            # Make 2d array from the fitting function.
            fit_array_2d = np.tile(fit_array[1], (img_eff.shape[1], 1)).T
            img_sub = img_eff - fit_array_2d
            return img_sub
            """

    def get_subext_fn(self, _fn: str) -> str:
        """ Get filename with sub-extension.
        :param _fn: File name.
        :return: File name with a sub-extension.
        """
        fns = os.path.splitext(_fn)
        fn_out = fns[0] + self.sub_extension + fns[1]
        return fn_out

    def write_file(self, img_bs: np.ndarray) -> None:
        """ Write result to a new fits file.
        """
        hdu_out = fits.PrimaryHDU(img_bs.astype('float32'),
                                  header=self.hdu.header)
        # hdu_out.header = hs.hdu.header
        hdul_out = fits.HDUList([hdu_out])
        self.append_log_to_fits_headers(hdul_out[0].header)
        if os.path.exists(self.fn_out) and self.overwrite is True:
            try:
                os.remove(self.fn_out)
            except PermissionError:
                raise PermissionError('Failed to overwrite file: {}\n'.format(self.fn_out))
            print('File {} exists. Overwrite.'.format(self.fn_out))
        try:
            hdul_out.writeto(self.fn_out)
        except OSError:
            raise OSError('Failed to write file: {}\n'.format(self.fn_out))
        if self.verbose:
            print('Fits file {} written.'.format(self.fn_out))

    def append_log_to_fits_headers(self, hdr: fits.header.Header) -> None:
        # History.
        hdr['history'] = '{}: {}'.format(
            os.path.basename(__file__),
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        )
        hdr['history'] = 'File names: {} -> {}'.format(self.fn, self.fn_out)
        # Readout noise.
        for port_n in range(1, len(self.stat['stddevs']) + 1):
            hdr['RFRNO_P{:1d}'.format(port_n)] = \
                (self.stat['stddevs'][port_n - 1],
                 'Reference pixel reg. stddev. (port #{:1d})'.format(port_n))
        # Mean.
        for port_n in range(1, len(self.stat['means']) + 1):
            hdr['RFAVE_P{:1d}'.format(port_n)] = \
                (self.stat['means'][port_n - 1],
                 'Reference pixel reg. mean (port #{:1d})'.format(port_n))
        # Check flag of processing.
        hdr['REFPIXTR'] = (True, 'Reference pixel region subtraction (Boolean)')

    def refpixtrim_all(self, median: bool = False) -> None:
        """ Do all processes of reference pixel region subtraction.
        """

        # Check output file existence.
        if os.path.exists(self.fn_out) and self.overwrite is False:
            print('File {} exists. Skip.'.format(self.fn_out))
            raise FileExistsError('File {} exists. Skip.'.format(self.fn_out))

        # File open.
        self.read_image()

        # Check the imred has been processed.
        if self.check_processed() is True:
            raise ValueError('Image has been already processed. Skip.')

        # Read y-range information from the original fits file.
        self.read_yrange()
        # Memory array for processed sub-images.
        img_subs = []

        # Process each port area.
        for port_n in range(1, RefPixTrim.__N_PORTS + 1):
            # Cut out the entire area of the port #port_n.
            img = self.get_port_area_data(port_n)
            # Cut out the effective area.
            img_eff = self.get_effective_area_data(port_n)
            # Cut out the pre-scan and over-scan regions.
            img_psos, img_ps, img_os = self.get_psos_data(img, port_n)
            # Calculate statistics.
            stddev, mean = self.calc_stddev_and_mean(img_os)
            self.stat['stddevs'].append(stddev)
            self.stat['means'].append(mean)
            # Subtract overscan level from the sub-imred and
            # store it to the memory array.
            img_subs.append(self.subtract_overscan_region(img_eff, img_os, median=median))

        # Concatenate all the processed sub-images.
        img_sub = np.concatenate(img_subs, axis=1)

        self.write_file(img_sub)


if __name__ == '__main__':

    # Option parser
    parser = argparse.ArgumentParser(prog='refpixtrim', fromfile_prefix_chars='@')
    parser.add_argument('-s', '--sub-extention', metavar='str', type=str,
                        default='_bt',
                        help='Sub-extention of output file (default: _bt)')
    parser.add_argument('-m', '--median', action='store_true',
                        default=False,
                        help='Median value is used as an reference pixel region level. '
                        'Default: not subtract a level of reference pixel region.')
#                             'Default: Fitting overscan region in y-direction by Legendre 2nd func.')
    parser.add_argument('-o', '--overwrite', action='store_true',
                        default=False,
                        help='Overwrite existing file(s)')
    parser.add_argument('files', metavar='fn', type=str, nargs='+',
                        help='Fits file names (\'@list.txt\': file list))')
    args = parser.parse_args()

    for fn in args.files:
        hs = RefPixTrim(fn, sub_extension=args.sub_extention, overwrite=args.overwrite)
        # Subtract overscan regions.
        hs.refpixtrim_all(args.median)
