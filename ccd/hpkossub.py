#!/usr/bin/env python3
#
#    HPK CCD REDUCTION SUB ROUTNINE
#       hpbossub.py: : Overscan subtraction and Merge 4 port images
#
#       Original code (howossub.cl: CL script) by K. S. Kawabata
#
#       Python version:
#         Ver 1.0  2021-10-26: H. Akitaya
#

import sys
import os
from datetime import datetime
import argparse

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy import stats

import astropy.io.fits as fits

class HPKOsSub(object):
    
    def __init__(self, fn, verbose=True, sub_extention='.bs', overwrite=False,
                 howpol=True):
        self.fn = fn  # Input file name.
        self.hdulist = None  # HDU list for the input image.
        self.hdu = None  # The first HDU of the input image.
        self.verbose = verbose  # Verbose mode (bool).
        self.y1 = None  # Y-start index of the effective area.
        self.y2 = None  # Y-end index of the effective area.
        self.sub_extention = sub_extention  # Sub-extention (xxxxSUBEXT.ext).
        self.overwrite = overwrite  # Overwrite mode (bool).
        self.howpol = howpol  # Instrument is HOWPol or not. (Default: HOWPol)
        self.stat = {'stddevs': [], 'means': []}  # Statistics tables.
        self.fn_out = self.get_subext_fn(fn)  # Output file name.

    def read_image(self):
        """ Read fits image.
        """
        try:
            self.hdulist = fits.open(self.fn)
        except OSError:
            sys.stderr.write('File {} read error.\n'.format(self.fn))
            return False
        if self.verbose:
            print('Fits file {} open.'.format(self.fn))
        self.hdu = self.hdulist[0]
        return True

    def check_processed(self):
        """ Check the image was processed or not.
        """
        if 'REDOSSUB' in self.hdu.header:
            if self.hdu.header['REDOSSUB'] is True:
                return True
        return False

    def read_yrange(self):
        """ Read grid numbers of yrange from fits headers.
        """
        self.y1 = self.hdu.header['EFPYMIN1'] - 1
        self.y2 = self.y1 + self.hdu.header['EFPYRNG1'] - 1

    def get_port_area_data(self, port_n=1):
        """ Get a port data array of number port_n.
        """
        # Define fits header keys of coordinates of the port range.
        kyd_img_x1 = 'PORT{:1d}X1'.format(port_n)
        kyd_img_x2 = 'PORT{:1d}X2'.format(port_n)

        # Index numbers of the image area,
        x1 = self.hdu.header[kyd_img_x1] - 1
        x2 = self.hdu.header[kyd_img_x2] - 1

        # Cut out the image area as ndarray.
        img = self.hdu.data[self.y1:self.y2+1, x1:x2+1]

        # Return image area ndarray of the selected port.
        return img

    def get_effective_area_data(self, port_n=1):
        """ Get data of an effective area of port # port_n.
        """
        # Define fits header keys of coordinates of the effective range
        # in the port no. port_n.
        kyd_eff_x1 = 'EFPXMIN{:1d}'.format(port_n)
        kyd_eff_rng = 'EFPXRNG{:1d}'.format(port_n)
        
        # Index numbers of the effective image area,
        x1 = self.hdu.header[kyd_eff_x1] - 1
        x_rng = self.hdu.header[kyd_eff_rng]

        # Cut out the image area as ndarray.
        img_eff = self.hdu.data[self.y1:self.y2+1, x1:x1+x_rng]

        # Return image area ndarray of the selected port.
        return img_eff

    def get_psos_data(self, img, port_n):
        """ Get pre/overscan regions data array.
        """
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
        img_ps = img[0:self.y2-self.y1+1, ps_x1:(ps_x2+1)]
        img_os = img[0:self.y2-self.y1+1, os_x1:(os_x2+1)]

        # Concantenate pre/overscan regions
        img_psos = np.concatenate((img_ps, img_os), axis=1)

        # Return concatenated array, prescan region array, and ovserscan array.
        return img_psos, img_ps, img_os

    def calc_stddev_and_mean(self, img, low=2.5, high=2.5):
        """ Calculate a standard deviation and mean in an image array 'img' with
        sigma-clipping methond.
        """
        # Get sigma-clipped array.
        clp, _, _ = stats.sigmaclip(np.ravel(img), low, high)
        # Calculate a standard deviation.
        stddev = np.std(clp)
        mean = np.mean(clp)
        return stddev, mean

    def subtract_overscan_region(self, img_eff, img_psos, median=True,
                                 plot=False):
        """ Subtract overscan region.
        If median is False, fit overscan region values in y-direction with
        a polynominal function. Otherwise, median value is subtracted.
        """
        # Median subtraction. (fit == False)
        if median is True:
            midpt = np.median(img_psos, axis=None)
            img_sub = img_eff - midpt
            return img_sub
        # Function fitting subtraction. (fit == True)
        else:
            psosr_1d = np.median(img_psos, axis=1)
            #  print(img_eff.shape, len(psosr_1d))  # Debug.
            x = np.arange(0, self.y2-self.y1+1, 1)
            # Fitting Legendre 2nd order function.
            ospsr_fit = Legendre.fit(x, psosr_1d, 2)
            # A y-direction discrete values from the fitting function.
            fit_array = ospsr_fit.linspace(self.y2-self.y1+1)

            # Debug: plot fitting
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

    def get_subext_fn(self, fn):
        """ Get filename with sub-extension.
        """
        fns = os.path.splitext(self.fn)
        fn_out = fns[0] + self.sub_extention + fns[1]
        return fn_out

    def write_file(self, img_bs):
        """ Write result to new fits file.
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
                sys.stderr.write('Failed to overwrite file: {}\n'.format(self.fn_out))
                return False
            print('File {} exists. Overwrite.'.format(self.fn_out))
        try:
            hdul_out.writeto(self.fn_out)
        except OSError:
            sys.stderr.write('Failed to write file: {}\n'.format(self.fn_out))
            return False
        return True

    def append_log_to_fits_headers(self, hdr):
        # History.
        hdr['history'] = '{}: {}'.format(
            os.path.basename(__file__),
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        )
        hdr['history'] = 'File names: {} -> {}'.format(self.fn, self.fn_out)
        # Readout noise.
        for port_n in range(1, len(self.stat['stddevs'])+1):
            hdr['OSRNO_P{:1d}'.format(port_n)] = \
                (self.stat['stddevs'][port_n-1],
                 'Pre/Ovserscan reg. stddev. (port #{:1d})'.format(port_n))
        # Mean.
        for port_n in range(1, len(self.stat['means'])+1):
            hdr['OSAVE_P{:1d}'.format(port_n)] = \
                (self.stat['means'][port_n-1],
                 'Pre/Ovserscan reg. mean (port #{:1d})'.format(port_n))
        # Check flag of processing.
        hdr['REDOSSUB'] = (True, 'Pre/over-scan region subtraction (Boolean)')


if __name__ == '__main__':
    
    # Option parser
    parser = argparse.ArgumentParser(prog='hpkossub', fromfile_prefix_chars='@')
    #parser.add_argument('-l', '--list', metavar='str', type=str,
    #                       default='HOW.list',
    #                       help = 'File list name')
    parser.add_argument('-s', '--sub-extention', metavar='str', type=str,
                           default='.bs',
                           help = 'Sub-extention of output file (default: .bs)')
    parser.add_argument('-m', '--median', action='store_true',
                        default=False,
                        help = 'Median value is uses as an overscan region level. (Default: Fitting overscan region in y-direction by Legendre 2nd func.')
    parser.add_argument('-o', '--overwrite', action='store_true', 
                        default=False,
                        help = 'Overwrite existing file(s)')
    parser.add_argument('files', metavar='fn', type=str, nargs='+',
                        help='Fits file names (\'@list.txt\': file list))')
    args = parser.parse_args()

    """ Debug.
    """
    #print(args)  # for debug
    #print(args.files)  # for debug

    for fn in args.files:
        hs = HPKOsSub(fn, sub_extention=args.sub_extention, overwrite=args.overwrite)
        # Check output file existence.
        if os.path.exists(hs.fn_out) and args.overwrite is False:
            print('File {} exists. Skip.'.format(hs.fn_out))
            continue
        # File open.
        if hs.read_image() is False:
            continue
        # Check the image has been processed.
        if hs.check_processed() is True:
            print('Image has been already processed. Skip.')
            continue
        hs.read_yrange()
        img_subs = []
        for port_n in range(1, 5):
            img = hs.get_port_area_data(port_n)
            img_eff = hs.get_effective_area_data(port_n)
            img_psos, img_ps, img_os = hs.get_psos_data(img, port_n)
            #stddev, mean = hs.calc_stddev_and_mean(img_psos)
            stddev, mean = hs.calc_stddev_and_mean(img_os)
            hs.stat['stddevs'].append(stddev)
            hs.stat['means'].append(mean)
            #img_subs.append(hs.subtract_overscan_region(img_eff, img_psos, fit=fit))
            img_subs.append(hs.subtract_overscan_region(img_eff, img_os, median=args.median))
        
        img_sub = np.concatenate(img_subs, axis=1)
        hs.write_file(img_sub)
