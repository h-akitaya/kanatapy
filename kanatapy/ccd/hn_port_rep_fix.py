#!/usr/bin/env python3
#
# HONIR CCD port replacement fix tool.
#   (for the case of one redundant pixel followed by real pixel data series)
#
#  Since 2022-11-24 H. Akitaya (PERC/CIT)
##

""" Revision history.
Ver. 0.1: 2022-11-24 H. Akitaya
Ver. 0.2: 2023-01-24 H. Akitaya; add argparse options.
Ver. 0.3: 2023-04-05 H. Akitaya; little fix.
"""

__version__ = '0.3'
__author__ = 'Hiroshi Akitaya'

import os
import sys

import numpy as np
from astropy.io import fits


def port_replacement_fix(fn_in, fn_out, overwrite=False):
    """
    Fix port replacement of an HPK CCD image.
    fn_in: input fits file.
    fn_out: output fits file.
    """
    # Output file and overwrite mode check.
    if os.path.exists(fn_out) and (overwrite is False):
        sys.stderr.write('File {} exists. Abort.\n'.format(fn_out))
        sys.exit(1)

    # Image reconstruction.
    with fits.open(fn_in) as hdul:
        imgdata = hdul[0].data
        hdr = hdul[0].header

        # Port 4 (Replacement, reflect, and 1-pix shift)
        data_p4_a = imgdata[:, 0:1]
        data_p4_b = imgdata[:, 1:536]
        data_p4 = np.fliplr(np.hstack((data_p4_b, data_p4_a)))

        # Port 1 (Replacement and reflect)
        data_p1_a = imgdata[:, 536:1072]
        data_p1 = np.fliplr(data_p1_a)

        # Port 2 (Replacement and reflect)
        data_p2_a = imgdata[:, 1072:1608]
        data_p2 = np.fliplr(data_p2_a)

        # Port 3 (Replacement and reflect)
        data_p3_a = imgdata[:, 1608:2144]
        data_p3 = np.fliplr(data_p3_a)

        # Merge
        hdul[0].data = np.hstack((data_p1, data_p2, data_p3, data_p4))

        # History.
        hdul[0].header['HISTORY'] = 'Port_replacement_fix done.'

        # Fits file write.
        hdul.writeto(fn_out, overwrite=overwrite)
        return hdul[0].data


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fn_in', type=str, help='Input file name (fits)')
    parser.add_argument('fn_out', type=str, help='Output file name (fits)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite mode')
    parser.add_argument('--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__),
                        help='Show version')
    args = parser.parse_args()

    file_name_in = args.fn_in
    file_name_out = args.fn_out

    if not os.path.isfile(file_name_in):
        print('File {} not found.'.format(file_name_in))
        sys.exit(1)

    if args.verbose:
        print('Overwrite mode {}'.format(args.overwrite))

    # Main routine.
    data = port_replacement_fix(file_name_in, file_name_out, overwrite=args.overwrite)

    if args.verbose:
        print('Fix. {} -> {}'.format(file_name_in, file_name_out))
