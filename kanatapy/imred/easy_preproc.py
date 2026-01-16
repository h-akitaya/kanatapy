#
# quick-look preprocess for HONIR
#    since 2024-04-18 H. Akitaya (ARC/CIT)
#
import sys

import numpy as np

from astropy.io import fits

from kanatapy.imred import redtools

class EasyPreProc(object):

    def __init__(self, fn):
        try:
            self.hdul = fits.open(fn)
        except IOError:
            sys.stderr.write('Cannot open file {}\n'.format(fn))
            sys.exit(1)

    def __del__(self):
        self.hdul.close()

    def do_preproc(self, fn):
        self.trim_image()
        self.flatten_image()
        self.wcs_resolve()
    def trim_image(self, fn):
        if self.hdul[0].header['HN-ARM'] == 'opt':
            pass

    def flatten_image(self):
        pass

    def wcs_resolve(self):
        pass
