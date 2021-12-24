#!/usr/bin/env python
"""
Sample code for HPKOsSub class.
2021-12-24  H. Akitaya
"""

import sys
from kanatapy.ccd import hpkossub

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: {} filename'.format(sys.argv[0]))
        sys.exit(1)

    fn = sys.argv[1]
    hs = hpkossub.HPKOsSub(fn)
    hs.ossub_all(median=False)
