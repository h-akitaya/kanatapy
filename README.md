# Name

kanatapy

# Features

Kanata reduction and operation tools.

# Requirement

* astropy
* numpy
* scipy

# Installation

```bash
$ git clone https://github.com/h-akitaya/kanatapy.git
$ cd kanatapy
$ pip install --user kanatapy
```
# Usage

* hpkossub.py
```
 (path)/kanatapy/kanatapy/ccd/hpkossub.py
 usage: hpkossub [-h] [-s str] [-m] [-o] fn [fn ...]

positional arguments:
  fn                    Fits file names ('@list.txt': file list))

optional arguments:
  -h, --help            show this help message and exit
  -s str, --sub-extention str
                        Sub-extention of output file (default: .bs)
  -m, --median          Median value is uses as an overscan region level. (Default: Fitting overscan region in y-direction by Legendre 2nd func.
  -o, --overwrite       Overwrite existing file(s)
```


# Author

* Hiroshi AKITAYA
* PERC, Chiba Institute of Technology
* akitaya _at_ perc.it-chiba.ac.jp


