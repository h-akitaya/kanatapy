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

## hpkossub.py
```
 (path)/kanatapy/kanatapy/ccd/hpkossub.py
usage: hpkossub [-h] [-s str] [-m] [-o] [-c] fn [fn ...]

positional arguments:
  fn                    Fits file names ('@list.txt': file list))

optional arguments:
  -h, --help            show this help message and exit
  -s str, --sub-extention str
                        Sub-extention of output file (default: .bs)
  -m, --median          Median value is uses as an overscan region level. Default: Fitting overscan region in y-direction by Legendre 2nd func.
  -o, --overwrite       Overwrite existing file(s)
  -c, --compat-hntrimccd
                        Compatible to hntrimccd.cl in honirred IRAF package.
```
### example 1
```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py HP0277904_0.fits
Fits file HP0277904_0.fits opened.
Fits file HP0277904_0.bs.fits written.
```
### example 2
```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -o HP0277904_0.fits
Fits file HP0277904_0.fits opened.
File HP0277904_0.bs.fits exists. Overwrite.
Fits file HP0277904_0.bs.fits written.
```
### example 3
```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py HP027790[12]_0.fits
Fits file HP0277901_0.fits opened.
Fits file HP0277901_0.bs.fits written.
Fits file HP0277902_0.fits opened.
Fits file HP0277902_0.bs.fits written.
```
### example 4
```
$ ls -1 HP027791[12]_0.fits > file.lst
$ cat file.lst 
HP0277911_0.fits
HP0277912_0.fits
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py @file.lst
Fits file HP0277911_0.fits opened.
Fits file HP0277911_0.bs.fits written.
Fits file HP0277912_0.fits opened.
Fits file HP0277912_0.bs.fits written.
```
### example 5
```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -s "_ossub" HP0277958_0.fits
Fits file HP0277958_0.fits opened.
Fits file HP0277958_0_ossub.fits written.
```
### example 6
Compatible mode to hntrimccd.cl in honirred cl-script IRAF package.
Output image with the same y-direction size as hntrimccd.cl.
```
# hntrimccd.cl compatible mode. (with option -c or --compat-hntrimcccd)
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -c -s '_bs' HN0123456opt00.fits
Fits file HN0123456opt00.fits opened.
Fits file HN0123456opt00_bs.fits written.
$ imsize HN0123456opt00_bs.fits 
HN0123456opt00_bs.fits 00:00:00.000 +00:00:00.00 J2000 10.035mx11.505m -0.2940/0.2940s/pix  2048x2348 pix

# (cf.) Not compatible mode.
~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -s '_notcompat' HN0123456opt00.fits
Fits file HN0123456opt00.fits opened.
Fits file HN0123456opt00_notcompat.fits written.
$ imsize HN0123456opt00_notcompat.fits 
HN0123456opt00_notcompat.fits 00:00:00.000 +00:00:00.00 J2000 10.035mx11.015m -0.2940/0.2940s/pix  2048x2248 pix

```
# Author

* Hiroshi AKITAYA
* PERC, Chiba Institute of Technology
* akitaya _at_ perc.it-chiba.ac.jp


