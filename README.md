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
$ pip install --user -e .
$ # or pip install --user .
```
# Usage

## ccd/hpkossub.py as a command tool.
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
   * example 1
   ```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py HP0277904_0.fits
Fits file HP0277904_0.fits opened.
Fits file HP0277904_0.bs.fits written.
   ```
   * example 2
   ```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -o HP0277904_0.fits
Fits file HP0277904_0.fits opened.
File HP0277904_0.bs.fits exists. Overwrite.
Fits file HP0277904_0.bs.fits written.
   ```
   * example 3
   ```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py HP027790[12]_0.fits
Fits file HP0277901_0.fits opened.
Fits file HP0277901_0.bs.fits written.
Fits file HP0277902_0.fits opened.
Fits file HP0277902_0.bs.fits written.
   ```
   * example 4
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
   * example 5
   ```
$ ~/iraf/kanatapy/kanatapy/ccd/hpkossub.py -s "_ossub" HP0277958_0.fits
Fits file HP0277958_0.fits opened.
Fits file HP0277958_0_ossub.fits written.
   ```
## ccd/hpkossub.py as a module.
```
In [1]: from kanatapy.ccd import hpkossub
In [2]: hs = hpkossub.HPKOsSub('HP0277959_0.fits', sub_extention='_ossub')
In [3]: hs.ossub_all()
Fits file HP0277959_0.fits opened.
Fits file HP0277959_0_ossub.fits written.
```


# Author

* Hiroshi AKITAYA
* PERC, Chiba Institute of Technology
* akitaya _at_ perc.it-chiba.ac.jp


