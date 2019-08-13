# Chromograph

## What is Chromograph?
Chromograph is a python package to create PNG images from genetics data
such as BED and WIG files. It is primarily build to integrate with
software form [Clinical Genomics](https://github.com/Clinical-Genomics) and
[UPD Tool](https://github.com/bjhall/upd) but is stand-alone.

Integrated with [Scout](https://github.com/Clinical-Genomics/scout)
Chromograph can be used to visualise chromosomes:

![screenshot](https://github.com/mikaell/chromograph/blob/master/scout_example.png)


## Usage, command line
Chromograph as used from the command line.

```
Usage: chromograph [COMMAND] [OPTIONS]

Command (one must be provided):
  --upd <filepath>
  --ideo <filepath>
  --roh <filepath>

Options:
  --help        Show help message
  --combine     Default one image/png
  --normalize   normalize to mean (roh file only)
  --step <int>  Default 5000 (roh file only)
  --outd <path> Default outdir is infile
```

### Example
```
$ ./chromograph.py -r roh.wig --combine --step 3000
```

## Usage, lib
Chromograph used as module. File must be provided, other arguments are
optional. Example:
```
>>> import chromograph as chrm
>>> chrm.plot_ideogram("path_file.bed", 'combine', outdir=<path_output>)
>>> chrm.plot_upd(<file>, outd=<path>)
>>> chrm.plot_roh(<file>, 'combine', 'normalize', outd=<path>, step=<int>)
>>> chrm.plot_ideogram(file, 'combine', outd=<path>)
>>> chr.plot_ideogram("~/Work/files_test/cytoband.bed", 'combine', outd="/Users/Mikael/Downloads")

```

## Configuration
Configuration file is `config_chromograph.yml`.

* `chromosome_str` used by roh and ideogram.
* `chromosome_int` used by upd.
* `wig_step` fixed size step in wig-files.




## Data Formats
### UPD BED
UPD bed files are supported on format:
```
'chrom', 'start', 'end', 'updType'
```

### Ideogram BED
Ideogram bed files are supported on format:
```
IDEOGRAM_FORMAT = ['chrom', 'start', 'end', 'name', 'gStain']
```

### ROH Coverage WIG
ROH wig files are supported on format:
```
['chrom', 'coverage', 'pos']
```

### Requirements
 * MatPlotLib
 * Numpy
 * Pandas

## Installation
```
pip install numpy; pip install matplotlib; pip install pandas
python setup.py install
```
