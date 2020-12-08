# Chromograph


---->
## THIS README IS NOT UP TO DATE WITH LATEST CHANGES AND NEEDS UPDATE
<----



## What is Chromograph?
Chromograph is a python package to create PNG images from genetics data
such as BED and WIG files. It is primarily build to integrate with
software form [Clinical Genomics](https://github.com/Clinical-Genomics) and
[UPD Tool](https://github.com/bjhall/upd) but is stand-alone.

UPD Tool is 
Integrated with [Scout](https://github.com/Clinical-Genomics/scout)
Chromograph can be used to visualise chromosomes:

![screenshot](https://github.com/mikaell/chromograph/blob/master/scout_example.png)


## Usage, command line
Chromograph as used from the command line.

```
usage: chromograph [-h] [-a FILE] [-c FILE] [-f FILE] [-i FILE] [-r FILE]
                   [-s FILE] [--step STEP] [--version] [-d FILE] [-e]
                   [-k FILE] [-n] [-x]

optional arguments:
  -h, --help            show this help message and exit
  -a FILE, --autozyg FILE
                        Plot regions of autozygosity from bed file [OPERATION]
  -c FILE, --coverage FILE
                        Plot coverage from fixed step wig file [OPERATION]
  -f FILE, --fracsnp FILE
                        Plot fraction of homozygous SNPs from wig file
                        [OPERATION]
  -i FILE, --ideogram FILE
                        Plot ideograms from bed-file on format ['chrom',
                        'start', 'end', 'name', 'gStain'] [OPERATION]
  -r FILE, --regions FILE
                        Plot UPD regions from bed file [OPERATION]
  -s FILE, --sites FILE
                        Plot UPD sites from bed file [OPERATION]
  --step STEP           fixed step size (default 5000)
  --version             Display program version (1.0.0) and exit.
  -d FILE, --outd FILE  output dir
  -e, --euploid         Always output an euploid amount of files -even if some
                        are empty
  -k FILE, --rgb FILE   Set color (RGB hex, only with --coverage option)
  -n, --norm            Normalize data (wig/coverage)
  -x, --combine         Write all graphs to one file (default one plot per
                        file)

One OPERATION Command is needed for Chromograph to produce output
```

### Example
```
$ ./chromograph.py --autozyg rhocall.bed --outd tmp/
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



## Data Formats
### UPD WIG, Regions



### UPD BED, Sites
UPD bed files are supported on format:
```
'chrom', 'start', 'end', 'updType'
```
UPD bed files are generated with the tool UPDtool. An example call: `TODO`


### Ideogram BED
Ideogram bed files are supported on format:
```
IDEOGRAM_FORMAT = ['chrom', 'start', 'end', 'name', 'gStain']
```

### Coverage WIG
ROH wig files are supported on format:
```
['chrom', 'coverage', 'pos']
```
These can be created using Tiddit. For example: `TODO`


### Regions of Autozygosity, BED


### Requirements
Chromograph runs in Python 3.

Package requirements:
 * MatPlotLib
 * Numpy
 * Pandas
 * Pyaml

## Installation
```
pip install -r requirements.txt -e .
```



## Deprecated Functionality 
* Configuration file removed in version 0.3.2.