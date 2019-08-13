# Chromograph

## What is Chromograph?
Chromograph is a python package to create PNG images from genetics data
such as BED and WIG files. It is primarily build to integrate with
software form [Clinical Genomics](https://github.com/Clinical-Genomics) and
[UPD Tool](https://github.com/bjhall/upd) but is stand-alone.

Integrated with [Scout] (https://github.com/Clinical-Genomics/scout):


![screenshot](https://github.com/mikaell/chromograph/blob/master/scout_example.png)




## Usage, command line
```
Usage: chromograph [COMMANT] [OPTIONS]

Command:
  --upd <filepath>
  --ideo <filepath>
  --roh <filepath>

Options:
  --help        Show help message
  --step <int>  Default 5000
  --combine     Default one image/png
```



## Usage, lib
```
>>> import chromograph as chrm
>>> chrm.plot_ideogram("path_file.bed", 'combine', outdir=<path_output>)



plot_upd(<file>, (combine), (outd=<path>))


plot_roh(<file>, 'combine', 'normalize', outd=<path>, step=<int>)

plot_ideogram(file, 'combine', outd=<path>)
chr.plot_ideogram("~/Work/files_test/cytoband.bed", 'combine', outd="/Users/Mikael/Downloads")

```




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
['chrom', 'coverage', 'pos']

### Requirements
 * MatPlotLib
 * Numpy
 * Pandas

## Example
Genereate ...
```
$ ./chromograph.py -r roh.wig --step 3000 --combine
```
## Installation
```
pip install numpy; pip install matplotlib; pip install pandas
python setup.py install
```




