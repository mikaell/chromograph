
"""CHR_UTILS

A collection of auxillary functions for Chromograph
"""

import os
import pkg_resources


def filter_dataframe(frame, list_of_chromosomes):
    """Delete dataframe entries where 'chrome' does not appear in the
    list `list_of_chromosomes`

        Args:
            frame(dataframe)
            list_of_chromosomes
    """
    return frame[frame.chrom.apply(lambda x: x in list_of_chromosomes)]


def png_filename(infile, label):
    """Return filename with 'label' and suffix 'png'"""
    (filename, _ending) = os.path.splitext(infile)
    return filename + "_" + label + ".png"


def outpath(outd, infile, label):
    """Create a filepath for Chromograph's export files with label"""
    stripped_file = os.path.basename(infile)
    file_path = png_filename(stripped_file, label)
    return os.path.join(outd, file_path)


def parse_wig_declaration(wigfile, separator=" "):
    """Parse a wig file's declaration line on format:

           fixedStep  chrom=chrN start=position  step=stepInterval

        Only supports fixed step
    """
    max_iterations = 12         # parse no more lines after
    i = 0
    fp = open(wigfile)
    while i < max_iterations:
        line = fp.readline()
        x, *xs = line.split(separator)
        if x.lower() == 'fixedstep' and 'chrM' not in xs[0]:
            declaration = make_dict(xs)   # split xs on '=' to get a dict
            return cast(declaration)
        i += 1
    raise Warning('declarationNotFound')


def parse_bed(filepath, separator="\t"):
    """Return 'int' or 'str' if chromosomes are on format 'chr12'. Lines
    starting with '#' are treated as comments. """
    max_iterations = 12         # parse no more lines after
    i = 0
    fp = open(filepath)
    while i < max_iterations:
        line = fp.readline()
        if line.startswith("#"):
            continue
        x, *xs = line.split(separator)
        print("line: {}".format(x))
        return chr_type_format(x)

    raise Warning('declarationNotFound')


def cast(decl):
    """ Cast dict's values from string to internally used types """
    decl['start'] = int(decl['start'])
    decl['step'] = int(decl['step'])
    decl['chrom'] = chr_type_format(decl['chrom'])
    return decl


def chr_type_format(chrom):
    """If chrom is an integer return 'int' else return 'str'
    """
    try:
        int(chrom)
        return 'int'
    except ValueError:
        return 'str'


def parse_upd_regions(line):
    """Parse sites upt file bed file"""
    [chrom, start, stop, desc] = line.split()
    # TODO: fix bug when line ends with ws
    # For example getting: ['chrom=17', 'start=1', 'step=5000', '\n']
    #
    return {'chr': chrom, 'start': start, 'stop': stop, 'desc': make_dict(desc.split(';'))}


def make_dict(kv_list):
    """Iterate a list and split every element on '=' returning
       a dictionary.
          * returned keys in low case
    """
    d = {}
    for i in kv_list:
        k, v = i.split("=")
        d[k.lower()] = v
    return d
