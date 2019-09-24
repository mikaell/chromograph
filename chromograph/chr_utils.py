import os
import yaml
import pkg_resources

def read_cfg():
    """Read Yaml config"""
    cfg_path = pkg_resources.resource_filename('chromograph', 'config_chromograph.yml')
    try:
        with open(cfg_path, 'r') as ymlfile:
            return yaml.safe_load(ymlfile)
    except FileNotFoundError:
        return \
       {'chromosome_int': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                           '11', '12', '13', '14', '15', '16', '17', '18',
                           '19', '20','21', '22', 'M', 'X', 'Y'],
        'chromosome_str': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                           'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                           'chr14' , 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                           'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY'],
        'wig_step': 5000}



def filter_dataframe(frame, list_of_chromosomes):
    """Delete chromosome names not read into 'frame'

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
    f = os.path.basename(infile)    # strip path from filename
    outfile = png_filename(f, label)
    outpath = outd + "/" + outfile
    print("outfile: {}".format(outpath))
    return outpath

