import os
import yaml

def read_cfg():
    """Read Yaml config"""
    with open("config_chromograph.yml", 'r') as ymlfile:
        return yaml.safe_load(ymlfile)


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

