"""
Plot chromosoms and misc variants. Can be invoked from commandline
or from imported.

    $ ./chromograph.py --cyt cytoBand.bed

Further info:
  https://matplotlib.org/3.1.1/api/collections_api.html#matplotlib.collections.BrokenBarHCollection


Project on Github:

    https://github.com/mikaell/chromograph

"""
import os
import re
import sys
from argparse import ArgumentParser
import pandas
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from chromograph import __version__
from .chr_utils import (
    read_cfg,
    filter_dataframe,
    outpath,
    parse_wig_declaration,
    parse_upd_regions,
    chr_type_format,
)

matplotlib.use("Agg")

# TODO: instead of padding look-ahead and contsrict if overlap
# TODO: combined ROH image

HELP_STR_RGB = "graph color in RGB hex (only in combination with --coverage)"
HELP_STR_COMBINE = "plot all graphs in one file, default one graph per file"
HELP_STR_NORM = "normalize data (only used for wig/coverage)"
HELP_STR_UPD_SITE = "input UPD sites file on format {}"
HELP_STR_UPD_REGIONS = "input UPD regions file"
HELP_STR_IDEO = "input ideogram (bed-file) on format {}"
HELP_STR_COV = "input fixed step wig-file"
HELP_STR_VSN = "Display program version ({}) and exit."
HELP_STR_EU = "always output an euploid amount of files -some may be empty PNGs"

PADDING = 200000
COVERAGE_END = 249255000  # write all coverage files to the same size canvas
HEIGHT = 1
YBASE = 0
SPACE = 1
FIGSIZE = (6, 8)
FIGSIZE_SINGLE = (8, 8)
UPD_SITES_FORMAT = ["chrom", "start", "end", "updType"]
IDEOGRAM_FORMAT = ["chrom", "start", "end", "name", "gStain"]
WIG_FORMAT = ["chrom", "coverage", "pos"]
WIG_ORANGE = "#e89f00"
WIG_MAX = 70.0
PNG_BYTES = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x01\x00\x00\x00\x007n\xf9$\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01H\xaf\xa4q\x00\x00\x00\x00IEND\xaeB`\x82"
CHROMOSOMES = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "M",
    "X",
    "Y",
]

DEFAULT_SETTING = {"combine": False, "normalize": False, "euploid": False}

get_color = {
    # Cytoband colors
    "gneg": "#f7f7f7",
    "gpos25": "#666666",
    "gpos50": "#606060",
    "gpos75": "#2e2e2e",
    "gpos100": "#121212",
    "acen": "#b56666",
    "gvar": "#777777",
    "stalk": "#444444",
    # UPD sites colors
    "PB_HOMOZYGOUS": "#222222",
    "ANTI_UPD": "#555555",
    "PB_HETEROZYGOUS": "#888888",
    "UNINFORMATIVE": "#333333",
    "UPD_MATERNAL_ORIGIN": "#aa2200",  # Red
    "UPD_PATERNAL_ORIGIN": "#0044ff",  # Blue
    # UPD region colors
    "PATERNAL": "#0044ff",  # Blue
    "MATERNAL": "#aa2200",
    
}


def assure_dir(outd):
    """Create directory 'outd' if it does not exist"""
    print("outd: {}".format(outd))
    if not os.path.exists(outd):
        os.makedirs(outd)


def bed_collections_generator_combine(dataframe, y_positions, height):
    """ Iterate dataframe

        Args:
            dataframe(pandas dataframe)
            y_positions()
            height
        Yields:
            BrokenBarHCollection
    """
    for chrom, group in dataframe.groupby("chrom"):
        print("chrom: {}".format(chrom))
        yrange = (y_positions[chrom], height)
        xranges = group[["start", "width"]].values
        yield BrokenBarHCollection(xranges, yrange, facecolors=group["colors"], label=chrom)


def bed_collections_generator(dataframe, height):
    """ Interate dataframe
    Yeilds:
    Brokenbarhcollection --from BED DF to be plotted
    """
    for chrom, group in dataframe.groupby("chrom"):
        yrange = (0, height)
        xranges = group[["start", "width"]].values
        yield BrokenBarHCollection(xranges, yrange, facecolors=group["colors"], label=chrom)


def coverage_generator(dataframe, data_state):
    """Iterate dataframe and yeild chromosome

    Args:
        df(dataframe)
        data_state('coverage'|'normalized_coverage')

    Yeilds:
        dict -- {'label', 'x', 'y'}
    """
    for chrom, group in dataframe.groupby("chrom"):
        yield {"label": chrom, "x": group["pos"].values, "y": group[data_state].values}


def coverage_generator_combine(dataframe, height):
    """Iterate dataframe and yeild per chromosome, like coverage_generator()
    -with additional positional

    Args:
        dataframe --
        height ---

    Yeilds:
        BrokenBarhcollection
    """
    for chrom, group in dataframe.groupby("chrom"):
        yrange = (0, height)
        xranges = group[["start", "width"]].values
        yield BrokenBarHCollection(xranges, yrange, facecolors=group["colors"], label=chrom)


def common_settings(axis):
    """Common Matplotlib settings to remove visible axis"""
    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)
    axis.set_axis_off()  # Remove black line surrounding pic.


def print_individual_pics(dataframe, infile, outd, euploid):
    """Print one chromosomes per image file"""
    fig = plt.figure(figsize=(10, 0.5))
    axis = fig.add_subplot(111)
    plt.rcParams.update({"figure.max_open_warning": 0})

    is_printed = []
    for collection in bed_collections_generator(dataframe, HEIGHT):
        axis.add_collection(collection)
        common_settings(axis)
        axis.set_xlim((-12159968, 255359341))  # try to mimic nice bounds
        outfile = outpath(outd, infile, collection.get_label())
        print("outfile: {}".format(outfile))
        fig.savefig(outfile, transparent=True, bbox_inches="tight", pad_inches=0)
        axis.cla()  # clear canvas before next iteration
        is_printed.append(collection.get_label())
    if euploid:
        print_empty_pngs(infile, outd, is_printed)


def print_combined_pic(dataframe, chrom_ybase, chrom_centers, infile, outd, chr_list):
    """Print all chromosomes in a single PNG picture"""
    fig = plt.figure(figsize=FIGSIZE)
    axis = fig.add_subplot(111)
    for collection in bed_collections_generator_combine(dataframe, chrom_ybase, HEIGHT):
        axis.add_collection(collection)

    axis.set_yticks([chrom_centers[i] for i in chr_list])
    axis.set_yticklabels(chr_list)
    axis.axis("tight")
    outfile = outpath(outd, infile, "combined")
    print("outfile: {}".format(outfile))
    fig.savefig(outfile, transparent=True, bbox_inches="tight", pad_inches=0)


def print_empty_pngs(file, outd, is_printed):
    """Write an empty png file to disk for every chromosome what has, always including Y.
    Motivated by auxilary software not being able to handle missing output
    chromosome are missing in the wig."""
    gene_build = chr_type_format(is_printed[0])

    for chrom in CHROMOSOMES:
        if chrom in is_printed:
            continue
        if "chr" + chrom in is_printed:
            continue

        prefix = "chr" if gene_build == "str" else ""
        outfile = outpath(outd, file, prefix + chrom)
        print("print empty: {}".format(outfile))
        filestream = open(outfile, "bw")
        filestream.write(PNG_BYTES)
        filestream.close()


def bed_to_dataframe(file, spec):
    """Read a bed file into a Pandas dataframe according to 'spec' """
    return pandas.read_csv(file, names=spec, sep="\t", skiprows=1)


def wig_to_dataframe(infile, step, col_format):
    """Read a wig file into a Pandas dataframe

    infile(str): Path to file

    Returns:
        Dataframe

    """
    filestream = open(infile, "r")
    coverage_data = []
    pos = 0
    chrom = ""
    for line in filestream.readlines():
        try:
            wig_value = float(line)
            wig_value = wig_value if wig_value < WIG_MAX else WIG_MAX
            coverage_data.append([chrom, wig_value, pos])
            pos += step
        except ValueError:
            reresult = re.search("chrom=(\w*)", line)  # find chromosome name in line
            if reresult:
                print(line)
                start_pos = [chrom, 0, 1]  # write 0 in beginning, works against bug
                stop_pos = [chrom, 0, pos + 1]  # write 0 at end removes linear slope
                last_pos = [chrom, 0, COVERAGE_END]  # write trailing char for scale when plotting
                #               coverage_data.insert(start_pos)
                coverage_data.append(stop_pos)
                coverage_data.append(last_pos)
                chrom = reresult.group(1)  # start working on next chromosome
                pos = 0
    filestream.close()
    dataframe = pandas.DataFrame(coverage_data, columns=col_format)
    return dataframe


def graph_coordinates(list_of_chromosomes):
    """Iterate through list of chromosomes and return X
    (as center for graph) and Y coordinates for plotting.

        Args:

        Returns:
            chrom_ybase(int)
            chrom_centers(int)
    """
    ybase = YBASE
    chrom_ybase = {}
    chrom_centers = {}
    for chrom in list_of_chromosomes:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + HEIGHT / 2.0
        ybase += HEIGHT + SPACE
    return chrom_ybase, chrom_centers


def plot_ideogram(file, *args, **kwargs):
    """Visualize chromosome ideograms from bed-file. Format:

    Args:
        file(string path)

    Optional Args:
        combine -- output all graphs in one png
        outd=<str> -- output directory

    Returns:
          None
    """
    cfg = read_cfg()
    settings = normalize_args(file, args, kwargs)
    print(
        "Plot ideograms with settings\ncombine:{}\noutd:{}".format(
            settings["combine"], settings["outd"]
        )
    )

    # try two different chromosome formats, if both result in
    # an empty dataframe, raise IdeoParseError
    for setting in ["chromosome_str", "chromosome_int"]:
        chromosome_list = cfg[setting]
        dataframe = bed_to_dataframe(file, IDEOGRAM_FORMAT)
        # delete chromosomes not in CHROMOSOME_LIST
        dataframe = filter_dataframe(dataframe, chromosome_list)
        if dataframe.size > 0:
            break
    if dataframe.size == 0:
        raise Exception("Ideogram parsing")
    dataframe["width"] = dataframe.end - dataframe.start
    dataframe["colors"] = dataframe["gStain"].apply(lambda x: get_color[x])
    chrom_ybase, chrom_centers = graph_coordinates(chromosome_list)
    if settings["combine"]:
        print_combined_pic(
            dataframe, chrom_ybase, chrom_centers, file, settings["outd"], chromosome_list
        )
    else:
        print_individual_pics(dataframe, file, settings["outd"], settings["euploid"])


def plot_roh(file, *args, **kwargs):
    """Plot ROH file for analysis of isodisomy or heterodisomy"""
    return "TODO"


def plot_upd_sites(filepath, *args, **kwargs):
    """Visualize UPD data from bed-file. Bed format as:

        Chromosome <tab> Start <tab> End <tab> Upd-type

    These can be generated with SV caller upd.py. It can be found at
    [https://github.com/bjhall/upd]

    Args:
        filepath<str> -- input file on wig format

    Optional Args:
        combine -- output all graphs in one png
        normalize -- normalize to mean
        outd=<str> -- output directory
    Returns: None

    """
    cfg = read_cfg()
    settings = normalize_upd_sites_args(filepath, args, kwargs)

    print(
        "Plot UPD Sites with settings\ncombine:{}\neuploid:{}".format(
            settings["combine"], settings["euploid"]
        )
    )
    dataframe = bed_to_dataframe(filepath, UPD_FORMAT)
    dataframe.chrom = dataframe.chrom.astype(str)  # Explicitly set chrom to string (read as int)
    chromosome_list = cfg["chromosome_int"]
    dataframe = filter_dataframe(
        dataframe, chromosome_list
    )  # delete chromosomes not in CHROMOSOME_LIST_UPD
    dataframe["width"] = (dataframe.end - dataframe.start) + PADDING
    dataframe["colors"] = dataframe["updType"].apply(lambda x: get_color[x])
    chrom_ybase, chrom_centers = graph_coordinates(chromosome_list)
    if settings["combine"]:
        print_combined_pic(
            dataframe, chrom_ybase, chrom_centers, filepath, settings["outd"], chromosome_list
        )
    else:
        print_individual_pics(dataframe, filepath, settings["outd"], settings["euploid"])


def plot_wig(filepath, *args, **kwargs):
    """Outputs png:s of data given on WIG format

    Args:
        filepath(str)

    Optional Args:
        combine -- output all graphs in one png
        normalize -- normalize to mean
        outd(str) -- output directory
    Returns: None

    """
    cfg = read_cfg()
    decl = parse_wig_declaration(filepath, " ")
    settings = normalize_wig_args(decl, filepath, args, kwargs)

    print(
        "Plot WIG with settings \nstep: {}\noutd:{}\ncombine:{}\nnormalize:{}\neuploid:{}".format(
            settings["fixedStep"],
            settings["outd"],
            settings["combine"],
            settings["normalize"],
            settings["euploid"],
        )
    )

    chromosome_list = list_type(decl["chrom"], cfg)
    dataframe = wig_to_dataframe(filepath, settings["fixedStep"], WIG_FORMAT)
    dataframe = filter_dataframe(
        dataframe, chromosome_list
    )  # delete chromosomes not in CHROMOSOME_LIST
    dataframe["normalized_coverage"] = (dataframe.coverage / dataframe.coverage.mean()).round(0)
    print_wig(
        dataframe,
        filepath,
        settings["outd"],
        settings["combine"],
        settings["normalize"],
        settings["color"],
        settings["euploid"],
    )


def print_wig(dataframe, file, outd, combine, normalize, color, euploid):
    """Print wig graph as PNG file"""
    data_state = "normalized_coverage" if normalize else "coverage"
    if not combine:  # Plot one chromosome per png
        is_printed = []
        for chrom_data in coverage_generator(dataframe, data_state):
            fig, axis = plt.subplots(figsize=(8, 0.7))
            common_settings(axis)
            axis.stackplot(chrom_data["x"], chrom_data["y"], colors=color)
            plt.ylim(0, 5)
            axis.set_ylim(bottom=0)
            axis.set_xlim((0, 255359341))  # try to mimic nice bounds
            fig.tight_layout()
            outfile = outpath(outd, file, chrom_data["label"])
            print("outfile: {}".format(outfile))
            fig.savefig(outfile, transparent=True, bbox_inches="tight", pad_inches=0)
            is_printed.append(chrom_data["label"])
            plt.close(fig)  # save memory
        if euploid:
            print_empty_pngs(file, outd, is_printed)
    else:
        # TODO:
        print("WARNING: Combined WIG png not implemented!")
        False


def plot_upd_regions(file, *args, **kwargs):
    """  Print region as PNG file"""
    # These can show isodisomy or hetereodisomy in description

    # Parse sites upd file to brokenbarcollection
    read_line = []
    settings = normalize_args(file, args, kwargs)

    print(
        "Plot UPD REGIONS with settings \noutd:{}\neuploid: {}".format(
            settings["outd"], settings["euploid"]
        )
    )
    with open(file) as filepointer:
        for line in filepointer:
            read_line.append(parse_upd_regions(line))
    # TODO: 'sites' in a function called 'plot_regions'
    broken_bar = [sites_to_brokenbar(i) for i in read_line]

    # Prepare canvas and plot
    fig = plt.figure(figsize=(10, 0.5))
    x_axis = fig.add_subplot(111)
    is_printed = []
    for collection in broken_bar:
        x_axis.add_collection(collection)
        common_settings(x_axis)
        x_axis.set_xlim((-12159968, 255359341))  # try to mimic nice bounds
        outfile = outpath(settings["outd"], file, collection.get_label())
        fig.savefig(outfile, transparent=False, bbox_inches="tight", pad_inches=0)
        x_axis.cla()  # clear canvas before next iteration
        is_printed.append(collection.get_label())

    # print each name only once
    for name in dict.fromkeys(is_printed):
        print("outfile: {}".format(outpath(settings["outd"], file, name)))
    if settings["euploid"]:
        print_empty_pngs(file, settings["outd"], is_printed)


def sites_to_brokenbar(region):
    """ Make a MathPlotLIb 'BrokenbarCollection' from upd sites data"""
    start = int(region["start"])
    width = int(region["stop"]) - int(region["start"])
    color = get_color[region["desc"]["origin"]]  
    color = get_color[region["desc"]["TYPE"]]    # <--- handle iso/hetereo here
    xranges = [[start, width]]
    yrange = (0, 1)
    return BrokenBarHCollection(xranges, yrange, facecolors=color, label=region["chr"])


def normalize_wig_args(decl, filepath, args, kwargs):
    """Override default settings if argument is given, return settings dict
    for coverage/wig"""
    settings = normalize_args(filepath, args, kwargs)
    settings["color"] = WIG_ORANGE  # set default color, override if rgb in kw
    settings["fixedStep"] = decl["step"]  # cfg['wig_step']
    if "rgb" in kwargs and kwargs["rgb"] is not None:
        settings["color"] = rgb_str(kwargs["rgb"])
    if "step" in kwargs and kwargs["step"] is not None:
        settings["fixedStep"] = kwargs["step"]
    return settings


def normalize_upd_sites_args(filepath, args, kwargs):
    """Override default settings if argument is given, return settings dict
    for upd sites"""
    settings = normalize_args(filepath, args, kwargs)
    settings["outd"] = os.path.dirname(filepath)
    if "combine" in args:
        settings["combine"] = True
    return settings


def normalize_args(filepath, args, kwargs):
    """Return a dict of args set to default if not given"""
    settings = DEFAULT_SETTING
    settings["outd"] = os.path.dirname(filepath)

    if "combine" in args:
        settings["combine"] = True
    if "norm" in args:
        settings["normalize"] = True
    if "euploid" in args:
        settings["euploid"] = True
    if "outd" in kwargs and kwargs["outd"] is not None:
        print("output directory:{}".format(kwargs["outd"]))
        assure_dir(kwargs["outd"])
        settings["outd"] = kwargs["outd"]
    return settings


def rgb_str(color):
    """Return #color"""
    head, *_tail = str(color)
    if head == "#":
        return color  # color was alread on format "#123456"
    return "#" + str(color)


def list_type(chr_type, cfg):
    """ Determine type and return list of possible
    field names in parsed data frame"""
    if chr_type == "int":
        return cfg["chromosome_int"]
    return cfg["chromosome_str"]


def main():
    """Main function for Chromograph

    Parse incoming args and call correct function"""
    parser = ArgumentParser()
    parser.add_argument(
        "-u", "--sites", dest="updfile", help=HELP_STR_UPD_SITE.format(UPD_SITES_FORMAT), metavar="FILE"
    )
    parser.add_argument(
        "-g", "--regions", dest="regionsfile", help=HELP_STR_UPD_REGIONS, metavar="FILE"
    )
    parser.add_argument("-e", "--ideo", dest="ideofile", help=HELP_STR_IDEO.format(IDEOGRAM_FORMAT), metavar="FILE")
    parser.add_argument("-w", "--coverage", dest="coverage_file", help=HELP_STR_COV, metavar="FILE")
    parser.add_argument("-y", "--roh", dest="roh", help="TODO: regions of homozygosity", metavar="FILE")
    parser.add_argument("-o", "--outd", dest="outd", help="output dir", metavar="FILE")
    parser.add_argument("-r", "--rgb", dest="rgb", help=HELP_STR_RGB, metavar="FILE")
    parser.add_argument("-n", "--norm", dest="norm", help=HELP_STR_NORM, action="store_true")
    parser.add_argument("--step", type=int, help="fixed step size (default 5000)")
    parser.add_argument("-c", "--combine", help=HELP_STR_COMBINE, action="store_true")
    parser.add_argument("--version", help=HELP_STR_VSN.format(__version__),action="version", version="chromograph {}".format(__version__) )
    parser.add_argument("-p", "--euploid", help= HELP_STR_EU, action="store_true")

    args = parser.parse_args()

    # Make command line and library interfaces behave identical regarding args
    args.norm = "norm" if args.norm else None
    args.combine = "combine" if args.combine else None
    args.euploid = "euploid" if args.euploid else None

    if args.ideofile:
        plot_ideogram(args.ideofile, args.combine, outd=args.outd)
    if args.updfile:
        plot_upd_sites(args.updfile, args.combine, args.euploid, outd=args.outd, step=args.step)
    if args.roh:
        plot_roh(args.roh, args.combine, args.euploid, outd=args.outd)
    if args.coverage_file:
        plot_wig(
            args.coverage_file,
            args.combine,
            args.norm,
            args.euploid,
            outd=args.outd,
            step=args.step,
            rgb=args.rgb,
        )
    if args.regionsfile:
        plot_upd_regions(args.regionsfile, args.euploid, outd=args.outd)
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()


if __name__ == "__main__":
    main()
