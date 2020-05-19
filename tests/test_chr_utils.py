#
#
#
import os

from chromograph.chr_utils import (chrFormat, cast, read_cfg, filter_dataframe,
                                   png_filename, outpath, parseWigDeclarationLine,
                                   makeDict)

def test_chrFormat():
    # GIVEN a string integer, i.e. '1'
    # THEN chrFormat returns 'int'
    assert chrFormat('1') == 'int'


def test_chrFormat_str():
    # GIVEN a string, not representing an integer, i.e. 'chr3'
    # THEN chrFormat returns 'str'
    assert chrFormat('chr4') == 'str'


def test_makeDict():
    # GIVEN a list where each element is on format: '<char>=<int>'
    test_list = ['a=1', 'b=2', 'c=3']
    # THEN the list is transformed to a key/value dict, <char>:<int>
    test_dict = {'a':'1', 'b':'2', 'c':'3'}
    #
    assert test_dict == makeDict(test_list)


def test_png_filename():
    # GIVEN filename "test.file" and label "myLabel"
    # THEN png_filename() will return "test_myLabel.png"
    assert "test_myLabel.png" == png_filename("test.file", "myLabel")


def test_outpath(tmpdir):
    # WHEN setting up temporary directory and a file
    outd = tmpdir.mkdir("directory")
    infile = "file.test"
    # GIVEN simulated result, on format `path/to/directory/file_LABEL.png`
    test_dirpath = os.path.join(outd, "file_LABEL.png")
    # THEN outpath will match
    assert test_dirpath == outpath(outd, infile, "LABEL")


def test_cast():
    # GIVEN a dict, simulated to parsing a Wig-file header
    wigHeader = {'chrom': '1', 'start': '1', 'step': '10000\n'}
    # THEN calling `cast` parsed values are typecast to internally used types
    assert cast(wigHeader) == {'chrom': 'int', 'start': 1, 'step': 10000}
