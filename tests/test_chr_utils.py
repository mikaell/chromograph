"""Pytests for Chromograph """
import os
import unittest.mock as mock
from unittest.mock import mock_open
import pandas as pd
from chromograph.chr_utils import (chr_type_format, cast, filter_dataframe,
                                   png_filename, outpath, parse_wig_declaration,
                                   make_dict)



WIG_HEAD="""123
312312
12321

fixedStep chrom=1 start=1 step=10000"""


def test_chr_type_format():
    # GIVEN a string integer, i.e. '1'
    # THEN chrFormat returns 'int'
    assert chr_type_format('1') == 'int'


def test_chr_type_format_str():
    # GIVEN a string, not representing an integer, i.e. 'chr3'
    # THEN chrFormat returns 'str'
    assert chr_type_format('chr4') == 'str'


def test_make_dict():
    # GIVEN a list where each element is on format: '<char>=<int>'
    test_list = ['a=1', 'b=2', 'c=3']
    # THEN the list is transformed to a key/value dict, <char>:<int>
    test_dict = {'a':'1', 'b':'2', 'c':'3'}
    #
    assert test_dict == make_dict(test_list)


def test_png_filename():
    # GIVEN filename "test.file" and label "myLabel"
    # THEN png_filename() will return "test_myLabel.png"
    assert png_filename("test.file", "myLabel") == "test_myLabel.png"


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
    wig_header = {'chrom': '1', 'start': '1', 'step': '10000\n'}
    # THEN calling `cast` parsed values are typecast to internally used types
    assert cast(wig_header) == {'chrom': 'int', 'start': 1, 'step': 10000}


def test_filter_dataframe():
    # GIVEN a small dataframe
    test_frame = pd.DataFrame({'chrom':['chr1', 'chrERR'], 'coverage':[1, 2], 'pos':[1, 2]})
    # THEN entries not matching the filter are removed
    clean_frame = pd.DataFrame({'chrom':['chr1'], 'coverage':[1], 'pos':[1]})
    filtered_frame = filter_dataframe(test_frame, ['chr1'])
    # CAST frames to dicts to avoid disambigous series
    assert filtered_frame.to_dict() == clean_frame.to_dict()


@mock.patch('builtins.open', new_callable=mock_open, read_data=WIG_HEAD)
def test_parse_wig_declaration(mock_file):
    # GIVEN a mockup wig file with fixed steps, nonsense lines and fixedStep line

    # THEN chrom, start and step are parsed to a dict
    assert {'chrom': 'int', 'start': 1, 'step': 10000} == parse_wig_declaration('filename', ' ')


# TODO: add test to catch warning('declarationNotFound') when declaration
# is missing, test_parse_wig_declaration_warn()
