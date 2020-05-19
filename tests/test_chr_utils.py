#
#
#
import os, pathlib

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
    s = ['a=1','b=2','c=3']
    # THEN the list is transformed to a key/value dict, <char>:<int>
    d = {'a':'1', 'b':'2', 'c':'3'}
    # 
    assert d == makeDict(s)



def test_png_filename():
    # GIVEN filename "test.file" and label "myLabel"
    # THEN png_filename() will return "test_myLabel.png"
    assert "test_myLabel.png" == png_filename("test.file", "myLabel")

    
    # TISDAG: skriv klart kommentarer
def test_outpath(tmpdir):
    # WHEN setting up temporary directory and a file
    outd = tmpdir.mkdir("directory")
    infile = "file.test"
    # GIVEN
    d = os.path.join(outd, "file_LABEL.png")
    # THEN outpath
    assert d == outpath(outd, infile, "LABEL")
    



# /private/var/folders/vc/5c_gq8jj3r5fgmwq7r5k8xk00000gp/T/pytest-of-Mikael/pytest-13/test_outpath0/testOPinfile_tester.png


    
# def test_read_ini(tmpdir):
#     print(tmpdir)      # /private/var/folders/ry/z60xxmw0000gn/T/pytest-of-gabor/pytest-14/test_read0
#     d = tmpdir.mkdir("subdir")
#     fh = d.join("config.ini")
#     fh.write("""
# [application]
# user  =  foo
# password = secret  
# """)
 
#     print(fh.basename) # data.txt
#     print(fh.dirname)  # /private/var/folders/ry/z60xxmw0000gn/T/pytest-of-gabor/pytest-14/test_read0/subdir
#     filename = os.path.join( fh.dirname, fh.basename )
