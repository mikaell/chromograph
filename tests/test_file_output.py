import chromograph.chromograph as chrom
from chromograph import (
    autozygosity_example,
    coverage_example,
    cytoband_example,
    homozygous_snp_example,
    upd_regions_example,
    upd_sites_example,
)
import hashlib

"""Test suite for Chromograph

   Note that different files work on different chromosomes.
"""
# TODO: use tmp-dir for files written


IDEOGRAM_MD5 = "83ba68e788730ada6ce77c7e2c0ebde0"
CHROM1_MD5 = "f5b32773a52d0975e6935072bd02e5d5"
COVERAGE_MD5 = "570acd69ed69d868e2f7b9f638a2e480"
UPD_REGIONS_MD5 = "001e43b6e85ebe6edfea59dd07139862"
UPD_SITES_MD5 = "620c70adcf5197f07647f2ec7e69ca6c"

def test_ideogram():
    result = chrom.plot_ideogram(cytoband_example)
    print("RESULT")
    print(result)
    with open("tests/example_files/cytoband_chr1.png", "rb") as f:
        assert f.read() # read file as bytes
        f.close()

def test_coverage():
    chrom.plot_coverage_wig(coverage_example)
    with open("tests/example_files/coverage_chr1.png", "rb") as f:
        assert f.read() # read file as bytes
        f.close()

def test_upd_regions():
    chrom.plot_upd_regions(upd_regions_example)
    with open("tests/example_files/upd_regions_12.png", "rb") as f:
        assert f.read() # read file as bytes
        f.close()

def test_upd_sites():
    chrom.plot_upd_sites(upd_sites_example)
    with open("tests/example_files/upd_sites_1.png", "rb") as f:
        assert f.read() # read file as bytes
        f.close()

    


    

    
