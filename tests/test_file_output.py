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

IDEOGRAM_MD5 = "84095d4ab8a64756c647e2248373a2da"
CHROM1_MD5 = "f5b32773a52d0975e6935072bd02e5d5"
COVERAGE_MD5 = "a736401fe0ccba6bf616ed92a4679da1"
UPD_REGIONS_MD5 = "01c8692e23b75af845c73d1e81afc30a"
UPD_SITES_MD5 = "f29a92edffd43c4d2a3879f7559233ef"

def test_ideogram():
    chrom.plot_ideogram(cytoband_example)
    with open("tests/example_files/cytoband_chr1.png", "rb") as f:
        bytes = f.read() # read file as bytes
        readable_hash = hashlib.md5(bytes).hexdigest();
        # THEN produced files checksom will be correct
        assert readable_hash == IDEOGRAM_MD5

def test_coverage():
    chrom.plot_coverage_wig(coverage_example)
    with open("tests/example_files/coverage_chr1.png", "rb") as f:
        bytes = f.read() # read file as bytes
        readable_hash = hashlib.md5(bytes).hexdigest();
        # THEN produced files checksom will be correct
        assert readable_hash == COVERAGE_MD5


def test_upd_regions():
    chrom.plot_upd_regions(upd_regions_example)
    with open("tests/example_files/upd_regions_12.png", "rb") as f:
        bytes = f.read() # read file as bytes
        readable_hash = hashlib.md5(bytes).hexdigest();
        assert readable_hash == UPD_REGIONS_MD5


def test_upd_sites():
    chrom.plot_upd_sites(upd_sites_example)
    with open("tests/example_files/upd_sites_1.png", "rb") as f:
        bytes = f.read() # read file as bytes
        readable_hash = hashlib.md5(bytes).hexdigest();
        assert readable_hash == UPD_SITES_MD5


    


    

    
