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


COVERAGE_WIG_PNG = "a736401fe0ccba6bf616ed92a4679da1"


def test_coverage():
    # GIVEN example coverage file
    # WHEN calling chromograph
    chrom.plot_coverage_wig(coverage_example)

    # Python program to find MD5 hash value of a file
 
    with open("tests/example_files/coverage_chr1.png", "rb") as f:
        bytes = f.read() # read file as bytes
        readable_hash = hashlib.md5(bytes).hexdigest();
        assert readable_hash == COVERAGE_WIG_PNG
    

    # THEN produced files checksom will be correct
    

    
