from pkg_resources import get_distribution

__version__ = get_distribution('chromograph').version

# Test files
autozygosity_example = "tests/example_files/autozygosity.bed"
coverage_example = "tests/example_files/coverage.wig"
cytoband_example = "tests/example_files/cytoband.bed"
homozygous_snp_example = "tests/example_files/fractionHomozygousSNP.wig"
upd_regions_example = "tests/example_files/upd_regions.bed"
upd_sites_example = "tests/example_files/upd_sites.bed"




