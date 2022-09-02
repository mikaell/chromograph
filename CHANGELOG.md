# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [1.1.5]
### [Fixed]
- Write transparent png:s instead of empty (which actually where not empty)

## [1.1.4]
### [Added]
- Add option `--chunk / -u` to set Matplotlib chunk size from commandline

## [1.1.3]
### [Changed]
- Exit(0) instead of exit(1) when empty bedfile given as input

## [1.1.2]
### [Fixed]
- Typo bug

## [1.1.1]
### [Fixed]
-Adjust chunksize to stop Matplotlib from throwing OverflowError

## [1.1]
### Added
-Test files of every supported type
## Changed
-Isodisomy to replace homodisomy in code, backwards compatible
### Fixed
-Bug fixed when parsing UPD regions
-Delete tmp files added to repo by mistake

## [1.0.2]
### Fixed
-Better check for 37 vs 38 name convention
-Handle empty string in parsing
-Check for empty dataframe

## [1.0.1]
### Fixed
-Don't crash on empty BED file

## [1.0.0]
### Added
-Add support for autozygosity bed:s.
-Add support for fraction of homozygosity wig:s
### Changed
-Command line arguments name and flags are changed. Not backwards compatible.
-Printed colors daker and uniform
### Fixed
-Bug that wrote over output
-Debug prints removed

## [0.3.6]
### Added
- Region UPD pictures are now split in the middle to visually indicate two
chromosomes.

## [0.3.5]
### Added
- Pixel perfect prints 7750 X 385 px

### [0.3.4]
### Added
- Picture DPI increase from 100 X 100 to 1000 X 1000
### Fixed
- Alignment of BED files was offset

## [0.3.3]
### Fixed
- Coverage plots from Tiddit data was printed strange

## [0.3.2]
### Added
- Add Support for Hetero/Isodisomy of UPD Region Files.

## [0.3.1]
### Fixed
- Bugfix: outd ignored for upd sites

## [0.3.0]
### Added
- Change flag "upd" to "sites
### Fixed
- Refactor code

## [0.2.1]
### Fixed
- When using `euploid` flag written empty png:s did not match correct genome prefix (37 vs 38)

## [0.2.0]
### Added
- Option `--euploid`
- More Pytest tests

### Fixed
- Better follows Pylint (though far from perfect)


## [0.1.3]
### Added
- Support for second guessing ideogram bed file chromosome format (i.e. 1 or chr1)
- Test file `test_chr_utils.py`

### Fixed
- Bug: Sloping graphs
- Bug: use `os.path` instead of concatenating string

## [0.1.2]

### Added
- Improved args handling

### Fixed
- Bug: Remove dependency on X display



## [0.1.1]
### Fixed
- Configured and parameterized WIG step was previously ignored
- Bug when sometimes malplaced linear slopes appeared in coverage pictures

## [0.1]
First version.
