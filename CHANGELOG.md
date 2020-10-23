# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
