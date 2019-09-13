import codecs
import os
import io

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, 'chromograph', '__version__.py')) as f:
    exec(f.read(), about)

def parse_reqs(req_path='./requirements.txt'):
    """Recursively parse requirements from nested pip files."""
    install_requires = []
    with io.open(os.path.join(here, 'requirements.txt'), encoding='utf-8') as handle:
        # remove comments and empty lines
        lines = (line.strip() for line in handle
                 if line.strip() and not line.startswith('#'))

        for line in lines:
            # check for nested requirements files
            if line.startswith('-r'):
                # recursively call this function
                install_requires += parse_reqs(req_path=line[3:])

            else:
                # add the line as a new requirement
                install_requires.append(line)

    return install_requires

# What packages are required for this module to be executed?
REQUIRED = parse_reqs()

setup(
    name='chromograph',
    version=about['__version__'],
    description='tool for plotting genetic data',
    author='mikaell',
    author_email='mikael.laaksonen@scilifelab.se',
    packages=find_packages(exclude=('tests*', 'docs', 'examples')),
    include_package_data=True,
    zip_safe=False,
    install_requires=REQUIRED,
    entry_points=dict(
        console_scripts=[
            'chromograph = chromograph.chromograph:main',
        ],
    ),

)
