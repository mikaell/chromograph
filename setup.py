import codecs
from setuptools import setup, find_packages

setup(
    name='chromograph',
    version='0.1',
    description='tool for plotting genetic data',
    author='mikaell',
    author_email='mikael.laaksonen@scilifelab.se',
    packages=find_packages(exclude=('tests*', 'docs', 'examples')),
    include_package_data=True,
    zip_safe=False,
    install_requires=['pandas', 'matplotlib', 'numpy', 'pyaml'],
    entry_points=dict(
        console_scripts=[
            'chromograph = chromograph.chromograph:main',
        ],
    ),

)
