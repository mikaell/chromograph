from setuptools import setup

setup(
   name='chromograph',
   version='0.1',
   description='tool for plotting genetic data',
   author='mikaell',
   author_email='mikael.laaksonen@scilifelab.se',
   packages=['chromograph'],  #same as name
   install_requires=['pandas', 'matplotlib', 'numpy'],
)
