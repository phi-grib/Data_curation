from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name='data_curation',
   version='1.0',
   description='Data curation package. Optimized for pandas Dataframes',
   license='GNU',
   long_description=long_description,
   author='Eric March Vila',
   author_email='eric.march@upf.edu',
   url='https://github.com/phi-grib/Data_curation',
   packages=['curate']
)
