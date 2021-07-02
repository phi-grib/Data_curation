from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='data-curation',
    version='2.0',
    description='Data curation package. Optimized for pandas Dataframes',
    license='GNU',
    long_description=long_description,
    author='Eric March Vila',
    author_email='eric.march@upf.edu',
    url='https://github.com/phi-grib/Data_curation',
    download_url='https://github.com/phi-grib/Data_curation.git',
    packages=find_packages(),
    # If any package contains *.txt or *.rst files, include them:
    # package_data={'': ['*.yaml', '*.yml']},
    package_data={'data-curation': ['curate/config.yaml','curate/children/*.yaml']},
    entry_points={
       'console_scripts': ['datacur=curate.curate_src:main'],
    }
)
