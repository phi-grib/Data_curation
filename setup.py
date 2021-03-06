from setuptools import setup, find_packages

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='data_curation',
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
    package_data={'': ['*.yaml', '*.yml']},
    entry_points={
       'console_scripts': ['datacur=curate.curate_src:main'],
    }
)
