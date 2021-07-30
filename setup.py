from setuptools import setup
from setuptools import find_packages
from vcfparser import __version__

setup(name='VCFParser',
      version=__version__,
      python_requires='>=3.7.0,<4',
      description='Updateable VCF parser for COVID-19 Variant Call Files (VCFs)',
      author='Kyryl Bessonov',
      license='GPLv3',
      install_requires=['pandas', 'matplotlib>=3.3', 'numpy', 'pysam', 'openpyxl'],
      entry_points={'console_scripts':['vcfparser=vcfparser.vcfparser:main']},
      python_require='>=3.7,<4',
      include_package_data=True,
      packages=find_packages(exclude=['tests']),
      tests_require=['pytest']

)
