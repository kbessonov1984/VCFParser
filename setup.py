from setuptools import setup
from setuptools import find_packages

setup(name='VCFParser',
      version='0.0.1',
      python_requires='>=3.7.0,<4',
      description='Updateable VCF parser for COVID-19 Variant call files',
      author='Kyryl Bessonov',
      license='GPLv3',
      install_requires=['pandas', 'matplotlib', 'numpy', 'pysam'],
      entry_points={'console_scripts':['vcfparser=vcfparser.vcfparser:main']},
      python_require='>=3.7,<4',
      include_package_data=True,
      packages=find_packages(exclude=['tests']),
      tests_require=['pytest']

)
