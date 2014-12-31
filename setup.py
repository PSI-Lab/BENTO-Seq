from setuptools import setup
from bento_seq.version import version

setup(name='BENTO-Seq',
      version=version,
      description='Estimate the percent spliced-in ratio (PSI) '
      'of alternative splicing events using bootstrap.',
      author='Hannes Bretschneider',
      author_email='hannes@psi.utoronto.ca',
      packages=['bento_seq'],
      scripts=['bin/bento-seq'],
      install_requires=['pysam', 'numpy']
    )
