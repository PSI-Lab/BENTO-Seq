from setuptools import setup
from bootstrap_psi.version import version

setup(name='Bootstrap PSI',
      version=version,
      description='Estimate the percent spliced-in ratio (PSI) '
      'of alternative splicing events using bootstrap.',
      author='Hannes Bretschneider',
      author_email='hannes@psi.utoronto.ca',
      packages=['bootstrap_psi'],
      scripts=['bin/bootstrap-psi'],
      install_requries=['pysam (>=0.7.8)', 'numpy (>=1.6)']
    )
