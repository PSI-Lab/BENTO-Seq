from setuptools import setup
from bootstrap_tool.version import version

setup(name='Bootstrap PSI',
      version=version,
      description='Estimate the percent spliced-in ratio (PSI) '
      'of alternative splicing events using bootstrap.',
      author='Hannes Bretschneider',
      author_email='hannes@psi.utoronto.ca',
      packages=['bootstrap_tool'],
      scripts=['bootstrap_psi.py'],
      requries=['pysam (>=0.7.8)', 'numpy (>=1.7)']
    )


      
