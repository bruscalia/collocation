from setuptools import setup

setup(
  name = 'collocation',
  packages = ['collocation'],
  version = '0.1.0',
  license='MIT',
  description = 'A repository with a python implementation of Orthogonal Collocation by Villadsen and Stewart (1967) for solving second order boundary value problems.',
  author = 'Bruno Scalia C. F. Leite',
  author_email = 'bruscalia12@gmail.com',
  url = 'https://github.com/bruscalia/collocation',
  download_url = 'https://github.com/bruscalia/collocation',
  keywords = ['Orthogonal Collocation',
              'Boundary Value Problems',
              'Differential Equations',
              'Numerical Methods'],
  install_requires=[
          'numpy>=1.19.0',
          'scipy>=1.7.0',
      ],
)