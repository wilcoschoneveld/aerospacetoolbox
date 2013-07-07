from setuptools import setup

setup(name = 'aerospacetoolbox',
      version = '0.6',
      description = 'Aerospace Toolbox',
      author = 'Wilco Schoneveld',
      author_email = 'schoneveld.wj@gmail.com',
      url = 'https://github.com/wilcoschoneveld/aerospacetoolbox',
      requires = ['scipy'],
      #package_dir = {'aerospacetoolbox' : 'src'},
      packages = ['aerospacetoolbox'],
     )
