from setuptools import setup, find_packages

setup(
    name = 'aerospacetoolbox',
    version = '0.8',
    description = 'Aerospace Toolbox',
    author = 'Wilco Schoneveld',
    author_email = 'schoneveld.wj@gmail.com',
    url = 'https://github.com/wilcoschoneveld/aerospacetoolbox',
    install_requires = ['scipy'],
    packages = find_packages(exclude=['tests']),
    package_data = {
        '': ['egm96.dac']
    }
)
