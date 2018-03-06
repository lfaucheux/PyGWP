from distutils.core import setup
from PyGWP import __version__ as v

setup(
    name='PyGWP',
    packages = ['PyGWP'],
    package_data={
        'PyGWP': [
            'test/*',
        ]
    },
    version=v,
    description="A CO2-equivalent computer based on static or dynamic CO2-relative global warming potentials.",
    author='Laurent Faucheux',
    author_email="laurent.faucheux@hotmail.fr",
    url='https://github.com/lfaucheux/PyGWP',
    download_url = 'https://github.com/lfaucheux/PyGWP/archive/{}.tar.gz'.format(v),
    keywords = ['scientific modelling', 'gwp', 'ghg'],
    classifiers=[],
)
