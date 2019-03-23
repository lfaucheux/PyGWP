# -*- coding: utf8 -*-
from __future__ import absolute_import
from setuptools import setup
import pathlib as pa
import codecs as cd

package_version = "0.0.30"
package_name    = 'PyGWP'

with pa.Path('requirements.txt').open() as requirements:
    requires = [l.strip() for l in requirements]

with cd.open('README.md', encoding='utf-8') as readme_f:
    readme = readme_f.read()

setup(
    license      = 'MIT',
    name         = package_name,
    version      = package_version,
    packages     = [package_name],
    package_data = {
        package_name: [
            'test/*',
        ]
    },
    description =(
        "A CO2-equivalent computer based on static or dynamic "
        "CO2-relative global warming potentials coded in Python27, %s."%package_name
    ),
    long_description              = readme,
    long_description_content_type = 'text/markdown',
    author       = 'Laurent Faucheux',
    author_email = "laurent.faucheux@hotmail.fr",
    url          = 'https://github.com/lfaucheux/%s'%package_name,
    download_url = 'https://github.com/lfaucheux/{}/archive/{}.tar.gz'.format(
        package_name, package_version
    ),
    classifiers  = [
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
    ],
    keywords = [
        'Global Warming Potential',
        'Static Global Warming Potential',
        'Dynamic Global Warming Potential'
    ],
    install_requires = requires,
)
