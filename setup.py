#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Kathryn Materna
# Copyright 2021. ALL RIGHTS RESERVED.
# United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import re
from setuptools import find_packages, setup

# Parameter defs
CWD = os.getcwd()


def get_version():
    with open('version.txt', 'r') as f:
        m = re.match("""version=['"](.*)['"]""", f.read())

    assert m, "Malformed 'version.txt' file!"
    return m.group(1)


setup(
    name='Elastic_stresses_py',
    version=get_version(),
    description='This is the Elastic_stresses_py package',
    package_dir={
        'PyCoulomb': 'PyCoulomb',
        '': 'PyCoulomb'
    },
    packages=['PyCoulomb'] + find_packages('PyCoulomb'),
    scripts=[
        'PyCoulomb/bin/disp_okada_driver.py',
        'PyCoulomb/bin/elastic_stresses_config_writer.py',
        'PyCoulomb/bin/elastic_stresses_driver.py'
    ],
    zip_safe=False,
)
