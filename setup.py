# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-alignment",
    version="2017.3.0.dev",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.3.*', 'q2-types == 2017.3.*',
                      'scikit-bio'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Create and work with alignments in QIIME 2.",
    license="BSD-3-Clause",
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins': ['q2-alignment=q2_alignment.plugin_setup:plugin']
    },
    package_data={'q2_alignment.tests': ['data/*']}
)
