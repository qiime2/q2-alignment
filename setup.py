# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-alignment",
    # TODO stop duplicating version string
    version="0.0.0-dev",
    packages=find_packages(),
    install_requires=['scikit-bio', 'qiime >= 2.0.0', 'q2-types'],
    #package_data={'q2_alignment': ['markdown/*md']},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Create and work with alignments in QIIME 2.",
    license="BSD",
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugins': ['q2-alignment=q2_alignment.plugin_setup:plugin']
    }
)
