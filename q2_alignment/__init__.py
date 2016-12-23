# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._mafft import mafft
from ._filter import mask

__version__ = pkg_resources.get_distribution('q2-alignment').version

__all__ = ['mafft', 'mask']
