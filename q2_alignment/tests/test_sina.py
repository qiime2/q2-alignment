# ----------------------------------------------------------------------------
# Copyright (c) 2018, Elmar Pruesse.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import unittest

from q2_alignment import sina

from q2_types.feature_data import AlignedDNAFASTAFormat, DNAFASTAFormat

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio

import skbio


if __name__ == "__main__":
    unittest.main()

class SINATests(TestPluginBase):
    package = 'q2_alignment.tests'
    
