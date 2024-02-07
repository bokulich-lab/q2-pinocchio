# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-long-reads-qc",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Christos Matzoros",
    author_email="christosmatzoros@gmail.com",
    description="QIIME 2 Plugin for quality control of long read sequences.",
    url="https://github.com/bokulich-lab/q2-long-reads-qc",
    entry_points={"qiime2.plugins": ["q2-long-reads-qc=q2_long_reads_qc.plugin_setup:plugin"]},
    zip_safe=False,
)
