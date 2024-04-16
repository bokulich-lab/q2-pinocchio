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
    name="q2-minimap2",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Christos Matzoros",
    author_email="christosmatzoros@gmail.com",
    description="QIIME 2 Plugin for quality control and taxonomic "
    "classification of long read sequences using Minimap2.",
    url="https://github.com/bokulich-lab/q2-minimap2",
    entry_points={"qiime2.plugins": ["q2-minimap2=q2_minimap2.plugin_setup:plugin"]},
    package_data={
        "q2_minimap2": ["citations.bib"],
        "q2_minimap2.tests": ["data/*", "data/consensus/*"],
        "q2_minimap2.types": ["data/*"],
    },
    zip_safe=False,
)
