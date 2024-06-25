# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-pinocchio",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Christos Matzoros",
    author_email="christosmatzoros@gmail.com",
    description="Plugin for quality control and taxonomic "
    "classification of long-read sequencing data.",
    url="https://github.com/bokulich-lab/q2-pinocchio",
    entry_points={"qiime2.plugins": ["q2-pinocchio=q2_pinocchio.plugin_setup:plugin"]},
    package_data={
        "q2_pinocchio": [
            "citations.bib",
            "assets/*/*",
        ],
        "q2_pinocchio.tests": [
            "data/*",
            "data/*/*",
            "data/*/*/*",
        ],
        "q2_pinocchio.types": ["data/*"],
    },
    zip_safe=False,
)
