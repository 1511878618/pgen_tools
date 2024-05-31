# coding:utf-8
# Copyright (c) 2023  Tingfeng Xu. All Rights Reserved.

from setuptools import find_packages
from setuptools import setup

import pgenlib_tools
import os 
from pathlib import Path
script_path = os.path.dirname(os.path.abspath(__file__)) + "/pgenlib_tools/scripts"
print(script_path)
scripts = [str(i) for i in Path(script_path).rglob("*.py") if "pgenlib_tools" in str(i)]
print(scripts)

with open("requirements.txt") as file:
    REQUIRED_PACKAGES = file.read()

setup(
    name='pgenlib_tools',
    version=pgenlib_tools.__version__.replace('-', ''),
    description=('pgenlib_tools'),
    long_description='',
    # url='https://github.com/1511878618/cadFace',
    author='Tingfeng Xu',
    author_email='xutingfeng@big.ac.cn',
    install_requires=REQUIRED_PACKAGES,
    packages=find_packages(),
    scripts=scripts

    )
