from setuptools import setup, find_packages
import os
import textwrap

setup(
    name='fmmgen',
    version='0.0.1',
    description='Code generator for Fast Multipole and Barnes Hut operators',
    long_description=textwrap.dedent("""\
         fmmgen: FMM and Barnes-Hut code generator

         fmmgen aims to provide simple, efficient implementations of Cartesian FMM 
         operators which can be used in various projects.
 
         The Cartesian operators can difficult to implement by hand because
         of the large number of terms which occur. This project uses symbolic
         algebra to derive them, and these symbolic representations can be
         used directly, or ignored completely with only functions generated.
         """),
    author="Ryan Alexander Pepper",
    author_email="ryan.pepper@soton.ac.uk",
    url="https://github.com/rpep/fmmgen",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Programming Language :: C",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Linux",
        "Operating System :: MacOS"
        ],
    install_requires = ["cython", "pytest", "sympy", "numpy"],
    python_requires=">= 3.10",
)
