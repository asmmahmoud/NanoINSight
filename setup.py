from setuptools import setup, find_packages
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

exec(open("tmpname/version.py").read())

setup(
    name='tmpname',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/AsmaaSamyMohamedMahmoud/tmpname',
    license='gpl-3.0',
    author='Asmaa',
    author_email='asmmahmoud@mun.ca',
    description='Repeat annotation tool for insertions called by NanoVar',
    keywords=['insertion structural variant repeat annotation'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['biopython>=1.74', 'scikit-allel>=1.3.7'],
    entry_points={
        "console_scripts": [
            "tmpname=tmpname.tmpname:main",
        ],
    },
    python_requires='>=3',
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)