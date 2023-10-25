from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = 'v0.3.0'
DESCRIPTION = 'Molecular Docking library for Python'
#LONG_DESCRIPTION = ''

# Setting up
setup(
    name="JupyterDock",
    version=VERSION,
    author="Angel J. Ruiz-Moreno",
    author_email="<angel.j.ruiz.moreno@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=['rdkit', 'datamol', 'pyopenms', 'PubChemPy','pandas', 'plotly', 'mols2grid', 'rxnmapper', 'distinctipy', 'numpy', 'matplotlib', 'tqdm'],
    keywords=['python', 'metabolite', 'prediction', 'microbiome', 'cheminformatics', 'metabolism'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
    #include_package_data=True,
    #package_data={'': ['DataBase/*.tsv.gz']},
)