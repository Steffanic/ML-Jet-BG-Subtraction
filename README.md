# ML-Jet-BG-Subtraction
C++ scripts for generating a realistic Jet event and python scripts for a machine learning solution to BG subtraction

# Installation
## Prerequisites
The generator code requires ROOT 5 which can be found [here](https://root.cern.ch/releases). Installation instructions are cryptic, but I would look [here](https://root.cern.ch/building-root) for info about that. Good luck!

The code also requires Pythia 6 which is included in the git, but there is no explicit gauruntee that it will properly connect with your system. We recommend looking at the maindriver_newton.C file in Analysis_Code/Generator to fix any user-specific filepaths.

The machine learning requires python3 and a couple of libraries that can be installed via pip

`pip3 install pandas matplotlib numpy sklearn`

