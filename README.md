# Meshfree Surface Integration

This repository contains code for the numerical tests in "High-Order Meshfree Surface Integration, Including Singular Integrants". A 4.8 million element mesh of the surface used for the tests of Subsection 4.1 is also included. Please see the cited paper listed under "Citation" for a more detailed description of the method and tests.

## File Names
The test files for each Subsection are named so that testXpYpZ.m corresponds to Subsection X.Y.Z.
fig1_code.m and fig4_code.m contain the data and scripts necessary to produce Figure 1 and 4, respectively.

place_b_points.m is a script for generating points on a surface described by a level set.

## Instructions
The test files are set up to generate a single value corresponding to one point distribution with one set of parameters. Each file has the parameters necessary to generate each value in a plot or table, with comments indicating which parameters need to be changed.

## Citation

    @article{vennruuth2025,
        title={High-Order Meshfree Surface Integration, Including Singular Integrands}, 
        author={Daniel R. Venn and Steven J. Ruuth},
        year={2025},
        eprint={},
        archivePrefix={arXiv},
        primaryClass={math.NA},
        url={}
    }
