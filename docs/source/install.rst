Installation
------------

The full python package can be installed using `pip`.

To make life easier, both the Python and other dependencies can be installed using conda.

First, ensure you have `anaconda <https://www.anaconda.com/products/individual/>`_ or `miniconda <https://docs.conda.io/en/latest/miniconda.html/>`_ installed.

.. tip::
    I strongly recommend using `mamba <https://github.com/mamba-org/mamba>`_ to install the dependencies. Mamba is essentially a compiled and parallelized version of conda. It will greatly speed up your installation.

After cloning this repo, please run: ::

    cd arrakis/

    conda env create
    # or - if you have mamba:
    mamba env create

This will install the python dependencies and the command-line scrips into a conda environment called `arrakis310`, which can be activated by: ::

    conda activate arrakis310

An installation of Singularity is also required.

.. attention::

   The version of Singularity that is required currently on conda-forge appears to be broken. For now, I will not include it in the conda environment. Please source it elsewhere for the time being.

That's it!
