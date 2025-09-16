Installation
------------

The full python package can be installed using `pip`. 

A small number of compiled packages are required (CASA and MongoDB). To make life easier, both the Python and other dependencies can be installed using `conda`.

First, ensure you have  `conda <https://github.com/conda-forge/miniforge>`_ installed (miniforge is recommended).

.. tip::
    I strongly recommend using `mamba <https://github.com/mamba-org/mamba>`_ to install the dependencies. Mamba is essentially a compiled and parallelized version of conda. It will greatly speed up your installation.

After cloning this repo, please run: ::

    cd arrakis/

    conda env create
    # or - if you have mamba:
    mamba env create -f environment.yml

This will install the non-dependencies, including `uv`, and create a virtual environment called `arrakis310`

Then you can install the Python dependencies and the command-line scripts into the environment using `pip` or `uv` e.g. ::

    conda activate arrakis310
    # Recommended - uses the provided lock file
    uv sync
    # or - flexible install
    uv pip install .
    # or - same as above but slower
    pip install .

An installation of Singularity is also required.

.. attention::

   The version of Singularity that is required currently on conda-forge appears to be broken. For now, I will not include it in the conda environment. Please source it elsewhere for the time being.

That's it!
