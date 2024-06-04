Installation
============
The Pseudobatch Python package can be install through PYPI using pip. Most of the functionality can be installed simply by calling

.. code-block:: bash

    pip install pseudobatch

The error propagation functionality requires installation of cmdstanpy and CmdStan. Please see the cmdstanpy documentation for how to install both cmdstanpy and CmdStan in your specific setup (https://mc-stan.org/cmdstanpy/installation.html). After successfully installing both proceed to install the remaining dependencies of the error propagation module through pip.

.. code-block:: bash

    pip install pseudobatch[error_propagation]

Now the error propagation module is installed and ready to use. Note that the first time you import the error propagation module CmdStan will compile the Stan model this will take several minutes.

