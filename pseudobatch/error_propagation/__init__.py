import shutil
import warnings
from pathlib import Path
from importlib.metadata import distribution

# The error propagation module requires to be installed separately.
# Cmdstanpy is here used as an indication of whether the error
# propagation module is installed or not. If cmdstanpy is installed,
# the error propagation module will be loaded.
try:
    distribution("cmdstanpy")
    import cmdstanpy

    cmdstanpy_installed = True
except:
    cmdstanpy_installed = False
    raise ImportError(
        "cmdstanpy is not installed. To use the error propagation module, "
        "please install cmdstanpy and CmdStan. You can find installation instructions "
        "for different setups at https://mc-stan.org/cmdstanpy/installation.html."
        "After installing cmdstanpy and CmdStan install the remaining dependencies "
        "for the error propagation module by running "
        "'pip install pseudobatch[error_propagation]'."
    )

if cmdstanpy_installed:
    from pseudobatch.error_propagation.error_propagation import (
        run_error_propagation,
    )
