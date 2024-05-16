import shutil
import warnings
from pathlib import Path
from importlib.metadata import distribution

try:
    distribution("cmdstanpy")
    import cmdstanpy

    cmdstanpy_installed = True
except:
    cmdstanpy_installed = False

if cmdstanpy_installed:
    from pseudobatch.error_propagation.error_propagation import (
        run_error_propagation,
    )

    STAN_FILES_FOLDER = Path(__file__).parent / "stan"
    CMDSTAN_VERSION = "2.31.0"

    # on Windows specifically, we should point cmdstanpy to the repackaged
    # CmdStan if it exists. This lets cmdstanpy handle the TBB path for us.
    local_cmdstan = STAN_FILES_FOLDER / f"cmdstan-{CMDSTAN_VERSION}"
    if local_cmdstan.exists():
        cmdstanpy.set_cmdstan_path(str(local_cmdstan.resolve()))

    def load_stan_model(name: str) -> cmdstanpy.CmdStanModel:
        """
        Try to load precompiled Stan models. If that fails,
        compile them.
        """
        try:
            model = cmdstanpy.CmdStanModel(
                exe_file=STAN_FILES_FOLDER / f"{name}.exe",
                stan_file=STAN_FILES_FOLDER / f"{name}.stan",
                compile=False,
            )
        except ValueError:
            warnings.warn(
                f"Failed to load pre-built model '{name}.exe', compiling"
            )
            model = cmdstanpy.CmdStanModel(
                stan_file=STAN_FILES_FOLDER / f"{name}.stan",
                stanc_options={"O1": True},
            )
            shutil.copy(
                model.exe_file,  # type: ignore
                STAN_FILES_FOLDER / f"{name}.exe",
            )

        return model

    ERROR_PROPAGATION = load_stan_model("error_propagation")

    # example: just print the info of the model
    print(ERROR_PROPAGATION.exe_info())
