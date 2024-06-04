# This Dockerfile is used to build the image in which the results
# of the paper can be reproduced. The image is based on the
# Jupyter Data Science Notebook image from Jupyter Docker Stack and 
# extends the image with CmdStan, and the Julia environment. 
# Furthermore the image uses the exact versions of the libraries, 
# that were used to produce the results in the paper.
# We want to thank the Jupyter Docker Stack Recepies as we got a 
# lot inspiration from the following recipe:
# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/recipes.html#add-a-custom-conda-environment-and-jupyter-kernel

# Start from a core stack version, locked the version to ensure reproducibility
# Base image: quay.io/jupyter/datascience-notebook:2024-05-20
# Licence: https://github.com/jupyter/docker-stacks/blob/main/LICENSE.md
############################## LICENCE FOR BASE IMAGE ##############################
# This project is licensed under the terms of the Modified BSD License (also known as New or Revised or 3-Clause BSD), as follows:

#     Copyright (c) 2001-2015, IPython Development Team
#     Copyright (c) 2015-, Jupyter Development Team

# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# Neither the name of the Jupyter Development Team nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#####################################################################################

FROM quay.io/jupyter/datascience-notebook:2024-05-20

# Set versions of core software
ENV CMDSTAN_VERSION=2.33.0
ENV CMDSTANPY_VERSION=1.1.0
# Name your environment and choose the Python version
ARG env_name=venv_pseudobatch
ARG py_ver=3.10.8

# Copy all files from directory to docker image
COPY --chown=${NB_UID}:${NB_GID} . .

# remove the precompile files for pseudobatch error_propagation
RUN rm -f \ 
    pseudobatch/error_propagation/stan/error_propagation \
    pseudobatch/error_propagation/stan/error_propagation.exe \
    pseudobatch/error_propagation/stan/error_propagation.hpp

# You can add additional libraries required for notebooks here
RUN mamba create --yes -p "${CONDA_DIR}/envs/${env_name}" \
    python=${py_ver} \
    'ipykernel' \
    'jupyterlab' && \
    mamba clean --all -f -y

# Install cmdstanpy and cmdstan inside the environment
RUN mamba run -p "${CONDA_DIR}/envs/${env_name}" mamba install --yes -c conda-forge cmdstanpy=="${CMDSTANPY_VERSION}" cmdstan=="${CMDSTAN_VERSION}" && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Create Python kernel and link it to jupyter
RUN "${CONDA_DIR}/envs/${env_name}/bin/python" -m ipykernel install --user --name="${env_name}" && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Installing the pseudobatch from local source using pip
RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install --no-cache-dir -e \
    '.[error_propagation]'

# Install specific requirements to match versions used when producing the paper
RUN "${CONDA_DIR}/envs/${env_name}/bin/pip" install --no-cache-dir -r article/requirements.txt

# This changes the custom Python kernel so that the custom environment will
# be activated for the respective Jupyter Notebook and Jupyter Console
# hadolint ignore=DL3059
RUN /opt/setup-scripts/activate_notebook_custom_env.py "${env_name}"

USER ${NB_UID}

# Setup Julia environment for simulations
RUN julia --project=julia-env -e 'using Pkg; Pkg.activate(); Pkg.instantiate()' && \
    chmod -R go+rx "${CONDA_DIR}/share/jupyter" && \
    fix-permissions "${JULIA_PKGDIR}" "${CONDA_DIR}/share/jupyter"
