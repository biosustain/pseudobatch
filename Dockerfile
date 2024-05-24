# Start from a core stack version
FROM quay.io/jupyter/datascience-notebook:2024-05-20

ENV CMDSTAN_VERSION=2.31.0

# Copy all files from directory to docker image
COPY --chown=${NB_UID}:${NB_GID} . .

# Install in the default python3 environment
# Install cmdstanpy
RUN pip install --no-cache-dir "cmdstanpy==1.0.4"

# Install cmdstan
RUN python -m cmdstanpy.install_cmdstan --version "${CMDSTAN_VERSION}" --cores 2

RUN pip install --no-cache-dir -e ".[error_propagation]" && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

RUN julia --project=julia-env -e 'using Pkg; pkg"activate"; pkg"precompile"' && \
    chmod -R go+rx "${CONDA_DIR}/share/jupyter" && \
    fix-permissions "${JULIA_PKGDIR}" "${CONDA_DIR}/share/jupyter"
