#!/bin/bash

# remove the old virtual environment if it exists
rm -rf .venv-docker

# Create a new virtual environment
python -m venv .venv-docker

# Activate the virtual environment
source .venv-docker/bin/activate

# Install the current directory
pip install -e '.[development]'

# Run the tests
pytest

# install error_propagation
pip install '.[error_propagation]'