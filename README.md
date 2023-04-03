# pseudobatch

## Building the documentation

Pseudobatch's documentation lives in the folder `docs`, written in Jupyter
notebooks and restructured text documents and built using
[Sphinx](https://www.sphinx-doc.org). 

The source files can be found in the folder `docs/source`.

In order to rebuild the documentation after editing, first make sure that you
have installed all the dependencies by running `pip install -e
.'[development]'` from the project root. Next, change directory to `docs` and
run the command `make html`. To view your changes run `open
build/html/index.html` or just click through to this file using your file
explorer.
