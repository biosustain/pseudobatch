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

## Building docker image 
This is a note to developers who what to rebuild/update the docker image. If you simply want to use the docker image see the description in the [article folder](./article/README.md). 

The Docker container relies on the [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook/tags/), see also [here](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html#jupyter-datascience-notebook). We used this image because it includes both Python and Julia out of the box, the down side is that it is quite large. To create the image you need a local copy of this repository and inside that folder run the following command:

```
docker build . -t pseudobatch:{version}
```

This recreates the docker image. Be aware that installing the Julia packages and cmdstan are both take quite some time, thus expect that it take ~ 15 - 30 min to build. For that exact reason the docker-compose specifies the current folder as a volume.