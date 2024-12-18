# Pseudo-batch transformation article
This folder contains all information that is required to reproduce the results of the article. We provide a Docker image that contains a copy of the full github repository. Inside this image the simulations and data analysis can be redone. Because the purpose of the Docker image is to reproduce the results in the article the version of the pseudobatch code will be locked at the state of publication.

## This folder
- **figures** - figures used for the publication and extra figures for presentations.
- **julia-env** - files that describe the required packages for Julia scripts.
- **notebooks** - jyputer notebooks that reproduces the results and plots described in the paper.
- **simulation_scripts** - Julia scripts used to generate the simulated data.
- **data** - contains the data both simulated and real world data used in the article
- **requirements.txt** - specifies the exact versions of all python packages required to reproduce the results of the paper. The main use is for building the Docker image.
- **run_all.sh** - Shell script to reproduce the simulated data. This script should be executed from within the Docker container.


## How to reproduce results
### 0. Install docker on your local machine
See how to install Docker on your machine at the [Docker website](https://docs.docker.com/get-docker/).

### 1. Building docker image 
Currently, the publication of the Docker image is under review in the data repository (data.dtu.dk). Until, the image can be made public you can build the image from yourself. 

The Docker container relies on the [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook/tags/), see also [here](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html#jupyter-datascience-notebook). We used this image because it includes both Python and Julia out of the box, the down side is that it is quite large. To create the image you need a local copy of the pseudobatch repository and inside the repository root folder run the following command (Be sure to have the docker client running):

```
docker build . -t pseudobatch:1.2
```

This recreates the docker image, expect that it take ~ 15 - 30 min to build.

### 2. Start the docker container
- Ensure that the image is correctly build by running `docker image ls`. An image named *pseudobatch:1.2* should be found on that list.
- Run `docker run --rm -it -p 8888:8888 "pseudobatch:1.2"`
- Open the container using the instructions printed in the terminal

### 3. Rerun the simulations
We used several simulated datasets to test and showcase the pseudobatch transformation. This simulated datasets can be recreated by running a series of Julia scripts. To simplify the process they can be all be run at once using the shell script `article/run_all.sh`. Please use the following commands

1. Open a terminal inside the Docker container
2. Navigate into the `article` folder, using `cd article`
3. Run the shell script using the command `./run_all.sh`

This will take 10 - 20 minutes because the all the simulation dependencies has to be installed.

### 4. Rerun analysis
To rerun the analysis simply open the notebooks, **CHANGE THE KERNEL** to `venv_pseudobatch` (do this in the upper right corner), and run all cells. The error propagation analysis will take a while because it needs first has to compile the stan model.


