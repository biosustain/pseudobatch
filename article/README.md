# Pseudo-batch transformation article
This folder contains all information that is required to reproduce the results of the article. We provide a Docker image that contains a copy of the full github repository. Inside this image the simulations and data analysis can be redone.

## This folder
- **figures** - figures used for the publication and extra figures for presentations.
- **julia-env** - files that describe the required packages for Julia scripts.
- **notebooks** - jyputer notebooks that reproduces the results and plots described in the paper.
- **scripts** - Julia scripts used to generate the simulated data.
- **simulated_data** - contains the data both simulated and real world data used in the article


## How to reproduce results
### 1. Download the docker image
- The docker image can be found here: XXX
- load the image to your local docker by running `docker load --input PATH-TO-IMAGE-TAR-FILE`

### 2. Start the docker container
- Ensure that the image is correctly loaded by running `docker image ls`. An image named *pseudobatch:1.0* should be found on that list.
- Run `docker run --rm -it -p 8888:8888 "pseudobatch:1.0"`
- Open the container using the instructions printed in the terminal

### 3. Rerun the simulations
The simulated data that is used to test and show case the pseudo-batch transformation can be recreated by running a series of Julia scripts. To simplify the process they can be all be run at once by opening a terminal and executing `./run_all.sh` (this should be done INSIDE the docker container).

### 4. Rerun analysis
To rerun the analysis simply open the notebooks and run them.
