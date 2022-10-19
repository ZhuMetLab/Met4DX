# Met4DX

## Introduction
`Met4DX` is an R package for timsPro 4D data processing.

The docker image zhulab/mettracer-r contains entire envorienment for running `Met4DX`. For convinience and taking fully use of `Met4DX`, users can pull it and run `Met4DX` just as following.

## What is Met4DX

`MetT4DX-r` is an Docker environment to processing isotope labelled metabolomics data with #Met4DX R package. It is based on the [`r-base`](https://hub.docker.com/_/r-base/) docker.

## Pulling image

Users can pull the MetTracer-r image with the following script

```bash
docker pull zhulab/met4dx-r
```

## Data preparation


### Preparing the R script

File run.R is an R source code file in  data folder. We provide a template here. Users only need to change the folder name in general. Other parameters are recommended parameters.

- Parameters

```R

```

### Overview of the data preparation

## Run data processing work with mettracer-r image

- go to your data folder (e.g., data)

```base
 cd data
```

- run docker using following code (*User should be permited to run docker service*)

```bash
# MUST keep the code exactly as it is!
docker run -it --rm -v "$PWD":/data -u $(id -u ${USER}):$(id -g ${USER}) zhulab/met4dx-r Rscript run.R
```

- wait till data processing work done

- Explaining `docker run` arguments

- `-v "$PWD":/home/${USER}`: mapping current dirctory as home directory in docker container

- `-u $(id -u ${USER}):$(id -g ${USER})`: using current user to run the container

- `Rscript ~/run.R`: run run.R in container home directory with `Rscript`  command

# License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a>

This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
