# Batch-effect-removal-benchmarking
Datasets and scripts used in article entitled "A benchmark of batch-effect correction methods for single-cell RNA sequencing data" (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9).

The same files are also available at https://hub.docker.com/repository/docker/jinmiaochenlab/batch-effect-removal-benchmarking


## Docker installation

Please confirm that:

1. You have installed docker as the instructions: [Install Docker Desktop on Mac](https://docs.docker.com/docker-for-mac/install/), [Install Docker Desktop on Windows](https://docs.docker.com/docker-for-windows/install/) and [Docker engine for Linux](https://docs.docker.com/install/) according to your OS.

2. You are in "docker" or "sudo" group if you are using Mac or Linux OS.

## Docker usage

To download the data, please follow these steps:

1. Go to target folder
2. run ```sudo docker create --name batcheffect  jinmiaochenlab/batch-effect-removal-benchmarking``` to download docker image.
3. run ```sudo docker cp batcheffect:/batch_effect/ ./``` to copy the data to current folder.
4. run  ```sudo docker container rm batch_effect; docker image rm jinmiaochenlab/batch-effect-removal-benchmarking``` to delete docker image.
4. If you are in "docker" group, then you don't need to use sudo.

### If you encounter any issue to install or run docker, please consult your system administrator. It will save you a lot of time.
