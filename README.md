[![codecov](https://codecov.io/gh/SABS-Group-2-2021-22/drug-discovery-game-backend/branch/main/graph/badge.svg?token=8521K2DNMB)](https://codecov.io/gh/SABS-Group-2-2021-22/drug-discovery-game-backend)
[![Documentation Status](https://readthedocs.org/projects/drug-discovery-backend/badge/?version=latest)](https://drug-discovery-backend.readthedocs.io/en/latest/?badge=latest)

# drug-discovery-game-backend
## Documentation 
Documentation can be found here: [https://drug-discovery-backend.readthedocs.io](https://drug-discovery-backend.readthedocs.io)

## Installation
You must have miniconda installed to create the python environment.
After cloning the repository, install the dependencies with the following commands:
```
cd drug-discovery-game-backend
conda env create -f environment.yml
```

## Running the backend locally

This assumes the conda environment has been created and necessary packages installed.

```
conda activate dd_game
export FLASK_APP=api/service
flask run -p 8000
```
