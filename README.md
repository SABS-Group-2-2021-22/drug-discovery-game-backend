[![codecov](https://codecov.io/gh/SABS-Group-2-2021-22/drug-discovery-game-backend/branch/main/graph/badge.svg?token=8521K2DNMB)](https://codecov.io/gh/SABS-Group-2-2021-22/drug-discovery-game-backend)
[![Documentation Status](https://readthedocs.org/projects/drug-discovery-backend/badge/?version=latest)](https://drug-discovery-backend.readthedocs.io/en/latest/?badge=latest)

# drug-discovery-game-backend
## Documentation 
Documentation can be found here: [https://drug-discovery-backend.readthedocs.io](https://drug-discovery-backend.readthedocs.io)

## Installation
You must have miniconda installed to create the python environment.
Download the repository and then execute install.sh file in the directory you want the drug discovery game to be installed into. To download and run the install.sh file, run the following code in your command terminal:

For Linux or Windows subsystem Linux:
```bash
curl -O https://raw.githubusercontent.com/SABS-Group-2-2021-22/drug-discovery-game-backend/225-adding-installation-file-for-easy-install/install_Linux.sh && chmod +x install_Linux.sh && ./install_Linux.sh
```

For OSX:
```bash
curl -O https://raw.githubusercontent.com/SABS-Group-2-2021-22/drug-discovery-game-backend/225-adding-installation-file-for-easy-install/install_OSX.sh && chmod +x install_OSX.sh && ./install_OSX.sh
```


This will clone both the backend and front end repositories into a drug-discovery-game repository and install the dependencies.

## Running the backend locally

To run the game, in separate terminals, run:

```bash
./run_backend.sh
```


in the backend directory and 

```bash
./run_frontend.sh
```

 in the frontend directory.

 To run the game with the LLM functionality, add the -LLM y flag to the run commands:

 ```bash
./run_backend.sh -LLM y
```


in the backend directory and 

```bash
./run_frontend.sh -LLM y
```

 in the frontend directory.




