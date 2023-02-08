# This script will run a local version of the backend.
# If you have named your conda environment differently please change the environment name from 'dd_game'to something else
# If you are using a different port for the backend API change the API_URL
# Do not use this script for deployment, only for local development

# export PERMITTED_FRONTEND_ORIGIN='http://127.0.0.1:3000'
export PERMITTED_FRONTEND_ORIGIN='http://localhost:3000'
export FLASK_APP=api/service
flask run
