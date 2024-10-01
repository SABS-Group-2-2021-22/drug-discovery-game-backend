
#!/bin/bash

# Define repositories
mkdir Drug-discovery-game
cd Drug-discovery-game

BACKEND_REPO="git@github.com:SABS-Group-2-2021-22/drug-discovery-game-backend.git"
FRONTEND_REPO="git@github.com:SABS-Group-2-2021-22/drug-discovery-game-app.git"

# Clone repositories
echo "Cloning repositories..."
git clone $BACKEND_REPO
git clone $FRONTEND_REPO

# Setup Backend
echo "Setting up the backend..."
cd drug-discovery-game-backend
conda env create -f environment.yml


# Create run_backend.sh
echo "Creating run_backend.sh..."
cat <<EOF > run_backend.sh
#!/bin/bash
# Activate backend environment and run
cd drug-discovery-game-backend
conda activate dd_game
export FLASK_APP=api/service
flask run -p 8000
EOF
chmod +x run_backend.sh

cd ..



# Setup Frontend
echo "Setting up the frontend..."
cd drug-discovery-game-app
npm install


# Create run_frontend.sh
echo "Creating run_frontend.sh..."
cat <<EOF > run_frontend.sh
#!/bin/bash
# Run frontend server
cd drug-discovery-game-app
npm start
EOF
chmod +x run_frontend.sh
cd ..



echo "Installation completed. You can now run './run_backend.sh' and './run_frontend.sh' in separate terminals within their respective directories."