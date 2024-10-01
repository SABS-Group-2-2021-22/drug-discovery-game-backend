
#!/bin/bash

# Define repositories
mkdir Drug-discovery-game
cd Drug-discovery-game

BACKEND_REPO="https://github.com/SABS-Group-2-2021-22/drug-discovery-game-backend.git"
FRONTEND_REPO="https://github.com/SABS-Group-2-2021-22/drug-discovery-game-app.git"

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
# Default set to no
LLM_FLAG="n"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -LLM) LLM_FLAG="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check the LLM flag
if [ "$LLM_FLAG" = "y" ]; then
    echo "Switching to main_LLM branch..."
    git checkout 222-openai-api-integration-BE
else
    echo "Using main branch..."
    git checkout main
fi

# Initialize conda
eval "$(conda shell.bash hook)"

# Activate the environment
conda activate dd_game

# Check if activation was successful
if [ $? -ne 0 ]; then
    echo "Failed to activate Conda environment. Please ensure the environment 'dd_game' exists."
    exit 1
fi

export FLASK_APP=api/service
flask run -p 8000
EOF
chmod +x run_backend.sh

cd ..



# Setup Frontend
echo "Setting up the frontend..."
cd drug-discovery-game-app
sudo apt install nodejs
sudo apt install npm
npm install


# Create run_frontend.sh
echo "Creating run_frontend.sh..."
cat <<EOF > run_frontend.sh
#!/bin/bash
# Run frontend server

# Default set to no
LLM_FLAG="n"

# Parse arguments
while [[ "\$#" -gt 0 ]]; do
    case \$1 in
        -LLM) LLM_FLAG="\$2"; shift ;;
        *) echo "Unknown parameter passed: \$1"; exit 1 ;;
    esac
    shift
done

# Check the LLM flag
if [ "\$LLM_FLAG" = "y" ]; then
    echo "Switching to main_LLM branch..."
    git checkout OPENCHAT_with_openai
else
    echo "Using main branch..."
    git checkout main
fi


npm start
EOF
chmod +x run_frontend.sh
cd ..



echo "Installation completed. You can now run './run_backend.sh' and './run_frontend.sh' in separate terminals within their respective directories."