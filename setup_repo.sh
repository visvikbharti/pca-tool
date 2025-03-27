#!/bin/bash

# Create necessary directories
mkdir -p .streamlit
mkdir -p .github/workflows
mkdir -p demo_images

# Move files to the right locations if they exist
if [ -f ".streamlit/config.toml" ]; then
    echo "Config file already exists"
else
    echo "Moving config.toml to .streamlit folder"
    # Add config.toml content
    cat > .streamlit/config.toml << EOL
[theme]
primaryColor = "#4682B4"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
font = "sans serif"

[server]
maxUploadSize = 20
EOL
fi

# Set up GitHub Actions workflow
if [ -f ".github/workflows/deploy.yml" ]; then
    echo "GitHub Actions workflow already exists"
else
    echo "Creating GitHub Actions workflow"
    # Add workflow content
    cat > .github/workflows/deploy.yml << EOL
name: Deploy to Streamlit Cloud

on:
  push:
    branches:
      - main
  workflow_dispatch:

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout repository
        uses: actions/checkout@v2
      
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt
      
      - name: Test dependencies
        run: |
          python -c "import pandas, numpy, matplotlib, seaborn, sklearn, streamlit; print('All dependencies successfully imported')"
      
      - name: Check file structure
        run: |
          ls -la
          echo "✅ Repository is ready for Streamlit Cloud deployment"
EOL
fi

# Take a screenshot for the demo if needed
echo "Repository setup complete! You may want to add a screenshot of your app to the demo_images folder."
echo "Run your app locally with 'streamlit run pca-streamlit-app.py' and take a screenshot for the README."

echo "✅ Setup complete. Next steps:"
echo "1. Initialize git repository: git init"
echo "2. Add all files: git add ."
echo "3. Commit changes: git commit -m \"Initial commit\""
echo "4. Create repository on GitHub named 'pca-tool'"
echo "5. Push to GitHub: git remote add origin https://github.com/visvikbharti/pca-tool.git && git push -u origin main"