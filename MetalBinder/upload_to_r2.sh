#!/bin/bash

# Configuration
BUCKET_NAME="metal-binder-data"  # Change this to your actual bucket name
DATA_DIR="/home/ryangustafson/Documents/GitHubProj/PhD-Research/Local/lanm_output"

# Check if Wrangler is installed
if ! command -v wrangler &> /dev/null
then
    echo "Wrangler could not be found. Please run 'npm install -g wrangler' first."
    exit 1
fi

echo "🚀 Starting upload of $DATA_DIR to R2 bucket: $BUCKET_NAME"

# Iterate over all files in the data directory
cd "$DATA_DIR" || exit 1

# Using 'find' to handle large numbers of files and subdirectories
find . -type f | while read -r file; do
    # Remove leading './' for the object key
    object_key="lanm_output/${file#./}"
    
    echo "Uploading: $object_key"
    
    # Run wrangler put with --remote flag to ensure it goes to the cloud
    wrangler r2 object put "$BUCKET_NAME/$object_key" --file "$file" --remote
    
    # Optional: add a small sleep if you hit rate limits
    # sleep 0.1
done

echo "✅ Upload complete!"
