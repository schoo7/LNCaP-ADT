#!/bin/bash

# Define source and destination directories
SOURCE_DIR="/media/simon/external1/Lin_RRBS/usftp21.novogene.com/01.RawData"
DEST_DIR="/media/simon/external1/RRBS_raw"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find and copy all .fq.gz files from the source directory to the destination directory
find "$SOURCE_DIR" -type f -name "*.fq.gz" -exec cp {} "$DEST_DIR" \;

echo "All .fq.gz files have been copied to $DEST_DIR"
