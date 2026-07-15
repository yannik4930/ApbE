#!/bin/bash

# zählt Ordner im aktuellen Verzeichnis, deren Name mit "subdirectory" beginnt
ntasks=$(find . -maxdepth 1 -type d -name "Subdirectory*" | wc -l | tr -d ' ')

echo "Found $ntasks subdirectories"

if [ "$ntasks" -eq 0 ]; then
    echo "ERROR: No subdirectory* folders found."
    exit 1
fi

sbatch --ntasks="$ntasks" fresco_job.sh
