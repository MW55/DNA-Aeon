#!/bin/sh

# if -h flag is passed, print help
if [ "$1" = "-h" ]; then
    echo "Usage: $0 [input_filename]"
    echo "You have to respect the paths structure"
    echo "NOREC4DNA_BASE_PATH should be ./libraries/norec4dna"
    echo "DNA_AEON_PATH should be ./"
    echo "Filename should be the name of the file that should be in ./data/ folder"
    echo "Input filename should be the name of the file that should be in ./data/ folder"
    echo "CONFIGS_PATH should be ./configs/"
    echo "ENCODED_PATH should be ./data/encoded/"
    echo "Evaluates a file using the default evaluator"
    exit 0
fi
# ../python/error_simulation.py -c $1 -m $2