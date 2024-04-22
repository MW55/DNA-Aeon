#!/bin/sh

# if -h flag is passed, print help
if [ "$1" = "-h" ]; then
    echo "Usage: $0 [input_file] [output_file]"
    echo "Evaluates a file using the default evaluator"
    exit 0
fi
# ../python/encode.py -c $1 -m $2