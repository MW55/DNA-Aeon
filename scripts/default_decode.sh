#!/bin/sh

# if -h flag is passed, print help
if [ "$1" = "-h" ]; then
    echo "Usage: $0 [input_file] [output_file]"
    echo "Decodes a file using the default decoder"
    exit 0
fi
#../python/decode.py -c $1 -m $2
