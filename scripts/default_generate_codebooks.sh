#!/bin/sh

# if -h flag is passed, print help
if [ "$1" = "-h" ]; then
    echo "Usage: ./default_generate_codebooks.sh"
    echo "Generate codebooks for the default configuration"
    echo "Default Path to Contrained Kaos"
    echo "GC is not defined, but should default to 40-60 as min-max"
    echo "HP is set to 4"
    echo "Motifs is not used"
    echo "Output is set to configs/codebooks/default_codebooks.json"
    echo "Length of output is set to ?"
    exit 0
fi
else
    # Generate codebooks for the default configuration
    ../python/generate_codebooks.py --Path "../libraries/ContrainedKaos/ConstrainedKaos.jar" --Output "../configs/codebooks/default_codebooks.json" -hp 4