#!/usr/bin/bash

# Define the path to the program you want to execute
PROGRAM_PATH="~/DNA-AEON/python/decode.py"

# Get the current date and time in YYYY-MM-DD-HHMM format
TIMESTAMP=$(date +%Y-%m-%d-%H%M)

# Define where to store output and error files with timestamp
STDOUT_FILE="~/DNA-AEON/logs/stdout_$TIMESTAMP.log"
STDERR_FILE="~/DNA-AEON/logs/stderr_$TIMESTAMP.log"

# Execute the program and redirect stdout and stderr
python3 ${PROGRAM_PATH} -c configs/config-files/config.json 1>${STDOUT_FILE} 2>${STDERR_FILE}
