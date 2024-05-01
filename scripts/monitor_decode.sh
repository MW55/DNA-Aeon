# Define the home dir
# We need two options to specify the mode of the script
# -m : mode of the script
# -c : hardcoded
# -i : to skip inner decode
# -o : to skip outer decode

if [ "$1" == "-i" ]; then
        MODE=true
else :
        MODE=false
fi

HOME_DIR="/Users/mguyot/Documents"

# Define the path to the program you want to execute
PROGRAM_PATH="${HOME_DIR}/DNA-Aeon/python/decode.py"

# Get the current date and time in YYYY-MM-DD-HHMM format
#TIMESTAMP=$(date +%Y-%m-%d-%H%M)

# Define where to store output and error files with timestamp
STDOUT_FILE="${HOME_DIR}/DNA-Aeon/logs/decode/debug/stdout_debug.log"
STDERR_FILE="${HOME_DIR}/DNA-Aeon/logs/decode/debug/stderr_debug.log"
#STDOUT_FILE="${HOME_DIR}/DNA-Aeon/logs/decode/${TIMESTAMP}/stdout_${TIMESTAMP}.log"
#STDERR_FILE="${HOME_DIR}/DNA-Aeon/logs/decode/${TIMESTAMP}/stderr_${TIMESTAMP}.log"

# Create the logs directory if it does not exist
#mkdir -p "${HOME_DIR}/DNA-Aeon/logs/decode/${TIMESTAMP}"

# Create the files if they do not exist
#touch ${STDOUT_FILE}
#touch ${STDERR_FILE}

#python3 /home/mathis/DNA-Aeon/scripts/json_print.py
# Execute the program and redirect stdout and stderr

# to-do : find a way to close my session w
#echo "this script use tmux session so you can log out and recover the session"
#echo "reminder of the control. you can detach from the session with Ctrl-B D"
#echo "to reattach : tmux attach -t DNA_Aeon"
#tmux new-session -d -s DNA_Aeon

if [ "$MODE" == "true" ]; then
    echo "Skipping inner decode"
    python3 ${PROGRAM_PATH} -m "sys"  -c "${HOME_DIR}/DNA-Aeon/configs/config-files/config.json" -i 1>${STDOUT_FILE} 2>${STDERR_FILE}
else
    echo "full decode"
    python3 ${PROGRAM_PATH} -m "sys"  -c "${HOME_DIR}/DNA-Aeon/configs/config-files/config.json" 1>${STDOUT_FILE} 2>${STDERR_FILE}
fi

if [ $? -eq 0 ]; then
        echo "script Succeed (returned 0) : Check logs to validate"
else
        echo "script Failed : Check logs to identify"
fi