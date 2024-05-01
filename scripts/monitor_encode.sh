# Define the path to the program you want to execute
PROGRAM_PATH="${HOME_DIR}/DNA-Aeon/python/encode.py"

# Get the current date and time in YYYY-MM-DD-HHMM format
TIMESTAMP=$(date +%Y-%m-%d-%H%M)

# Define where to store output and error files with timestamp
STDOUT_FILE="${HOME_DIR}/DNA-Aeon/logs/encode/${TIMESTAMP}/stdout_${TIMESTAMP}.log"
STDERR_FILE="${HOME_DIR}/DNA-Aeon/logs/encode/${TIMESTAMP}/stderr_${TIMESTAMP}.log"

# Create the logs directory if it does not exist
mkdir -p "${HOME_DIR}/DNA-Aeon/logs/encode/${TIMESTAMP}"

# Create the files if they do not exist
touch ${STDOUT_FILE}
touch ${STDERR_FILE}

# Execute the program and redirect stdout and stderr
echo "Scripts has to be run from "./DNA-Aeon""
echo " JSON file content : "
python3 /home/mathis/DNA-Aeon/scripts/json_print.py
echo "reminder of tmux"
echo "(done) tmux new -s session_name to create a new session"
echo "       tmux ls to list all sessions"
echo "       tmux attach -t session_name to attach to a session when detached"
# I want to create a tmux session to run the script
tmux has-session -t encode_session 2>/dev/null
if [ $? == 0 ]; then
    tmux kill-session -t encode_session
fi
# Start the session and send the command to execute your script, then terminate
tmux new-session -d -s encode_session "python3 ${PROGRAM_PATH} -c '/home/mathis/DNA-Aeon/configs/config-files/config.json' 1>${STDOUT_FILE} 2>${STDERR_FILE}; tmux kill-session -t encode_session"

# Check the return code of the program
if [ $? -eq 0 ]; then
        echo "script Succeed (returned 0) : Check logs to validate"
else
        echo "script Failed : Check logs to identify"
fi