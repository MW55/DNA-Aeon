# Notes 17.04.24

## eval.py program
- [ ] Goal of the function
    calling arithmetic_modulator_error_correction
        - with '-d config_file_path' (eval-config.json)
    res = run_eval (True or False)
    
## list of ideas to do
- [ ] put order by grouping things
    - group the python files
    - [ ] put the cmake files in the build no in ./
    - [ ] generate the binaries in the bin folder
    - [ ] group the libraries under the libraries folder
    - [ ] create a docker folder 

----

# Update29.04.24
- [x] validation using the Dorn file (without error) (4.8kb) => OK
- [x] validation using a JPEGXL (.jxl) file (150 kb) => error
- [x] validation using a .jpeg file (50 kb) => error ?

## Actions ?
- observe interaction between encode and decode


---
# Update 1.05.24 :
- [x] monitor_encode.sh and monitor_decode.sh : scripts to run on a other computer to collect the outputs and and assert program
## Problems
- encoding el.jpg seems to "works"
- but not decoding (or it takes a lot of time)
## Observations
- the fsanitize showed there's a stack depth crash (meaning too many recursions)
    - the ECdecoding::checkpointCheck is recursive and seemes to never ends
    - so the error somewhat comes from the CRC generations
    - so we have like invalid .zip for the NOREC4DNA 
    - there's a warning in stdout from /DNA-Aeon/logs/decode/debug/stdout_debug.log
        - "Warning Config-entry 'checksum' not known!"
        - No Packet was correctly decoded. Check your configuration.
### decode.py
- decode.py calls 
    - ConfigWorker.py in lib folder
        - the config file is the .ini file created
        - we then create a class and call.execute() on it
            -  
    
```markdown
[/Users/mguyot/Documents/DNA-Aeon/data/decoded.txt.zip]
algorithm = RU10
error_correction = crc
insert_header = True
savenumberofchunks = False
mode_1_bmp = False
upper_bound = 0.5
number_of_chunks = 12
config_str = USE_HEADER_CHUNK: True, NUMBER_OF_CHUNKS: 12 NUMBER_OF_CHUNKS_LEN_FORMAT: I ID_LEN_FORMAT: H ERROR_CORRECTION: crc32 CRC_LEN_FORMAT(Optional): I FILE: data/00005_560x888_69.jxl OUT_FILE:  Distribution: RaptorDistribution_S=12
id_len_format = H
number_of_chunks_len_format = I
packet_len_format = I
crc_len_format = I
master_seed = 2147483648
distribution = RaptorDistribution_S=12
rules = []
chunk_size = 3618
dropped_packets = 0
created_packets = 17
checksum = 9 
checksum_len_str = B
repair_symbols = 2
asdna = False
number_of_splits = 1
read_all = True
```

### Solve
- 1) it's written ERROR_CORRECTION being crc32 (not sure) => I stand for 32
the line checksum gives a warning about checksum being not used ?
in our case number of split if not equal to 0
- 2) 
- 3) the .zip seems empty

### Results of test1
- rewritting of the recursive call because of stack track pb
- it seems that we reached a max queue size (3/4)
    - how is it handled ?
    - programs doesn't seems to stop 

#### Analysis of ECDecoding.cpp
- so in queueCheck
    - ...
    - at the end
        - we check is queue.size is bigger than a config number (200'000 in our case)
        - we increase queueCounter++
            - if queueCounter is > than a config number (1)
                we throw an error
            - else 
                we like advance and iter
                (meaning we cut in half the queue size)


---

# Code Analysis of the CRC
- for each thread (4), we call do decode
do_decode consists of :
1) creating a ECdecoding object
2) creating a SeqEntry by method .decode()
3) updating metric_str
4) get data
then we store data in the results array
results is a list<tuple<string,vector<unsigned char>>>
result.emplace_back(metric_str, data) // a tuple
    metric_str is a double !
    get_data returns a out.data which store an array of symbols

## details of step 2 : calling decode
### preparation steps
-

### try main loop
- 
### check CRC
- setup with queue 
- SeqEntry 
    we call res.ac.finish()
    we checkCandidate()


---

# 08052024
- how to control the length of oligos to be in the 200-300 range rather than 2000
- impact of changing the packet redundancy 
- impact of chunk size ?
