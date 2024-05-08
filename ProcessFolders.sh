#!/bin/bash

# function to sum all files
function sumMatchedFiles {
    mkdir Sum 2>/dev/null
    for i in {0..7}; do
        echo "Processing Board $i"
        ls Source/*/Processed/Board$i-Aligned.root
        hadd -f Sum/Board$i-Aligned.root Source/*/Processed/Board$i-Aligned.root
    done
}

# Function of Generate Calibration file
function generateCaliFile {
    echo "Processing Calibration"
    root -l -b -q "Codes/Calibration.cpp(\"Cali\",\"Sum/\")"
}

# Function of sum all GetPos.root files
function sumPosFiles {
    mkdir Sum 2>/dev/null
    echo "Processing Board $i"
    ls Source/*/Processed/GetPos.root
    hadd -f Sum/GetPos.root Source/*/Processed/GetPos.root
}

# # function to link file to directory
# function linkCaliFileToDir {
#     file=$1
#     dir=$2
#     ln -s $file $dir/$(basename $file)
# }

# for dir in $(ls -d Source/*); do
#     echo "Processing $dir"
#     cd $dir
#     for i in {0..7}; do
#         echo "Linking Board $i to " $dir
#         linkCaliFileToDir ../Cali/Board$i-Cali.root $dir/Processed/Cali
#     done
#     cd ..
# done

# Function of getposition through root executable script GetPos.cpp
function getPosition {
    # Iterator
    Iterator=0

    for dir in $(ls -d Source/*); do
        # Judge whether number of threads is larger than 4
        # If so, wait until the number of threads is less than 4

        echo "Processing $dir"
        root -l -b -q "Codes/GetPos.cpp(\"Cali/\",\"$dir/Processed/\")" >/dev/null 2>&1 &
        # echo $(ps -ef | grep "root -l -b -q Codes/GetPos.cpp(\"Cali/\",\"$dir/Processed/\")")
        Iterator=$((Iterator + 1))
        # if Iterator mod 8 == 0, wait until all threads are finished
        if [ $((Iterator % 8)) -eq 0 ]; then
            # while [ $(ps -ef | grep "root -l -b -q Codes/GetPos.cpp(\"Cali/\",\"$dir/Processed/\")" | wc -l) -gt 4 ]; do
            wait
            # done
        fi
    done
}

# Function of correct position from GetPos.root
function correctPosition {
    echo "Processing correction"
    root -l -b -q "Codes/Correct1.cpp(\"Sum/\")"
}

# Function of extracting imaging data
function extractImagingData {
    echo "Processing extraction"
    root -l -b -q "Codes/ExtractToImaging.cpp(\"Sum/\")"
}

# sumMatchedFiles #2> /dev/null
# generateCaliFile   #2> /dev/null
getPosition        #2> /dev/null
sumPosFiles        #2> /dev/null
correctPosition    #2> /dev/null
extractImagingData #2> /dev/null
# root -l -b -q "Codes/GetPos.cpp(\"Cali/\",\"Source/2024-04-01-23-32-31/Processed/\")"
