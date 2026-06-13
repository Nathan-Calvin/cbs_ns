#!/bin/bash

echo "-------------------------------------------------"
echo "-------------------------------------------------"
echo "-------------------------------------------------"
echo "WARNING: THIS IS A DRY RUN. It will NOT run MESA!" 
echo "         This will generate tiny text files that act as placeholders"
echo "         for the real MESA outputs."
echo "The 'real' MESA rn file is presently called rn.bak."

mkdir -p photos
echo "Test run data here." > photos/x100
echo "Test run data here." > photos/x101
echo "Test run data here." > photos/x102
echo "Test run data here." > photos/x103

mkdir -p LOGS
echo "Test run file." > LOGS/history.data
echo "Test run file." > LOGS/profiles.index
echo "Test run file." > LOGS/profile1.data
echo "Test run file." > LOGS/profile2.data
echo "Test run file." > LOGS/profile3.data

echo "-------------------------------------------------"
echo "-------------------------------------------------"
echo "-------------------------------------------------"

