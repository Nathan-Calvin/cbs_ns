#!/bin/bash

# Simple script that cleans up the directories on all of the remote systems for testing this
# program

# made by Ryan Deaton


echo "Cleaning all remote and local directories!"

ssh ryan@153.106.201.40 "./server_del_script.sh"
echo "All remote directories cleared on Xeon Phi"

rm -r ~/AutoMESA_Storage/*
rm -r /storage/wumas/AutoMESARuns/*
echo "All Local Directories Cleared"

ssh rad53@borg.calvin.edu "rm -r /storage/wumas/AutoMESA/*"
echo "All Borg Directories Cleared."
