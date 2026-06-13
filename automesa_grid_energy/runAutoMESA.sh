#!/usr/bin/env bash 
# 
# This script runs a MESAModel on a given host. 
# Author: Jess Vriesema 
# Date: 10 July 2025
# ------------------------------------------------------ 

AUTOMESA=/storage/wumas/CntrlAutoMESA/AutoMESA.py
RUNFILE=$1

LOGFILE="run_log.txt" 
# For testing purposes only:
#LOGFILE="/dev/null"


# Output a header 
echo "===========================================" |& tee -a $LOGFILE 
date >> $LOGFILE   # print the date when this run started 
pwd >> $LOGFILE    # print the current working directory (for future debugging?) 
tree >> $LOGFILE     # output the contents of the current directory 
echo "===========================================" |& tee -a $LOGFILE 
 
# Does the restart file exist?
if [ ! -f "$RUNFILE" ]; then
	echo "The file \"$RUNFILE\" does not exist!" |& tee -a $LOGFILE
	exit 1
fi

# Run AutoMESA (the time command measures and outputs runtime)
echo "$0: About to run AutoMESA...." |& tee -a $LOGFILE
time python3 $AUTOMESA $RUNFILE |& tee -a $LOGFILE 
echo "$0: Finished running AutoMESA." |& tee -a $LOGFILE
 
# Output a footer 
echo "===========================================" |& tee -a $LOGFILE 
date >> $LOGFILE   # print the date when this run started 
pwd >> $LOGFILE    # print the current working directory (for future debugging?) 
tree >> $LOGFILE     # output the contents of the current directory 
echo "===========================================" |& tee -a $LOGFILE 

