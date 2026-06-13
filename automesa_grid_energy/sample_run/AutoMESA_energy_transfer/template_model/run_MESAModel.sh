#!/usr/bin/env bash 
# 
# This script runs a MESAModel on a given host. 
# Author: Jess Vriesema 
# Date: 10 July 2025
# ------------------------------------------------------ 


#Needed parameters for BORG hosts. leave commented out.
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -J AUTOMESAII

LOGFILE="model_log.txt" 
 
# Output a header 
echo "============================================================" >> $LOGFILE 
echo "Running $0 on $(hostname)...."
date >> $LOGFILE   # print the date when this run started 
pwd >> $LOGFILE    # print the current working directory (for future debugging?) 
tree >> $LOGFILE     # output the contents of the current directory 
echo "============================================================" >> $LOGFILE 


echo "============================================================" >> $LOGFILE 
echo "Initializing the enviornment variables."
MESA_DIRECTORY=$(grep "export MESA_DIR" ~/.bashrc)
MESASDK_DIRECTORY=$(grep "export MESASDK_ROOT" ~/.bashrc)
THREADS=$(grep "export OMP_NUM_THREADS" ~/.bashrc)


cmd1=$MESA_DIRECTORY
cmd2=$MESASDK_DIRECTORY
cmd3=$THREADS

eval $MESA_DIRECTORY
eval $MESASDK_DIRECTORY
eval $THREADS

echo $MESA_DIR, $MESA_DIRECTORY >> $LOGFILE
echo $MESASDK_ROOT, $MESASDK_DIRECTORY >> $LOGFILE
echo $OMP_NUM_THREADS, $THREADS >> $LOGFILE

MESA_SDK_INIT="source $MESASDK_ROOT/bin/mesasdk_init.sh"

eval $MESA_SDK_INIT

echo $MESA_SDK_INIT >> $LOGFILE

echo "============================================================" >> $LOGFILE 

# Compile the model code
echo "Compiling the MESA model..." >> $LOGFILE
time ./mk >> $LOGFILE 2>&1
echo "============================================================" >> $LOGFILE 


echo "Starting a MESA model run..." >> $LOGFILE
# Are we using a restart file?
if [ "$#" -ne 1 ]; then
	# No restart file
	# Run AutoMESA (the time command measures and outputs runtime)
	time ./rn >> $LOGFILE 2>&1
else
	RESTARTFILE=$1
	# Does the restart file exist?
	if [ ! -f "$RESTARTFILE" ]; then
		echo "The restart file \"$RESTARTFILE\" does not exist on $(hostname)!" |& tee -a $LOGFILE
		exit 1
	else
		# Restart this MESA model from the specified RESTARTFILE
		# (The time command measures and outputs runtime.)
		time ./re $RESTARTFILE >> $LOGFILE 2>&1
	fi
fi


# Output a footer 
echo "============================================================" >> $LOGFILE 
date >> $LOGFILE   # print the date when this run started 
pwd >> $LOGFILE    # print the current working directory (for future debugging?) 
tree >> $LOGFILE     # output the contents of the current directory 
echo "============================================================" >> $LOGFILE 

