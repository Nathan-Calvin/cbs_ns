#!/usr/bin/env bash 
# 
# This script runs a MESAModel on a given host. 
#
# It takes one optional command-line argument, which is the nickname of this host.
# If that argument is missing, it uses the output of `hostname` instead.
#
# Author: Jess Vriesema 
# Date: 3 June 2026
# ------------------------------------------------------ 

echo "Running Bash script \"${BASH_SOURCE[0]##*/}\"..."

# Usage statement function
usage() {
	cat << EOF
	Usage: $(basename "$0") [hostname]

	The optional [hostname] argument specifies the name of this
	computer for when it outputs to the logfile. If not present,
	the hostname is determined from the output of the hostname
	command.
EOF
	exit 1
}


# Get the hostname
if [ $# -gt 1 ]; then
	echo "Too many input arguments."
	usage
	exit 1
elif [ $# -eq 1 ]; then
	HOSTNAME=$1
else
	HOSTNAME=$(hostname)
fi


LOGFILE="model_log.txt"
# The INIT_MESA script should set MESA_DIR and MESASDK_ROOT
# as well as source $MESASDK_ROOT/bin/mesasdk_init.sh.
INIT_MESA="/storage/wumas/init_mesa.sh"
 
# Output headers 
echo "============================================================" >> "$LOGFILE" 
echo "$HOSTNAME:  Running $0..." >> "$LOGFILE"
echo "============================================================" >> "$LOGFILE" 
# Set environment variables
# Does the INIT_MESA script exist? If not, tell user now:
if [ -f "$INIT_MESA" ]; then
	#echo "$HOSTNAME:  The INIT_MESA script '$INIT_MESA' exists." >> "$LOGFILE"
	#echo "$HOSTNAME:  Setting environment variables for MESA by running '$INIT_MESA' in the current shell..." >> "$LOGFILE"
	# Try running the script
	source "$INIT_MESA" >> "$LOGFILE"
	#echo "? == $?" >> "$LOGFILE"
	if [ $? -eq 0 ]; then
		: #echo "$HOSTNAME:     ...successfully ran '$INIT_MESA'!" >> "$LOGFILE"
	else
		echo "$HOSTNAME:     ...FAILED when running '$INIT_MESA'!" >> "$LOGFILE"
		exit 1
	fi
else
	echo "$HOSTNAME:  ERROR: The INIT_MESA script file '$INIT_MESA' could not be found." >> "$LOGFILE"
	exit 1
fi
#echo "$HOSTNAME:  TESTING MESA ENVIRONMENT VARIABLES:" >> "$LOGFILE"
#echo "$HOSTNAME:     MESA_DIR == '$MESA_DIR'" >> "$LOGFILE"
#echo "$HOSTNAME:     MESASDK_ROOT == '$MESASDK_ROOT'" >> "$LOGFILE"

echo "============================================================" >> "$LOGFILE" 
date >> "$LOGFILE"   # print the date when this run started 
pwd >> "$LOGFILE"    # print the current working directory (for future debugging?) 
tree >> "$LOGFILE"     # output the contents of the current directory 
echo "============================================================" >> "$LOGFILE" 


# Compile the model code
echo "$HOSTNAME:  Compiling the MESA model..." >> "$LOGFILE"
time ./mk >> "$LOGFILE" 2>&1
if [ $? -eq 0 ]; then
	echo "$HOSTNAME:     ...successfully compiled the MESA model!" >> "$LOGFILE"
else
	echo "$HOSTNAME:     ...FAILED when compiling the MESA model!" >> "$LOGFILE"
fi
echo "============================================================" >> "$LOGFILE" 


echo "$HOSTNAME:  Starting a MESA model run..." >> "$LOGFILE"
# Are we using a restart file?
if [ "$#" -ne 1 ]; then
	# No restart file
	# Run AutoMESA (the time command measures and outputs runtime)
	time ./rn >> "$LOGFILE" 2>&1
	if [ $? -eq 0 ]; then
		echo "$HOSTNAME:     ...successfully ran the MESA model!"
	else
		echo "$HOSTNAME:     ...FAILED when running the MESA model!"
	fi
else
	RESTARTFILE=$1
	# Does the restart file exist?
	if [ ! -f "$RESTARTFILE" ]; then
		echo "$HOSTNAME:  The restart file \"$RESTARTFILE\" does not exist on $(hostname)!" |& tee -a "$LOGFILE"
		exit 1
	else
		# Restart this MESA model from the specified RESTARTFILE
		# (The time command measures and outputs runtime.)
		time ./re $RESTARTFILE >> "$LOGFILE" 2>&1
		if [ $? -eq 0 ]; then
			echo "$HOSTNAME:     ...successfully ran a restarted MESA model!"
		else
			echo "$HOSTNAME:     ...FAILED when running a restarted MESA model!"
		fi
	fi
fi


# Output a footer 
echo "============================================================" >> "$LOGFILE" 
date >> "$LOGFILE"   # print the date when this run started 
pwd >> "$LOGFILE"    # print the current working directory (for future debugging?) 
tree >> "$LOGFILE"     # output the contents of the current directory 
echo "============================================================" >> "$LOGFILE" 

echo "$HOSTNAME:  Finished running $0."

