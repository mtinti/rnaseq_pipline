#!/bin/bash

N_PROC=$1

LOG=$(mktemp jupyter.log.XXXXX)
NAME="nb_${LOG##*.}"

if [ -z $N_PROC ]; then
    qrsh -cwd -V -N $NAME \
        jupyter notebook \
	--no-browser \
        --ip=\$\(hostname --fqdn\) \
	&>${LOG} &disown;
        
else
    qrsh -cwd -V -N $NAME -pe smp $N_PROC \
        jupyter notebook \
	--no-browser \
        --ip=\$\(hostname --fqdn\) \
	&>${LOG} &disown;
fi

PID=$!

echo "Server started with PID ${PID}"
echo "Logging to ${LOG}"





