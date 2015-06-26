#!/bin/bash

#for i in {1..10000}

LIP=true
if [ `which qstat ` ]
then
    LIP=false
fi

while true ;
	do
         if [ ${LIP} ]
	 then
	     echo "running: " `qstat -u $1 | grep $1 | grep " r " | wc -l` "    total: " `qstat -u $1 | grep $1 | wc -l` 
	 else
	     echo "RUNNING: `bjobs | grep RUN | wc -l`            PENDING: `bjobs | grep PEND | wc -l`"
	 fi
	 sleep 5	
	 done
