#!/bin/bash

source ~/.bashrc

n_Evts=1
i=0

JOB=${1} 

if [ ! -d "./output" ]; then
    mkdir "output"
fi

if [ ! -d "./log" ]; then
    mkdir "log"
fi

while [ $i -le $n_Evts ]
do
    echo "ROOT: ./output/event${JOB}_${i}.root"
    echo "LOG: ./log/event${JOB}_${i}.log"
    ./sim "./output/event${JOB}_${i}.root" "./log/event${JOB}_${i}.root"
    echo "Event $i complete"
    echo " "
    let i=i+1
done