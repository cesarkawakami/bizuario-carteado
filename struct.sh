#!/bin/bash
for str in "$@"
do
    mkdir -p $str
    echo n=$str > $str/Makefile
    cat Makefile >> $str/Makefile
    cp modelo.cpp $str/$str.cpp
done
