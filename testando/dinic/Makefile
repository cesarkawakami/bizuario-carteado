n=fastflow
debug:
	g++ -Wall -g -o $n $n.cpp
release:
	g++ -Wall -O2 -o $n $n.cpp
tests: debug
	./$n < in1 | diff -ysN - out1
