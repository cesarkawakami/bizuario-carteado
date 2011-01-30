for i in a b c; do mkdir $i; echo n=$i > $i/Makefile; cat Makefile >> $i/Makefile;
cp modelo.cpp $i/$i.cpp; done
