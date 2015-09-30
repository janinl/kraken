CXX = /illumina/thirdparty/gcc/gcc-4.9.2/bin/g++
CXXENHANCEDFLAGS = -Wall -w -I/usr/local/include -I/home/avoicu/localperl/lib/5.22.0/x86_64-linux/CORE -L/home/avoicu/localperl/lib/5.22.0/x86_64-linux/CORE -fopenmp -O3 -std=c++11
CXXFLAGS = -Wall -fopenmp -O3 -std=c++11
PROGS = db_sort set_lcas classify make_seqid_to_taxid_map db_shrink

.PHONY: all install clean
	
all: $(PROGS) 
	
install: $(PROGS)
	cp $(PROGS) /home/avoicu/github/kraken/

clean:
	rm -f $(PROGS) *.o

db_shrink: krakendb.o quickfile.o

db_sort: krakendb.o quickfile.o

set_lcas: krakendb.o quickfile.o krakenutil.o seqreader.o

classify: krakendb.o quickfile.o krakenutil.o seqreader.o
	 $(CXX) $(CXXENHANCEDFLAGS) classify.cpp krakendb.o quickfile.o krakenutil.o seqreader.o -o classify
	 

make_seqid_to_taxid_map: quickfile.o

krakenutil.o: krakenutil.cpp krakenutil.hpp
	$(CXX) $(CXXFLAGS) -c krakenutil.cpp 

krakendb.o: krakendb.cpp krakendb.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c krakendb.cpp 

seqreader.o: seqreader.cpp seqreader.hpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c seqreader.cpp  

quickfile.o: quickfile.cpp quickfile.hpp
	$(CXX) $(CXXFLAGS) -c quickfile.cpp 
