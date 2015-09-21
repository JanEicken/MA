METHODS= init load save get douglas_peucker
COMMON_CXX = -g -Wall -Wextra -std=c++11 
OCT_CXX = `mkoctfile -p CXXFLAGS` $(COMMON_CXX)

all: links
	CXXFLAGS="$(OCT_CXX)" \
         mkoctfile trajcomp-octave.cpp 
links:
	@echo $(METHODS)
	for f in $(METHODS); do if ! [ -h trajcomp-$$f.oct ]; then ln trajcomp-octave.oct trajcomp-$$f.oct -s; fi; done
install:
	@echo "Installation for octave is not yet there"
	@echo "Please put all oct files (including symlinks) into"
	@echo "your octave project directory (or use addpath(...))"


clean:
	rm -f octave-core *.oct *.o
