all:
	@if ! [ -d /usr/include/trajcomp ]; then echo "Did not find includes, run make install first (from the main directory of this distribution).";exit 2; fi
	g++ -std=c++11 -o bloom bloom-sample.cpp 
	g++ -std=c++11 -o counting demo-counting.cpp
	g++ -std=c++11 -o timedecay demo-timedecay.cpp
clean:
	rm -f bloom counting timedecay
