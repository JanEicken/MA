tictoc: A simple timer object
===================

In order to analyze the runtime performance of your own programs, you can use as many instances of the class `tictoc` as you like in your project. Each of these classes stores the current time (as returned by the `ticks()` method, which can be implemented using `GetTickCount()` on Microsoft Windows.

Quick Reference
-------------------
	uint64_t ticks(void)
	{	
		struct timeval tv;
		gettimeofday(&tv, 0);
		return uint64_t( tv.tv_sec ) * 1000 + tv.tv_usec / 1000;
	}
	class tictoc
	{
		public:
		uint64_t t;
		void tic() {t = ticks();};
		uint64_t toc(){return ticks() - t; };
	};	// Function
	
	// Typical Use
	tictoc c;
	c.tic(); do_hard_work(); cout << "Hard Work took: " << c.toc() << endl;


Sample Code
-------------
The following code sample shows how to use the tictoc class:

[tictoc.cpp](tictoc.cpp)

	#include<trajcomp/trajcomp.hpp>
	#include<iostream>
	#include<string>
	#include<unistd.h> // for usleep
	using namespace std;
	using namespace trajcomp::tools;

	int main(void)
	{
	  tictoc c;
	  for (size_t i=0; i <= 100; i++)
	  {
		c.tic();
		usleep (1E5);
		cout << "usleep(1E5) took " << c.toc() << "ms" << endl;
	  }
	  
	  return 0;
	}

This example uses usleep to wait a given amount of time and tictoc to measure the overall time taken for that operation. 

Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample tictoc.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
This sample outputs the following:

	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 101ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms
	usleep(1E5) took 100ms


