progress: A simple progress bar for Linux terminals
===================

When writing data-centric software, it is often useful to be able to watch complex algorithms progress. Otherwise, it is difficult to estimate, when a given program will be ready. For this purpose, libtrajcomp provides a small template function using escape sequences for the linux terminal in order to overwrite a progress text. Currently, the progress is expected to be linear, other runtime-behavious might be added later, when we need them.

Quick Reference
-------------------
	// Function
	void progress(long pos, long size,std::string what="")
		
	
	// Typical Use
	progress(20,100,"Hard Work")


Sample Code
-------------
The following code sample shows how to use the progress function.

[progress.cpp](progress.cpp)

	#include<trajcomp/trajcomp.hpp>
	#include<string>
	#include<unistd.h> // for usleep
	using namespace std;
	using namespace trajcomp::tools;

	int main(void)
	{	
	  for (size_t i=0; i <= 100; i++)
	  {
		progress(i,100,"Sleeping");
		fflush(stdout);
		usleep (1E5);
	  }
	  printf("\n");
	  return 0;
	}

This example uses usleep to wait a given amount of time. Two decisions, we made, should be
mentioned:

*	progress does not flush the output buffer, hence, without calling 'fflush(stdout)' the 
	terminal gets updated only from time to time. It is definitely no good idea to flush every 
	call to progress, when the waiting time between progress steps is too short
*	progress overwrites the current line without adding a newline character. Therefore, after
	using progress, you can / should start manually on a new line. In this way, different progress 
	bars can overlay each other in the same output line and make your program less verbose.


Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample progress.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
This sample outputs the following text string changing over time:

	Progress(Sleeping): 35/100(99%)

