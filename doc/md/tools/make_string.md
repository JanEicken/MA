makestring: A flexible vector to string converter
===================

In some programming languages, more or less all objects get instantly converted to strings and are ready for being output or shown in debuggers, on a console or somewhere else. For C++, however, no such concept exists. Instead, every class can be used in stream outputs, when it implements
the stream operator "<<". For classes from foreign libraries including STL vector classes, this is not feasible as you would have to extend them just for this functionality.

For vectors made of stream-compatible types (all basic types have a << operator implementation), we took another common way of using a string stream for the data elements of our type. As this is quite
tedious to write anywhere in your source code, where you want to output a vector, we provide a template
make_string in the namespace tools, that does exactly this.
ean Distance Formula")

Quick Reference
-------------------
	// Class template
		template<class datatype, char delim=' '>
		std::string make_string(datatype data)
	
	// Typical Use
	trajcomp::tools::make_string(v)


Sample Code
-------------
The following code sample shows how to use the make_string function template.

[make\_string.cpp](make_string.cpp)

	#include<trajcomp/trajcomp.hpp>
	#include<vector>
	#include<iostream>
	using namespace std;
	int main(void)
	{
	  std::vector<double> A={1,2,3,4,5};
	  cout << "make_string: " << trajcomp::tools::make_string(A) << endl;
	  cout << "make_string (csv): " << trajcomp::tools::make_string<std::vector<double>,','>(A) << endl;
	  return 0;
	}

This example first creates a vector A and then outputs it via cout using the default make_string.

Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample make_string.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
This sample outputs the following text string:

	make_string: 1 2 3 4 5 
	make_string (csv): 1,2,3,4,5,

Concepts and Custom Types
------------------------------

make_string works for any type that supports STL-style iterators. In detail, it accesses 

*	the type datatype::iterator to instantiate an iterator,
*	the functions begin() and end() returning iterators,
*	the != operation on iterators,
*	needs the dereferenced iterator (e.g., the underlying type) to be stream-compatible.