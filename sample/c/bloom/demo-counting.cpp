#include<iostream>
#include<stdexcept>
#include<string>

#include<trajcomp/trajcomp.hpp> 	// tools::make_string
#include<trajcomp/murmur.hpp>		// murmur
#include<trajcomp/bloomfilter.hpp>	// BloomFilter

#include <trajcomp/countingbloomfilter.hpp>


using namespace std;
using namespace trajcomp::bloom;	// for BloomFilter
using namespace trajcomp::murmur;	// for murmur

int main(int argc, char **argv)
{
	// Create a bloom filter with 5 hash functions and 1000 binary slots
	CountingBloomFilter bf(5,1000);
	
	// Add even numbered strings
	for (size_t i=0; i <5; i++)
	{
		stringstream s("");
		int multiplicity = rand() % 15;
		s << "String" << i;
		bf.add(s.str(),multiplicity);
		cout << "Adding  "<< s.str() << " with multiplicity " << multiplicity << endl;
	}
	int t=0;
	while(!bf.empty())
	{
		// Remove random element
		stringstream s("");
		s << "String" << rand() % 5;
		cout << "\tCounting " << s.str() << " returned " << bf.count(s.str()) << endl;
		size_t mul = rand() % 15;
		cout << "\tRemoving " << s.str() << " " << mul << " times was possible? " << 
			(bf.remove(s.str(),mul)?"yes":"no") << endl;
		cout << "\t Counting again: " << bf.count(s.str()) << endl;
		
		//  query all and show the results
		for (size_t i=0; i <5; i++)
		{
			stringstream s("");
			s << "String" << i;
			cout << "\t" << s.str() << ": " << 
				(bf.contains(s.str())?"yes":"no") << 
				"[" << bf.count(s.str()) << "]" << endl;
		}
		t++;
	}
}
