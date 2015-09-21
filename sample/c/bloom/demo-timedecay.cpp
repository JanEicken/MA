#include<iostream>
#include<stdexcept>
#include<string>

#include<trajcomp/trajcomp.hpp> 	// tools::make_string
#include<trajcomp/murmur.hpp>		// murmur
#include<trajcomp/bloomfilter.hpp>	// BloomFilter

#include <trajcomp/timedecayingbloomfilter.hpp>


using namespace std;
using namespace trajcomp::bloom;	// for BloomFilter
using namespace trajcomp::murmur;	// for murmur

int main(int argc, char **argv)
{
	// Create a bloom filter with 5 hash functions and 1000 binary slots
	TimeDecayingBloomFilter bf(5,1000);
	
	// Add even numbered strings
	for (size_t i=0; i <5; i++)
	{
		stringstream s("");
		int epoch = rand() % 15;
		s << "String" << i;
		bf.add(s.str(),epoch);
		cout << "Adding  "<< s.str() << " for " << epoch <<  " rounds." << endl;
	}
	int t=0;
	while(!bf.empty())
	{
		cout << "Showing after t=" << t << endl;
		
		
		//  query all and show the results
		for (size_t i=0; i <5; i++)
		{
			stringstream s("");
			s << "String" << i;
			cout << "\tTesting " << s.str() << " returned " << 
				(bf.contains(s.str())?"yes":"no") << endl;
		}
		bf.decay(1);t++;
	}
}
