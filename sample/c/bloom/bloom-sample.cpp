#include<iostream>
#include<stdexcept>
#include<string>

#include<trajcomp/trajcomp.hpp> 	// tools::make_string
#include<trajcomp/murmur.hpp>		// murmur
#include<trajcomp/bloomfilter.hpp>	// BloomFilter

using namespace std;
using namespace trajcomp::bloom;	// for BloomFilter
using namespace trajcomp::murmur;	// for murmur

void main_bloom(void)
{
	// Create a bloom filter with 5 hash functions and 1000 binary slots
	BloomFilter<bool> bf(5,1000);
	// Add even numbered strings
	for (size_t i=0; i <= 10; i+=2)
	{
		stringstream s("");
		s << "String" << i;
		bf.add(s.str());
		cout << "Just put " << s.str() << " into the filter" << endl;
	}
	// Now query all and show the results
	for (size_t i=0; i <= 10; i++)
	{
		stringstream s("");
		s << "String" << i;
		cout << "Testing " << s.str() << " returned " << 
		(bf.contains(s.str())?"yes":"no") << endl;
	}
		
}


void main_murmur(std::string s)
{
	cout << "Calculating Murmur Hash for " << s << endl;
	std::vector<uint32_t> hash = murmur(s);
	std::string str = trajcomp::tools::make_string(hash);
	cout << str << endl;
}


int main(int argc, char **argv)
{
	switch(argc)
	{
		case 2:
			main_murmur(argv[1]);
			break;
		default:
			main_bloom();
		
	}
	
}
