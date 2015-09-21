#ifndef BF_HPP_INC
#define BF_HPP_INC
#include<vector>
#include<iostream>
#include<sstream>
#include<string>

#include<math.h> 

#include "murmur.hpp"


using namespace std;


namespace trajcomp{
	namespace bloom {




template <class element_type> 
class BloomFilter{
public:
  std::vector<element_type> filter;
  int d;
  int gSeed;
  unsigned int numElements;
  
  BloomFilter()
  {
	  configure(0,0,-1);
  }
  
  BloomFilter(int d, int size, int seed =-1) 
  { 
	resize(size); 
	this->d = d;
	gSeed = seed;
	numElements = 0;
  }
  
  void configure(int d, int size, int seed=-1)
  {
	  resize(size); 
	this->d = d;
	gSeed = seed;
	numElements = 0;
  }
  
  
  
  
  
  void summary()
  {
	  cout << "Filter Summary:\t"<< endl;
	  cout << "Size: \t" << filter.size() << endl;
	  int ones = 0;
	  for(size_t i=0; i < filter.size(); i++)
		if (filter[i] == 1) 
			ones ++;
	  cout << "Ones: \t" << ones << endl;
  }
  void resize(int size)
  {
	  numElements = 0;
      filter.clear(); for (size_t i =0; i < size; i++) filter.push_back(0);
  }
  
  double foz()
  {
	  size_t o=0;
	  for (size_t i=0; i < filter.size(); i++)
	    if (filter[i] == 0)
	      o++;
	  return  (double) o / (double) filter.size();
  }

  inline double _esize()
  {
	  return -log(foz())*filter.size()/d;
  }

  unsigned int esize()
  {
	  double num = _esize();
	  return (unsigned int) num;
  }
  inline double _eUnion(BloomFilter &b)
  {
	  //sets. Swamidass & Baldi (2007) show
	  double dot = 0;
	  for (size_t i =0; i < filter.size(); i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  // dot 
	  return -log(1- dot / (double) filter.size())*filter.size()/d;
	  
  }
  unsigned int eUnion(BloomFilter &b)
  {
	    
	  return (unsigned int) _eUnion(b);
  }
  inline double _eIntersect(BloomFilter &b)
  {
	  double dot = 0;
	  for (size_t i =0; i < filter.size(); i++)
	    if (filter[i] == 1 || b.filter[i] == 1)
	      dot ++;
	  double union_size = _eUnion(b);
	  double e = _esize() + b._esize() - union_size;
	  
	  return e;
	  
  }
  unsigned int eIntersect (BloomFilter &b)
  {
	    
	  return (unsigned int) _eIntersect(b);
  }



friend std::ostream& operator<< (std::ostream& lhs, const BloomFilter & p)
{
   for (size_t i=0; i < p.filter.size(); i++)
     lhs << (unsigned int) p.filter[i] << " ";
   return lhs;
}

   void dump()
   {
	for (size_t i=0; i < filter.size(); i++) cout << (unsigned int) filter[i] <<" ";
   cout << endl;   
   }


unsigned int hash(int index,const char *s, char orient='+')
{
    std::string st(s);
    return hash(index,st, orient);
}

unsigned int hash(int index,std::string elem, char orient='+')
{
  std::stringstream stream;
  stream << index << elem << orient;

  std::vector<uint32_t> hash = 
	//trajcomp::murmur::murmur(stream.str(),this->gSeed);
	trajcomp::murmur::murmur(elem,index);
  
  return hash[1] % filter.size();
}

void add(std::string elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::add(" << code.size() <<"," << elem<< ");" << std::endl;
	#endif
    for (size_t i=0; i < d; i++)
    {
		filter[hash(i,elem)] = 1;
    }
    numElements ++;
}

bool contains(std::string elem)
{
	#ifdef LOG_OPS
	std::cout << "BF::contains(" << code.size() <<"," << elem<< ")==";
	#endif
	for (size_t i=0; i < d; i++)
    {
		if (filter[hash(i,elem)] ==0) 
	         return false;
    }
    return true;
}
		bool empty() 
		{
			for (size_t i=0; i < filter.size(); i++)
			  if (filter[i] != 0)
				return false;
			return true;
		}


}; // class BloomFilter


} // bloom
} // trajcomp
#endif
