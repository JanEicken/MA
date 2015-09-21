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
