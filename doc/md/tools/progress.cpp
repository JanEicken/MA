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
