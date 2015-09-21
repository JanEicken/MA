#include<trajcomp/trajcomp.hpp>
#include<vector>

using namespace std;
using namespace trajcomp;

typedef trajcomp::trajectory<double> traj;


int main() 
{
	traj t(2);
	t.push_back({1.0d,1.0d});
	t.push_back({2.0d,2.0d});
	t.push_back({3.0d,3.0d});
	t.push_back({4.0d,4.0d});
	t.push_back({5.0d,5.0d});
	t.push_back({6.0d,6.0d});
	t.push_back({7.0d,7.0d});
	t.push_back({8.0d,8.0d});
	cout << "Reducing from:" << endl;
	t.dump();
	
	traj q(2);
	cout << "Every second and the last:" << endl;
	q = trajcomp::uniform_select(t,2);
	q.dump();
	cout << "Every third and the last:" << endl;
	q = trajcomp::uniform_select(t,3);
	q.dump();
	cout << "Every third:" << endl;
	q = trajcomp::uniform_select(t,3,false);
	q.dump();
	
    
    return 0;
}
