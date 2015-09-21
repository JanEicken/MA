#include<trajcomp/trajcomp.hpp>
#include<vector>

using namespace std;
using namespace trajcomp;

class Point
{
    public:
        double x;
        double y;

    Point(double thex,double they):x(thex),y(they)
    {
    }
    // provide some interface for vector-style access

    const size_t size() { return 2; };
    double operator[] (size_t where)
    {
        return (where == 0)?x:y;
    }


    // provide a double-valued substraction
    double operator- (Point &p)
    {
        return sqrt( (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y));
    }

    friend std::ostream & operator<<(std::ostream &os, const Point& p)
    {
        os << p.x << " " << p.y;
        return os;
    }

};


int main() {
    Point p(1,1);
    Point q(3,3);

    trajcomp::default_element_distance<Point> d;

    cout << "The distance between "
         << p
         << " and "
         << q   
         << " is "
         << d(p,q) << endl; 

    return 0;
}

