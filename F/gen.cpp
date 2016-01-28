#include <iostream>
#include <random>
#include <vector>

using namespace std;

int main()
{
    int z = 1;
    cout << z << '\n';
    while (z--) {
        const int n = 2000000; // parzyste
        const int side = n / 2;
        const int dimension = 500000000;
        const int spacing = dimension * 2 / side;
        const int top = dimension;
        vector<int> ys;
        uniform_int_distribution<int> left(-dimension, -100);
        uniform_int_distribution<int> right(100, dimension);
        uniform_int_distribution<int> tilt(-5, 5);
        mt19937 engine(time(0));
        for (int i = 0; i < side; ++i)
            ys.push_back(top - i * spacing + tilt(engine));
        cout << n << '\n';
        for (auto i = ys.rbegin(), e = ys.rend(); i != e; ++i) {
            cout << right(engine) << ' ' << *i << '\n';
        }
        for (auto y : ys)
            cout << left(engine) << ' ' << y << '\n';
    }
}