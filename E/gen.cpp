#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <unordered_map>
#include <map>
#include <set>
#include <random>

using namespace std;

int main()
{
    int z = 1;
    cout << z << '\n';
    while (z--) {
        int n = 15;
        cout << n << '\n';
        mt19937 engine(time(0));
        uniform_int_distribution<int> xDist(0, 20);
        uniform_int_distribution<int> yDist(0, 20);
        for (int i = 0; i < n; ++i) {
            int x1 = xDist(engine);
            int x2 = xDist(engine);
            int y1 = xDist(engine);
            int y2 = xDist(engine);
            if (x1 == x2) {
                --i;
                continue;
            }
            cout << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << '\n';
        }
    }
}