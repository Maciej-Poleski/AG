#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <random>

using namespace std;

int main()
{
    int z = 1;
    cout << z << "\n";
    mt19937 engine(time(0));
    long limit = 1000000000;
    uniform_int_distribution<int> dist(-limit, limit);
    while (z--) {
        int n = 3500;
        long gen = 0;
        cout << n << "\n";
        for (int i = 0; i < n; ++i) {
            long x, y, r;
            do {
                x = dist(engine);
                y = dist(engine);
                r = x * x + y * y;
                gen += 1;
            } while (r > limit*limit || r < limit*limit/100*99);

            cout << x << " " << y << "\n";
        }
        //clog << "Wygenerowano " << gen << " punktów\n";
    }
    return 0;
}