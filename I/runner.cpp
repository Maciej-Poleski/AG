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
    for (; ;) {
        if (system("./gen > killer2.in")) {
            cout << "Generator zdechł\n";
            break;
        }
        int result = system("./solution < killer2.in > killer2.tout");
        if (result) {
            cout << "Wzorcówka zdechła\n";
            break;
        } else {
            cout << "OK?\n";
        }
    }
    return 0;
}