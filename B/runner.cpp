#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <random>

using namespace std;

int main()
{
    for (; ;) {
        system("./gen > killer.in");
        int result = system("./solution < killer.in");
        if (result) {
            cout << "Uśmiercony\n";
            break;
        }
    }
    return 0;
}