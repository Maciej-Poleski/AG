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
        if (system("./gen killer")) {
            cout << "Generator zdechł\n";
            break;
        }
        int result = system("./solution < killer.in > killer.tout");
        if (result) {
            cout << "Wzorcówka zdechła\n";
            break;
        } else {
            cout << "OK?\n";
        }
        if (system("diff killer.out killer.tout")) {
            cout << "Błędna odpowiedź\n";
            break;
        }
    }
    return 0;
}