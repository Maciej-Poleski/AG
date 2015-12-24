#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <unistd.h>
#include <fcntl.h>
#include <fstream>

namespace std {
    template<>
    struct hash<pair<int, int>>
    {
        size_t operator()(const pair<int, int> &p) const noexcept
        {
            return p.first | (static_cast<size_t>(p.second) << 32);
        }
    };
};

using namespace std;


int main(int argc, char **argv)
{
    if (argc != 2) {
        cerr << argv[0] << "[name]\n";
        return 1;
    }
    ofstream input(string() + argv[1] + ".in");
    ofstream output(string() + argv[1] + ".out");
    int z = 1;
    input << z << '\n';
    while (z--) {
        int n = 10;
        input << n << '\n';
        int xdim = 5;
        int ydim = 5;
        int queryMargin = 3;
        uniform_int_distribution<int> setK(-8, 12);
        uniform_int_distribution<int> setX(queryMargin, xdim + queryMargin);
        uniform_int_distribution<int> setY(queryMargin, ydim + queryMargin);
        uniform_int_distribution<int> queryX(0, xdim + 2 * queryMargin);
        uniform_int_distribution<int> queryY(0, ydim + 2 * queryMargin);
        uniform_int_distribution<int> command(0, 1);
        uniform_int_distribution<int> coordSelect(0, n - 1);
        mt19937 engine(time(0));
        vector<pair<int, int>> coords(n);
        unordered_set<pair<int, int>> base;
        vector<vector<int>> MAP(xdim + 2 * queryMargin + 1);
        for (int i = 0; i <= xdim + 2 * queryMargin; ++i) {
            MAP[i].resize(ydim + 2 * queryMargin + 1);
        }
        for (auto &p : coords) {
            do {
                p.first = setX(engine);
                p.second = setY(engine);
            } while (base.count(p));
            base.insert(p);
            input << p.first << " " << p.second << "\n";
        }
        int q = 6;
        input << q << '\n';
        while (q--) {
            if (command(engine) == 0) {
                // set
                int i = coordSelect(engine);
                int k = setK(engine);
                input << "SET " << coords[i].first << " " << coords[i].second << " " << k << '\n';
                MAP[coords[i].first][coords[i].second] = k;
            } else {
                // query
                int x1 = queryX(engine);
                int x2 = queryX(engine);
                int y1 = queryY(engine);
                int y2 = queryY(engine);
                if (x1 > x2) {
                    swap(x1, x2);
                }
                if (y1 > y2) {
                    swap(y1, y2);
                }
                input << "QUERY " << x1 << " " << y1 << " " << x2 << " " << y2 << "\n";
                int result = 0;
                for (int i = x1; i <= x2; ++i)
                    for (int j = y1; j <= y2; ++j)
                        result += MAP[i][j];
                output << result << '\n';
            }
        }
    }
}