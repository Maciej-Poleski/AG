#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>

using namespace std;

struct Rectangle
{
    int xLow, yLow, xHigh, yHigh;
};

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        cin >> n;
        int xLow, xHigh, yLow, yHigh;
        xLow = yLow = numeric_limits<int>::max();
        xHigh = yHigh = numeric_limits<int>::min();
        vector<Rectangle> input(n);
        for (auto &rect : input) {
            cin >> rect.xLow >> rect.yLow >> rect.xHigh >> rect.yHigh;
        }
        for (auto rect : input) {
            xLow = min(xLow, rect.xLow);
            yLow = min(yLow, rect.yLow);
            xHigh = max(xHigh, rect.xHigh);
            yHigh = max(yHigh, rect.yHigh);
        }
        assert(xLow < xHigh);
        assert(yLow < yHigh);
        const unsigned width = xHigh - xLow;
        const unsigned height = yHigh - yLow;
        vector<vector<bool>> matrix(height);
        for (auto &row : matrix) {
            row.resize(width, false);
        }
        for (auto rect : input) {
            for (int i = rect.yLow - yLow; i < rect.yHigh - yLow; ++i) {
                for (int j = rect.xLow - xLow; j < rect.xHigh - xLow; ++j) {
                    matrix[i][j] = true;
                }
            }
        }
        long answer = 0;
        for (auto row : matrix) {
            for (auto cell : row) {
                if (cell) {
                    answer += 1;
                }
            }
        }
        cout << answer << "\n";
    }
}