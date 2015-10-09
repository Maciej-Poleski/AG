#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;

struct Point;
struct Vector;

struct Point
{
    int x, y;

    Point() = default;

    Point(int x, int y) : x(x), y(y)
    { }

    Point(const pair<int, int> &p) : Point(p.first, p.second)
    { }
};

struct Vector
{
    int x, y;

    Vector() = default;

    Vector(int x, int y) : x(x), y(y)
    { }

    explicit Vector(const Point &p) : Vector(p.x, p.y)
    { }
};

static Vector operator-(Point a, Point b)
{
    return {b.x - a.x, b.y - a.y};
}

// >0 - lewo
// =0 - prosto
// <0 - prawo
static long ccw(Point a, Point b, Point c)
{
//    Kalibracja:
//    cout<<ccw({0,0},{0,1},{0,1})<<'\n';
//    cout<<ccw({0,0},{0,1},{0,2})<<'\n';
//    cout<<ccw({0,0},{0,2},{0,1})<<'\n';
//    cout<<ccw({0,0},{0,1},{1,1})<<'\n';
//    cout<<ccw({0,0},{0,1},{1,2})<<'\n';
//    cout<<ccw({0,0},{0,1},{1,0})<<'\n';
//    cout<<ccw({0,0},{0,1},{-1,1})<<'\n';
//    cout<<ccw({0,0},{0,1},{-1,2})<<'\n';
//    cout<<ccw({0,0},{0,1},{-1,0})<<'\n';
    long ax = a.x;
    long ay = a.y;
    long bx = b.x;
    long by = b.y;
    long cx = c.x;
    long cy = c.y;
    return ax * by + bx * cy + cx * ay - cx * by - ax * cy - bx * ay;
}

static bool turnsLeft(Point a, Point b, Point c)
{
    return ccw(a, b, c) > 0;
}

static bool turnsStraight(Point a, Point b, Point c)
{
    return ccw(a, b, c) == 0;
}

static bool turnsRight(Point a, Point b, Point c)
{
    return ccw(a, b, c) < 0;
}


/******************* ROZWIĄZANIE *****************/

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        vector<pair<int, int>> upper, lower;
        cin >> n;
        for (int i = 0; i < n; ++i) {
            int x, a, b;
            cin >> x >> a >> b;
            upper.push_back({x, b});
            lower.push_back({x, a});
        }
        if (n == 1) {
            cout << "TAK\n";
            continue;
        }
        sort(upper.begin(), upper.end());
        sort(lower.begin(), lower.end());
        vector<Point> upperHull;
        for (auto p : upper) {
            while (upperHull.size() > 1 &&
                   !turnsLeft(upperHull[upperHull.size() - 2], upperHull[upperHull.size() - 1], p)) {
                upperHull.pop_back();
            }
            upperHull.push_back(p);
        }
        vector<Point> lowerHull;
        for (auto p : lower) {
            while (lowerHull.size() > 1 &&
                   !turnsRight(lowerHull[lowerHull.size() - 2], lowerHull[lowerHull.size() - 1], p)) {
                lowerHull.pop_back();
            }
            lowerHull.push_back(p);
        }
        int upperIterator = 1, lowerIterator = 1;
        bool ok = true;
        for (; upperIterator < upperHull.size() && lowerIterator < lowerHull.size();) {
            auto currentLower = lowerHull[lowerIterator - 1];
            auto nextLower = lowerHull[lowerIterator];
            auto currentUpper = upperHull[upperIterator - 1];
            auto nextUpper = upperHull[upperIterator];
            if (nextUpper.x <= nextLower.x) {
                if (turnsLeft(currentLower, nextLower, nextUpper)) {
                    upperIterator += 1;
                } else {
                    ok = false;
                    break;
                }
            } else {
                if (turnsRight(currentUpper, nextUpper, nextLower)) {
                    lowerIterator += 1;
                } else {
                    ok = false;
                    break;
                }
            }
        }
        cout << (ok ? "TAK\n" : "NIE\n");
    }
    return 0;
}