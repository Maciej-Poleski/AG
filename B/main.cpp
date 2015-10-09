#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
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
    return {a.x - b.x, a.y - b.y};
}

static long crossProduct(Vector a, Vector b)
{
    long ax = a.x;
    long ay = a.y;
    long bx = b.x;
    long by = b.y;
    return ax * by - ay * bx;
}

// (+/-)2*pole trójkąta
// b, c - podstawa; a - wierzchołek
static long areaTangled(Point a, Point b, Point c)
{
    return crossProduct(b - a, c - a);
}

static long distance2(Point a, Point b)
{
    long x = b.x - a.x;
    long y = b.y - a.y;
    return x * x + y * y;
}

static double distance(Point a, Point b)
{
    return hypot(b.x - a.x, b.y - a.y);
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

// potrzeba przynajmniej 2 punktów
// oszczędna - punkty tylko w wierzchołkach
// punkty muszą być unikalne i posortowane
// pierwszy punkt jest też ostatnim
static vector<Point> convexHull(const vector<pair<int, int>> &points)
{
    // przy przerabianiu na wersje z punktami "na wprost" uwaga na dwukrotne przetworzenie ostatnigo punktu
    vector<Point> result;
    for (auto p : points) {
        while (result.size() > 1 && !turnsRight(result[result.size() - 2], result[result.size() - 1], p)) {
            result.pop_back();
        }
        result.push_back(p);
    }
    for (auto i = points.crbegin(), e = points.crend(); i != e; ++i) {
        // copy-paste
        while (result.size() > 1 && !turnsRight(result[result.size() - 2], result[result.size() - 1], *i)) {
            result.pop_back();
        }
        result.push_back(*i);
    }
    return result; // copy elision
}

/******************* ROZWIĄZANIE *****************/

static int cyclicNext(int i, size_t size)
{
    i += 1;
    if (i == size) {
        return 0;
    }
    return i;
}

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        cin >> n;
        vector<pair<int, int>> points(n);
        for (auto &p : points) {
            cin >> p.first >> p.second;
        }
        sort(points.begin(), points.end());
        points.erase(unique(points.begin(), points.end()), points.end());
        if (points.size() < 3) {
            cout << "0.0\n";
            continue;
        }
        const auto hull = convexHull(points);
        int j = 1;
        double dist = numeric_limits<double>::max();
        for (int i = 1; i < hull.size(); ++i) {
            while (areaTangled(hull[cyclicNext(j, hull.size())], hull[i], hull[i - 1]) >=
                   areaTangled(hull[j], hull[i], hull[i - 1])) {
                j = cyclicNext(j, hull.size());
            }
            dist = min(dist, areaTangled(hull[j], hull[i], hull[i - 1]) / distance(hull[i - 1], hull[i]));
        }
        cout << fixed << setprecision(7) << dist << '\n';
    }
    return 0;
}