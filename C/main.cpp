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

    Point(const Point &) = default;

    Point(const pair<int, int> &p) : Point(p.first, p.second)
    { }
};

static istream &operator>>(istream &in, Point &p)
{
    return in >> p.x >> p.y;
}

struct Vector
{
    int x, y;

    Vector() = default;

    Vector(int x, int y) : x(x), y(y)
    { }

    Vector(const Vector &) = default;

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
    const auto limit = result.size();
    for (auto i = ++points.crbegin(), e = points.crend(); i != e; ++i) {
        // copy-paste
        while (result.size() > limit && !turnsRight(result[result.size() - 2], result[result.size() - 1], *i)) {
            result.pop_back();
        }
        result.push_back(*i);
    }
    return result; // copy elision
}

// Składa górną pół-otoczke z dwóch górnych pół-otoczek
static vector<Point> buildUpperHull(const vector<Point> &left, const vector<Point> &right)
{
    vector<Point> result(left);
    for (auto p : right) {
        while (result.size() > 1 && !turnsRight(result[result.size() - 2], result[result.size() - 1], p)) {
            result.pop_back();
        }
        result.push_back(p);
    }
    return result;
}

/******************* ROZWIĄZANIE *****************/

static vector<vector<vector<Point>>> dict;

// l - obecna warstwa, l+1 - nowa warstwa
static void buildNextLayer(int l)
{
    if (dict[l].size() == 1) {
        return;
    }
    assert(dict.size() == l + 1);
    dict.emplace_back();
    for (auto i = dict[l].begin(), e = dict[l].end(); i < e; i += 2) {
        if (i + 1 == e) {
            dict[l + 1].emplace_back(*i);
            continue;
        }
        dict[l + 1].push_back(buildUpperHull(*i, *(i + 1)));
    }
    buildNextLayer(l + 1);
}

static bool intersectsWithHull(Point a, Point b, const vector<Point> &hull)
{
    auto r = hull.size();
    decltype(r) l = 0;
    while (r - l > 1) {
        auto s = (l + r) / 2;
        auto d1 = -areaTangled(hull[s - 1], a, b) / distance(a, b);
        auto d2 = -areaTangled(hull[s], a, b) / distance(a, b);
        if (min(d1, d2) < 0) {
            return true;
        }
        if (d1 == d2) {
            return d1 < 0;
        }
        if (d1 < d2) {
            r = s;
        } else {
            l = s;
        }
    }
    if (-areaTangled(hull[l], a, b) / distance(a, b) < 0) {
        return true;
    }
    if ((r > l + 1) && (-areaTangled(hull[l + 1], a, b) / distance(a, b)) < 0) {
        return true;
    }
    return false;
}

static pair<int, int> takeNext(pair<int, int> hullId)
{
    hullId.second += 1;
    while (hullId.second % 2 == 0) {
        hullId.first += 1;
        hullId.second /= 2;
    }
    return hullId;
};

static bool isHullId(pair<int, int> hullId)
{
    return dict.size() > hullId.first && dict[hullId.first].size() > hullId.second;
}

static void findEdge(Point a, Point b, pair<int, int> hullId)
{
    if (hullId.first == 0) {
        cout << hullId.second + 1 << " ";
        return;
    }
    if (intersectsWithHull(a, b, dict[hullId.first - 1][hullId.second * 2])) {
        findEdge(a, b, make_pair(hullId.first - 1, hullId.second * 2));
    } else {
        findEdge(a, b, make_pair(hullId.first - 1, hullId.second * 2 + 1));
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        cin >> n;
        vector<Point> input(n);
        for (auto &p : input) {
            cin >> p;
        }
        dict.emplace_back();
        for (auto i = input.begin(), e = input.end() - 1; i != e; ++i) {
            dict[0].emplace_back(i, i + 2);
        }
        buildNextLayer(0);
        for (int i = 0, e = input.size() - 2; i < e; ++i) {
            Point a = input[i];
            Point b = input[i + 1];
            for (auto hullId = takeNext(make_pair(0, i)); ; hullId = takeNext(hullId)) {
                if (!isHullId(hullId)) {
                    cout << "0 ";
                    break;
                } else {
                    if (intersectsWithHull(a, b, dict[hullId.first][hullId.second])) {
                        findEdge(a, b, hullId);
                        break;
                    }
                }
            }
        }
        cout << "0\n";

        dict.clear();
    }
    return 0;
}