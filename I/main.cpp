//#define DEBUG
#ifndef DEBUG
#define NDEBUG
#else
#define _GLIBCXX_CONCEPT_CHECKS
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <unordered_map>
#include <map>
#include <set>
#include <unistd.h>
#include <fcntl.h>
#include <functional>

using namespace std;

template<typename Numeric>
struct Point;
template<typename Numeric>
struct Vector;

template<typename Numeric>
struct Point
{
    Numeric x, y;

    Point() = default;

    Point(Numeric x, Numeric y) : x(x), y(y)
    { }

    Point(const Point &) = default;

    template<typename Numeric2>
    Point(const Point<Numeric2> &o) : x(o.x), y(o.y)
    { }

    Point(const pair<Numeric, Numeric> &p) : Point(p.first, p.second)
    { }

    template<typename Numeric2>
    Point(const pair<Numeric2, Numeric2> &p) : Point(p.first, p.second)
    { }
};

template<typename Numeric>
static bool operator==(Point<Numeric> lhs, Point<Numeric> rhs)
{
    return (lhs.x == rhs.x) && (lhs.y == rhs.y);
}

namespace std {
    template<typename Numeric>
    struct hash<Point<Numeric>>
    {
        typedef Point<Numeric> argument_type;
        typedef size_t result_type;

        size_t operator()(Point<Numeric> p) const
        {
            hash<Numeric> hashn;
            auto a = hashn(p.x);
            size_t b = hashn(p.y);
            return a ^ (b << 1);
        }
    };
}


template<typename Numeric>
static istream &operator>>(const istream &in, Point<Numeric> &p)
{
    return in >> p.x >> p.y;
}

template<typename Numeric>
struct Vector
{
    Numeric x, y;

    Vector() = default;

    Vector(Numeric x, Numeric y) : x(x), y(y)
    { }

    Vector(const Vector &) = default;

    template<typename Numeric2>
    Vector(const Vector<Numeric2> &o) : x(o.x), y(o.y)
    { }

    explicit Vector(const Point<Numeric> &p) : Vector(p.x, p.y)
    { }
};

template<typename Numeric>
static Vector<Numeric> operator-(Point<Numeric> a, Point<Numeric> b)
{
    return {a.x - b.x, a.y - b.y};
}

template<typename Numeric>
static Numeric crossProduct(Vector<Numeric> a, Vector<Numeric> b)
{
    Numeric ax = a.x;
    Numeric ay = a.y;
    Numeric bx = b.x;
    Numeric by = b.y;
    return ax * by - ay * bx;
}

// (+/-)2*pole trójkąta
// b, c - podstawa; a - wierzchołek
template<typename Numeric>
static Numeric areaTangled(Point<Numeric> a, Point<Numeric> b, Point<Numeric> c)
{
    return crossProduct(b - a, c - a);
}

static long distance2(Point<int> a, Point<int> b)
{
    long x = b.x - a.x;
    long y = b.y - a.y;
    return x * x + y * y;
}

static double distance(Point<int> a, Point<int> b)
{
    return hypot(b.x - a.x, b.y - a.y);
}

static double distance(Point<double> a, Point<double> b)
{
    return hypot(b.x - a.x, b.y - a.y);
}

// >0 - lewo
// =0 - prosto
// <0 - prawo
template<typename Numeric>
static Numeric ccw(Point<Numeric> a, Point<Numeric> b, Point<Numeric> c)
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
    Numeric ax = a.x;
    Numeric ay = a.y;
    Numeric bx = b.x;
    Numeric by = b.y;
    Numeric cx = c.x;
    Numeric cy = c.y;
    return ax * by + bx * cy + cx * ay - cx * by - ax * cy - bx * ay;
}

template<typename Numeric>
static bool turnsLeft(Point<Numeric> a, Point<Numeric> b, Point<Numeric> c)
{
    return ccw(a, b, c) > 0;
}

template<typename Numeric>
static bool turnsStraight(Point<Numeric> a, Point<Numeric> b, Point<Numeric> c)
{
    return ccw(a, b, c) == 0;
}

template<typename Numeric>
static bool turnsRight(Point<Numeric> a, Point<Numeric> b, Point<Numeric> c)
{
    return ccw(a, b, c) < 0;
}

// potrzeba przynajmniej 2 punktów
// oszczędna - punkty tylko w wierzchołkach
// punkty muszą być unikalne i posortowane
// pierwszy punkt jest też ostatnim
static vector<Point<int>> convexHull(const vector<pair<int, int>> &points)
{
    // przy przerabianiu na wersje z punktami "na wprost" uwaga na dwukrotne przetworzenie ostatnigo punktu
    vector<Point<int>> result;
    for (auto p : points) {
        while (result.size() > 1 && !turnsRight<long>(result[result.size() - 2], result[result.size() - 1], p)) {
            result.pop_back();
        }
        result.push_back(p);
    }
    const auto limit = result.size();
    for (auto i = ++points.crbegin(), e = points.crend(); i != e; ++i) {
        // copy-paste
        while (result.size() > limit && !turnsRight<long>(result[result.size() - 2], result[result.size() - 1], *i)) {
            result.pop_back();
        }
        result.push_back(*i);
    }
    return result; // copy elision
}

// Składa górną pół-otoczke z dwóch górnych pół-otoczek
static vector<Point<int>> buildUpperHull(const vector<Point<int>> &left, const vector<Point<int>> &right)
{
    vector<Point<int>> result(left);
    for (auto p : right) {
        while (result.size() > 1 && !turnsRight<long>(result[result.size() - 2], result[result.size() - 1], p)) {
            result.pop_back();
        }
        result.push_back(p);
    }
    return result;
}

// a,b,c - współliniowe
// czy b leży pomiędzy a i c
// 17_geometria1/4
template<typename Numeric>
static bool onSegment(Point<Numeric> a, Point<Numeric> c, Point<Numeric> b)
{
    if ((min(a.x, c.x) <= b.x) && (b.x <= max(a.x, c.x)) && (min(a.y, c.y) <= b.y) && (b.y <= max(a.y, c.y))) {
        return true;
    }
    return false;
}

// Czy odcinki p0-p1 i p2-p3 przecinają się (również pojedyńcze punkty)
template<typename Numeric>
static bool segmentIntersect(Point<Numeric> p0, Point<Numeric> p1, Point<Numeric> p2, Point<Numeric> p3)
{
    auto I1 = ccw(p2, p3, p0);
    auto I2 = ccw(p2, p3, p1);
    auto I3 = ccw(p0, p1, p2);
    auto I4 = ccw(p0, p1, p3);
    if (I1 * I2 < 0 && I3 * I4 < 0) {
        return true;
    } else if (I1 == 0 && onSegment(p2, p3, p0)) {
        return true;
    } else if (I2 == 0 && onSegment(p2, p3, p1)) {
        return true;
    } else if (I3 == 0 && onSegment(p0, p1, p2)) {
        return true;
    } else if (I4 == 0 && onSegment(p0, p1, p3)) {
        return true;
    } else {
        return false;
    }
}

static long gcd(long a, long b)
{
    if (a <= 1 || b <= 1) {
        return 1;
    } else {
        return __gcd(a, b);
    }
}

static __int128 gcd(__int128 a, __int128 b)
{
    if (a <= 1 || b <= 1) {
        return 1;
    } else {
        return __gcd(a, b);
    }
}

static __int128 abs(__int128 n)
{
    if (n >= 0) {
        return n;
    } else {
        return -n;
    }
}

class Rational
{
public:
    __int128 numerator;
    __int128 denominator;

    Rational(__int128 numerator = 0, __int128 denominator = 1) : numerator(numerator), denominator(denominator)
    {
        assert(this->denominator > 0);
    }

    Rational reciprocal() const
    {
        return Rational(denominator, numerator);
    }

    Rational normalized() const
    {
        auto d = gcd(denominator, abs(numerator));
        return Rational(numerator / d, denominator / d);
    }
};

static Rational operator+(Rational lhs, Rational rhs)
{
    return Rational(lhs.numerator * rhs.denominator + rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

static Rational operator-(Rational lhs, Rational rhs)
{
    return Rational(lhs.numerator * rhs.denominator - rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

static Rational operator*(Rational lhs, Rational rhs)
{
    return Rational(lhs.numerator * rhs.numerator, lhs.denominator * rhs.denominator);
}

static Rational operator/(Rational lhs, Rational rhs)
{
    assert(rhs.numerator != 0);
    if (rhs.numerator >= 0) {
        return Rational(lhs.numerator * rhs.denominator, lhs.denominator * rhs.numerator);
    } else {
        return Rational(-(lhs.numerator * rhs.denominator), -(lhs.denominator * rhs.numerator));
    }
}

static bool operator==(Rational lhs, Rational rhs)
{
    return lhs.numerator * rhs.denominator == rhs.numerator * lhs.denominator;
}

static bool operator!=(Rational lhs, Rational rhs)
{
    return lhs.numerator != rhs.numerator || lhs.denominator != rhs.denominator;
}

static bool operator<(Rational lhs, Rational rhs)
{
    return lhs.numerator * rhs.denominator < rhs.numerator * lhs.denominator;
}

static bool operator>(Rational lhs, Rational rhs)
{
    return rhs < lhs;
}

static bool operator<=(Rational lhs, Rational rhs)
{
    return lhs.numerator * rhs.denominator <= rhs.numerator * lhs.denominator;
}

static ostream &operator<<(ostream &out, Rational n)
{
    return out << (long) n.numerator << '/' << (long) n.denominator;
}

static Point<Rational> crossLine(Rational a, Rational b, Rational c, Rational d)
{
    auto x = (d - b) / (a - c);
    auto y = a * x + b;
    return {x.normalized(), y.normalized()};
}

// Punkt przecięcia odcinków p1-p2 i p3-p4
static Point<Rational> crossSegment(const Point<Rational> &p1, const Point<Rational> &p2, const Point<Rational> &p3,
                                    const Point<Rational> &p4)
{
    Rational b = ((p2.y * p1.x - p1.y * p2.x) / (p1.x - p2.x)).normalized();
    // Nie ma odcinków pionowych
    Rational a = ((p1.x != 0) ? (p1.y - b) / p1.x : (p2.y - b) / p2.x).normalized();
    Rational d = ((p4.y * p3.x - p3.y * p4.x) / (p3.x - p4.x)).normalized();
    Rational c = ((p3.x != 0) ? (p3.y - d) / p3.x : (p4.y - d) / p4.x).normalized();
    return crossLine(a, b, c, d);
}

static Rational getY(const Point<Rational> &p1, const Point<Rational> &p2, Rational x)
{
    Rational b = ((p2.y * p1.x - p1.y * p2.x).normalized() / (p1.x - p2.x)).normalized();
    // Nie ma odcinków pionowych
    Rational a = ((p1.x != 0) ? (p1.y - b) / p1.x : (p2.y - b) / p2.x).normalized();
    return (a * x + b).normalized();
}

namespace std {
    template<>
    struct hash<Rational>
    {
        typedef Rational argument_type;
        typedef size_t result_type;

        size_t operator()(Rational r) const
        {
            size_t a = r.denominator;
            return r.numerator ^ (a << 1);
        }
    };
}

namespace std {
    template<>
    struct hash<Point<int>>
    {
        size_t operator()(const Point<int> &p) const noexcept
        {
            return p.x | static_cast<size_t>(p.y) << 32;
        }
    };
};


/******************* ROZWIĄZANIE *****************/

// w pozycji ogólnej
static vector<Point<double>> input;
// -1 -> n
// -2 -> n+1
// -3 -> n+2
static int n;

// typy danych:
// trójkąt - 3 punkty?
// krawędź - 2 punkty

struct Triangle
{
    int a, b, c;

    Triangle() : a(-1), b(-1), c(-1)
    { }

    Triangle(int a, int b, int c);
};

struct Edge
{
    static_assert(sizeof(int) == 4, "int must be 32bits wide");
    int a, b;

    Edge(int a, int b);
};

Triangle::Triangle(int a, int b, int c)
{
    assert(a >= 0);
    assert(b >= 0);
    assert(c >= 0);
    using std::swap;
    if (a > b) {
        swap(a, b);
    }
    if (b > c) {
        swap(b, c);
    }
    if (a > b) {
        swap(a, b);
    }
    assert(a < b);
    assert(b < c);
    assert(c < (1 << 21));
    this->a = a;
    this->b = b;
    this->c = c;
    assert(this->a >= 0);
    assert(this->b >= 0);
    assert(this->c >= 0);
}

static bool operator==(const Triangle &lhs, const Triangle &rhs)
{
    return lhs.a == rhs.a && lhs.b == rhs.b && lhs.c == rhs.c;
}

namespace std {
    template<>
    struct hash<Triangle>
    {
        size_t operator()(const Triangle &t) const noexcept
        {
            return (static_cast<size_t>(t.a) << 42) | static_cast<size_t >(t.b << 21) | t.c;
        }
    };
}


Edge::Edge(int a, int b)
{
    using std::swap;
    if (a > b) {
        swap(a, b);
    }
    assert(a < b);
    this->a = a;
    this->b = b;
}

static bool operator==(const Edge &lhs, const Edge &rhs)
{
    return lhs.a == rhs.a && lhs.b == rhs.b;
}

namespace std {
    template<>
    struct hash<Edge>
    {
        size_t operator()(const Edge &e) const noexcept
        {
            return static_cast<size_t>(e.a) << 32 | e.b;
        }
    };
}

static bool isCcw(Point<double> p1, Point<double> p2, Point<double> p3)
{
    return turnsLeft(p1, p2, p3) && turnsLeft(p2, p3, p1) && turnsLeft(p3, p1, p2);
}

static bool isIn(int pointId, const Triangle &triangle)
{
    // możliwy przypadek zdegenerowany!
    auto point = input[pointId];
    auto a = input[triangle.a];
    auto b = input[triangle.b];
    auto c = input[triangle.c];
    using std::swap;
    if (!isCcw(a, b, c)) {
        swap(b, c);
    }
    return !(turnsLeft(a, point, b) || turnsLeft(b, point, c) || turnsLeft(c, point, a));
}

static double det(const double M[][3])
{
    return M[0][0] * M[1][1] * M[2][2] + M[0][1] * M[1][2] * M[2][0] + M[0][2] * M[1][0] * M[2][1] -
           M[0][2] * M[1][1] * M[2][0] - M[0][1] * M[1][0] * M[2][2] - M[0][0] * M[1][2] * M[2][1];
}

static bool isLegal(int i, int j, int k, int l)
{
    // i-j - krawędź
    // r - wierzchołek
    // k - wierzchołek po przeciwnej stronie
    assert(i < j);
    assert(i >= 0);
    assert(j >= 0);
    assert(k >= 0);
    assert(l >= 0);
    if (!segmentIntersect(input[i], input[j], input[k], input[l])) {
        return true;
    }
    if (i >= n) {
        return true;
    } else if ((j < n) && (max(k, l) < n)) {
        // przypadek ogólny...
        Triangle t = {i, j, k};
        if (!isCcw(input[t.a], input[t.b], input[t.c])) {
            using std::swap;
            swap(t.b, t.c);
        }
        assert(isCcw(input[t.a], input[t.b], input[t.c]));
        const double Ax = input[t.a].x;
        const double Ay = input[t.a].y;
        const double Bx = input[t.b].x;
        const double By = input[t.b].y;
        const double Cx = input[t.c].x;
        const double Cy = input[t.c].y;
        const double Dx = input[l].x;
        const double Dy = input[l].y;
        const double M[3][3] = {{Ax - Dx, Ay - Dy, (Ax * Ax - Dx * Dx) + (Ay * Ay - Dy * Dy)},
                                {Bx - Dx, By - Dy, (Bx * Bx - Dx * Dx) + (By * By - Dy * Dy)},
                                {Cx - Dx, Cy - Dy, (Cx * Cx - Dx * Dx) + (Cy * Cy - Dy * Dy)}};
        bool result = det(M) < 0;
        if (!result) {
            assert(segmentIntersect(input[i], input[j], input[k], input[l]));
        }

        return result;
    } else {
        int negCount = 0;
        if (i >= n) {
            negCount += 1;
        }
        if (j >= n) {
            negCount += 1;
        }
        if (k >= n) {
            negCount += 1;
        }
        if (l >= n) {
            negCount += 1;
        }
        if (negCount == 1) {
            return j < n;
        } else {
            assert(negCount == 2);
            return max(i, j) > max(k, l);
        }
    }
}

struct WeightedEdge
{
    int a, b;
    double w;
};

static bool operator<(const WeightedEdge &lhs, const WeightedEdge &rhs)
{
    return lhs.w < rhs.w;
}

class FindUnion
{
    vector<int> parent;
    vector<int> rank;

public:
    explicit FindUnion(int n);

    bool tryUnion(int a, int b);

private:
    int find(int a);
};

FindUnion::FindUnion(int n) : parent(n), rank(n)
{
    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }
}

bool FindUnion::tryUnion(int a, int b)
{
    a = find(a);
    b = find(b);
    if (a == b) {
        return false;
    }
    if (rank[a] == rank[b]) {
        parent[b] = a;
        rank[b] = rank[a] + 1;
        return true;
    }
    using std::swap;
    if (rank[a] > rank[b]) {
        swap(a, b);
    }
    parent[a] = b;
    return true;
}

int FindUnion::find(int a)
{
    if (parent[a] == parent[parent[a]]) {
        return parent[a];
    }
    return parent[a] = find(parent[a]);
}

class EdgeMapper
{
    unordered_map<Edge, array<int, 2>> edgeTo2Triangles;

public:
    void initEdge(Edge e, int a, int b);

    void update(Edge e, int from, int to);

    bool has(Edge e) const;

    int oppositePoint(Edge e, int p) const;

    void remove(Edge e);

    vector<WeightedEdge> getEdges() const;

private:
    bool sidesAreCorrect(Edge e) const;
};

void EdgeMapper::initEdge(Edge e, int a, int b)
{
    assert(a != b);
    assert(edgeTo2Triangles.find(e) == edgeTo2Triangles.end());
    edgeTo2Triangles[e] = {a, b};
    assert(sidesAreCorrect(e));
}

void EdgeMapper::update(Edge e, int from, int to)
{
    edgeTo2Triangles[e][edgeTo2Triangles[e][0] != from] = to;
    assert(sidesAreCorrect(e));
}

bool EdgeMapper::has(Edge e) const
{
    return edgeTo2Triangles.find(e) != edgeTo2Triangles.end();
}

int EdgeMapper::oppositePoint(Edge e, int p) const
{
    assert(has(e));
    auto i = edgeTo2Triangles.find(e);
    int k = (i->second)[(i->second)[0] == p];
    return k;
}

void EdgeMapper::remove(Edge e)
{
    assert(has(e));
    edgeTo2Triangles.erase(e);
}

bool EdgeMapper::sidesAreCorrect(Edge e) const
{
    assert(has(e));
    auto i = edgeTo2Triangles.find(e);
    int p1 = i->second[0];
    int p2 = i->second[1];
    if (p1 == -1 || p2 == -1) {
        return true; // nie można sprawdzić
    }
    auto ea = input[e.a];
    auto eb = input[e.b];
    auto a = input[p1];
    auto b = input[p2];
    if (turnsStraight(ea, eb, a)) {
        return false;
    }
    if (turnsStraight(ea, eb, b)) {
        return false;
    }
    if (turnsLeft(ea, eb, a)) {
        return turnsRight(ea, eb, b);
    } else {
        return turnsRight(ea, eb, a) && turnsLeft(ea, eb, b);
    }
}

vector<WeightedEdge> EdgeMapper::getEdges() const
{
    vector<WeightedEdge> edges;
    for (auto edge : edgeTo2Triangles) {
        int a = edge.first.a;
        int b = edge.first.b;
        assert(a < b);
        if (b < n) {
            edges.push_back({a, b, distance(input[a], input[b])});
        }
    }
    return edges;
}

// czy można przekręcić i-j na r-k (oba odcinki przecinają się)
static bool canBeSwapped(int i, int j, int r, int k)
{
    return segmentIntersect(input[i], input[j], input[r], input[k]);
}

static double getMstValue(vector<WeightedEdge> edges)
{
    sort(edges.begin(), edges.end());
    FindUnion findUnion(n);
    double result = .0;
    for (const auto &e : edges) {
        if (findUnion.tryUnion(e.a, e.b)) {
            result += e.w;
        }
    }
    return result;
}

static vector<WeightedEdge> generateEdgesFromInput()
{
    vector<WeightedEdge> result;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            result.push_back({i, j, distance(input[i], input[j])});
        }
    }
    return result;
}

// struktury danych:
// punkt wejściowy -> trójkąt go zawierający
// trójkąt -> zbiór punktów w nim zawartych
// krawędź -> 2 trójkąty

int main()
{
    // TODO: sprawdzić kalibracje isLegal (szczególnie wyznacznika)
    //close(0);
    //open("/home/local/AG/I/testy/t00.in", O_RDONLY);
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        cin >> n;
        input.resize(n);
        for (auto &p : input)
            cin >> p.x >> p.y;
        mt19937 engine(404);
        //shuffle(input.begin(), input.end(), engine);
        double M = .0;
        for (auto p : input) {
            M = max(M, max(abs(p.x), abs(p.y)));
        }
        uniform_real_distribution<double> eps(0, 1);
        input.push_back({4 * M + eps(engine), eps(engine)});
        input.push_back({eps(engine), 4 * M + eps(engine)});
        input.push_back({-4 * M - eps(engine), -4 * M - eps(engine)});
        vector<Triangle> pointIdToTriangle(n);
        unordered_map<Triangle, vector<int>> triangleToPointsIds;
        EdgeMapper edgeTo2Triangles;
        function<void(int, int, int)> legalizeEdge;
        legalizeEdge = [&pointIdToTriangle, &triangleToPointsIds, &edgeTo2Triangles, &legalizeEdge](int r, int i,
                                                                                                    int j) {
            using std::swap;
            if (i > j) {
                swap(i, j);
            }
            assert(i < j);
            assert(triangleToPointsIds.find({r, i, j}) != triangleToPointsIds.end());
            assert(edgeTo2Triangles.has({i, j}));
            int k = edgeTo2Triangles.oppositePoint({i, j}, r);
            if (k != -1)
                assert(triangleToPointsIds.find({k, i, j}) != triangleToPointsIds.end());
            assert(edgeTo2Triangles.has({i, r}));
            assert(edgeTo2Triangles.has({j, r}));
            if (k != -1) {
                assert(edgeTo2Triangles.has({i, k}));
                assert(edgeTo2Triangles.has({j, k}));
            }
            if ((k != -1) && (!isLegal(i, j, r, k))) {
                vector<int> part1ToReassign = move(triangleToPointsIds[{i, j, r}]);
                vector<int> part2ToReassign = move(triangleToPointsIds[{i, j, k}]);
                assert(canBeSwapped(i, j, k, r));
                triangleToPointsIds.erase({i, j, r});
                triangleToPointsIds.erase({i, j, k});
                assert(r != k);
                edgeTo2Triangles.remove({i, j});
                edgeTo2Triangles.initEdge({r, k}, i, j);
                edgeTo2Triangles.update({r, i}, j, k);
                edgeTo2Triangles.update({i, k}, j, r);
                edgeTo2Triangles.update({k, j}, i, r);
                edgeTo2Triangles.update({j, r}, i, k);

                Triangle e1 = {i, j, r};
                Triangle e2 = {i, j, k};
                //cerr << "Kasowany trójkąt (" << -1 << "): {" << e1.a << ", " << e1.b << ", " << e1.c << "}\n";
                //cerr << "Kasowany trójkąt (" << -1 << "): {" << e2.a << ", " << e2.b << ", " << e2.c << "}\n";
                Triangle c1 = {r, k, j};
                Triangle c2 = {r, k, i};
                //cerr << "Tworzony trójkąt (" << -1 << "): {" << c1.a << ", " << c1.b << ", " << c1.c << "}\n";
                //cerr << "Tworzony trójkąt (" << -1 << "): {" << c2.a << ", " << c2.b << ", " << c2.c << "}\n";

                // touch
                triangleToPointsIds[{r, k, j}];
                triangleToPointsIds[{r, k, i}];
                auto assignPoint = [&pointIdToTriangle, &triangleToPointsIds, r, i, j, k](int p) {
                    Triangle t = {r, k, j};
                    if (!isIn(p, t)) {
                        t = {r, k, i};
                        assert(isIn(p, t));
                    }
                    pointIdToTriangle[p] = t;
                    triangleToPointsIds[t].push_back(p);
                };
                for (int p : part1ToReassign) {
                    assignPoint(p);
                }
                for (int p :part2ToReassign) {
                    assignPoint(p);
                }
                legalizeEdge(r, i, k);
                legalizeEdge(r, k, j);
            }
        };

        for (int i = 0; i < n; ++i) {
            pointIdToTriangle[i] = {n, n + 1, n + 2};
            triangleToPointsIds[{n, n + 1, n + 2}].push_back(i);
        }
        edgeTo2Triangles.initEdge({n, n + 1}, n + 2, -1);
        edgeTo2Triangles.initEdge({n + 1, n + 2}, n, -1);
        edgeTo2Triangles.initEdge({n + 2, n}, n + 1, -1);
        for (int i = 0; i < n; ++i) {
            Triangle t = pointIdToTriangle[i];
            //cerr << "Kasowany trójkąt (" << i << "): {" << t.a << ", " << t.b << ", " << t.c << "}\n";
            assert(isIn(i, t));
            // stwórz szkielet nowych trójkątów
            Triangle t1 = {t.a, t.b, i};
            Triangle t2 = {t.b, t.c, i};
            Triangle t3 = {t.c, t.a, i};
            //cerr << "Tworzony trójkąt (" << i << "): {" << t1.a << ", " << t1.b << ", " << t1.c << "}\n";
            //cerr << "Tworzony trójkąt (" << i << "): {" << t2.a << ", " << t2.b << ", " << t2.c << "}\n";
            //cerr << "Tworzony trójkąt (" << i << "): {" << t3.a << ", " << t3.b << ", " << t3.c << "}\n";

            // zmień strukture
            vector<int> pointToReassign = move(triangleToPointsIds[t]);
            triangleToPointsIds.erase(t);
            edgeTo2Triangles.update({t.a, t.b}, t.c, i);
            edgeTo2Triangles.update({t.b, t.c}, t.a, i);
            edgeTo2Triangles.update({t.c, t.a}, t.b, i);
            edgeTo2Triangles.initEdge({t.a, i}, t.b, t.c);
            edgeTo2Triangles.initEdge({t.b, i}, t.a, t.c);
            edgeTo2Triangles.initEdge({t.c, i}, t.a, t.b);

            // touch
            triangleToPointsIds[t1];
            triangleToPointsIds[t2];
            triangleToPointsIds[t3];
            for (int pId : pointToReassign) {
                assert(isIn(pId, t));
                if (pId == i) {
                    continue;
                }
                Triangle ts;
                if (isIn(pId, t1)) {
                    ts = t1;
                } else if (isIn(pId, t2)) {
                    ts = t2;
                } else {
                    assert(isIn(pId, t3));
                    ts = t3;
                }
                pointIdToTriangle[pId] = ts;
                triangleToPointsIds[ts].push_back(pId);
            }

            // usuń nielegalne krawędzie
            legalizeEdge(i, t.a, t.b);
            legalizeEdge(i, t.b, t.c);
            legalizeEdge(i, t.a, t.c);
        }

        // MST
        vector<WeightedEdge> edges = edgeTo2Triangles.getEdges();
#ifdef DEBUG
        for (const auto &e : edges) {
            cerr << "DE: " << e.a << " " << e.b << "\n";
        }
#endif
        double result = getMstValue(move(edges));
#ifdef DEBUG
        edges = generateEdgesFromInput();
        double forceResult = getMstValue(move(edges));
        if (abs(result - forceResult) >= 0.000000001) {
            cerr << result << " != " << forceResult << "\n";
            assert(false);
        }
#endif

        cout << fixed << setprecision(10) << result << '\n';
    }
    return 0;
}