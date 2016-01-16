//#define DEBUG
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
#include <stack>

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

class Rational
{
public:
    long numerator;
    long denominator;

    Rational(long numerator = 0, long denominator = 1) : numerator(numerator / __gcd(denominator, abs(numerator))),
                                                         denominator(denominator / __gcd(denominator, abs(numerator)))
    {
        assert(this->denominator > 0);
    }

    Rational reciprocal() const
    {
        return Rational(denominator, numerator);
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
    return lhs.numerator == rhs.numerator && lhs.denominator == rhs.denominator;
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
    return out << n.numerator << '/' << n.denominator;
}

static Point<Rational> crossLine(Rational a, Rational b, Rational c, Rational d)
{
    auto x = (d - b) / (a - c);
    auto y = a * x + b;
    return {x, y};
}

// Punkt przecięcia odcinków p1-p2 i p3-p4
static Point<Rational> crossSegment(const Point<Rational> &p1, const Point<Rational> &p2, const Point<Rational> &p3,
                                    const Point<Rational> &p4)
{
    Rational b = (p2.y * p1.x - p1.y * p2.x) / (p1.x - p2.x);
    // Nie ma odcinków pionowych
    Rational a = (p1.x != 0) ? (p1.y - b) / p1.x : (p2.y - b) / p2.x;
    Rational d = (p4.y * p3.x - p3.y * p4.x) / (p3.x - p4.x);
    Rational c = (p3.x != 0) ? (p3.y - d) / p3.x : (p4.y - d) / p4.x;
    return crossLine(a, b, c, d);
}

static Rational getY(const Point<Rational> &p1, const Point<Rational> &p2, Rational x)
{
    Rational b = (p2.y * p1.x - p1.y * p2.x) / (p1.x - p2.x);
    // Nie ma odcinków pionowych
    Rational a = (p1.x != 0) ? (p1.y - b) / p1.x : (p2.y - b) / p2.x;
    return a * x + b;
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

/******************* ROZWIĄZANIE *****************/

static bool operator<(const Point<long> &lhs, const Point<long> &rhs)
{
    return lhs.y > rhs.y || (lhs.y == rhs.y && lhs.x < rhs.x);
}

static bool isCCW(Point<long> p1, Point<long> p2, Point<long> p3)
{
    return turnsLeft<long>(p1, p2, p3) && turnsLeft<long>(p2, p3, p1) && turnsLeft<long>(p3, p1, p2);
}

// DEBUG
static void checkCCW(Point<long> p1, Point<long> p2, Point<long> p3)
{
    assert(isCCW(p1, p2, p3));
}

static void printCCW(const vector<Point<long>> &input, int p1, int p2, int p3)
{
    if (isCCW(input[p1], input[p2], input[p3])) {
        cout << p1 << ' ' << p2 << ' ' << p3 << '\n';
    } else {
        assert(isCCW(input[p1], input[p3], input[p2]));
        cout << p1 << ' ' << p3 << ' ' << p2 << '\n';
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
        cout << n - 2 << '\n';
        vector<Point<long>> input(n);
        for (auto &p : input) {
            cin >> p.x >> p.y;
        }
        const int topIndex = distance(input.begin(), min_element(input.begin(), input.end()));
        const int bottomIndex = distance(input.begin(), max_element(input.begin(), input.end()));
        auto isOnLeftSide = [topIndex, bottomIndex](int index) {
            if (topIndex < bottomIndex) {
                return index > topIndex && index <= bottomIndex;
            } else {
                return index > topIndex || index <= bottomIndex;
            }
        };
        vector<int> leftSide, rightSide;
        if (topIndex < bottomIndex) {
            for (int i = topIndex + 1; i <= bottomIndex; ++i)
                leftSide.push_back(i);
            for (int i = bottomIndex + 1; i < n; ++i)
                rightSide.push_back(i);
            for (int i = 0; i <= topIndex; ++i)
                rightSide.push_back(i);
        } else {
            for (int i = topIndex + 1; i < n; ++i)
                leftSide.push_back(i);
            for (int i = 0; i <= bottomIndex; ++i)
                leftSide.push_back(i);
            for (int i = bottomIndex + 1; i <= topIndex; ++i)
                rightSide.push_back(i);
        }
        vector<int> orderedPoints;
        merge(leftSide.begin(), leftSide.end(), rightSide.rbegin(), rightSide.rend(), back_inserter(orderedPoints),
              [&input](int lhs, int rhs) {
                  return input[lhs] < input[rhs];
              });
        assert(orderedPoints.size() == input.size());
#ifdef DEBUG
        for (int i = 0; i < orderedPoints.size() - 1; ++i) {
            assert(input[orderedPoints[i]] < input[orderedPoints[i + 1]]);
        }
#endif
        vector<int> S;
        auto isOnStackSide = [&S, n](int index) {
            assert(!S.empty());
            return (index == (S.back() + 1) % n) || (index == (n + S.back() - 1) % n);
        };
        int count = 0;
        S.push_back(orderedPoints[0]);
        S.push_back(orderedPoints[1]);
        for (int i = 2; i < n - 1; ++i) {
            assert(!S.empty());
            int ui = orderedPoints[i];
//            cout << "Stos:";
//            for (auto p : S)
//                cout << " " << p;
//            cout << "\nPunkt: " << ui << " " << (isOnLeftSide(ui) ? "L" : "R") << "\n";
            if (!isOnStackSide(ui)) {
                assert(!S.empty());
                const int o = S.back();
                while (S.size() > 1) {
                    int t = S.back();
                    S.pop_back();
                    printCCW(input, ui, t, S.back());
                    count += 1;
                }
                S.pop_back();
                assert(S.empty());
                S.push_back(ui);
                S.push_back(o);
            } else {
                if (isOnLeftSide(ui)) {
                    assert(isOnLeftSide(S.back()));
                    while (S.size() > 1 && turnsRight<long>(input[ui], input[S.back()], input[S[S.size() - 2]])) {
                        int t = S.back();
                        S.pop_back();
                        printCCW(input, ui, t, S.back());
                        count += 1;
                    }
                } else {
                    assert(!isOnLeftSide(S.back()));
                    while (S.size() > 1 && turnsLeft<long>(input[ui], input[S.back()], input[S[S.size() - 2]])) {
                        int t = S.back();
                        S.pop_back();
                        printCCW(input, ui, t, S.back());
                        count += 1;
                    }
                }
                S.push_back(ui);
            }
        }
        int un = orderedPoints[n - 1];
        assert(isOnStackSide(un));
        assert(S.size() > 1);
        for (int i = 0; i < S.size() - 1; ++i) {
            printCCW(input, un, S[i], S[i + 1]);
            count += 1;
        }
        assert(count == n - 2);
    }
    return 0;
}