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

static bool operator<(const Point<Rational> &lhs, const Point<Rational> &rhs)
{
    return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}

struct Event
{
    vector<int> begins;
//    vector<int> intersections;
//    vector<int> ends;
//
//    int anything() const;
};

//int Event::anything() const
//{
//    if (!begins.empty()) {
//        return begins.front();
//    }
//    if (!intersections.empty()) {
//        return intersections.front();
//    }
//    assert(!ends.empty());
//    return ends.front();
//}

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        cin >> n;
        vector<pair<Point<Rational>, Point<Rational>>> input;
        for (int i = 0; i < n; ++i) {
            int ax, ay, bx, by;
            cin >> ax >> ay >> bx >> by;
            if (ax > bx) {
                swap(ax, bx);
                swap(ay, by);
            }
            input.push_back(make_pair(Point<Rational>(ax, ay), Point<Rational>(bx, by)));
        }
        input.emplace_back();
        map<int, long> output;
        map<Point<Rational>, Event> Q;
        for (int i = 0; i < n; ++i) {
            Q[input[i].first].begins.push_back(i);
            Q[input[i].second]; //.ends.push_back(i);
        }

        Rational currentX;   // for sweepper

        auto sweepComp = [&input](Rational x, int lhs, int rhs) {
            return getY(input[lhs].first, input[lhs].second, x) < getY(input[rhs].first, input[rhs].second, x);
        };
        auto sweepCompCurrentX = [&currentX, sweepComp](int lhs, int rhs) {
            return sweepComp(currentX, lhs, rhs);
        };
        multiset<int, decltype(sweepCompCurrentX)> M(sweepCompCurrentX);
        while (!Q.empty()) {
            auto point = Q.begin()->first;
            auto event = Q.begin()->second;
            Q.erase(Q.begin());
            clog << Q.size() << '\r';
//            clog << "Zdarzenie w " << point.x << " " << point.y << "\n";
//            clog << "Początek:\n";
//            for (auto l : event.begins)
//                clog << l << "\n";
//            clog << "Przecięcie:\n";
//            for (auto l : event.intersections)
//                clog << l << "\n";
//            clog << "Koniec:\n";
//            for (auto l : event.ends)
//                clog << l << "\n";
//            clog << "Miotła:";
//            for (auto p : M) {
//                clog << " " << p;
//            }
//            clog << "\n";
            currentX = point.x;
            auto &U = event.begins;
            input[n] = make_pair(Point<Rational>(currentX - 1, point.y), Point<Rational>(currentX + 1, point.y));
            auto equalOnSweep = M.equal_range(n);
            vector<int> L, C;
            for (auto i = equalOnSweep.first; i != equalOnSweep.second; ++i) {
                int index = *i;
                if (input[index].first == point || input[index].second == point) {
                    L.push_back(index);
                } else {
                    C.push_back(index);
                }
            }
            if (L.size() + U.size() + C.size() > 1) {
                output[L.size() + U.size() + C.size()] += 1;
//                clog << "Punkt przecięcia " << point.x << " " << point.y << "\n";
            }
            const auto hint = M.erase(equalOnSweep.first, equalOnSweep.second);
            const auto nextX = currentX + 1; // any x greater than current
            auto sweepCompNextX = [sweepComp, nextX](int lhs, int rhs) {
                return sweepComp(nextX, lhs, rhs);
            };
            sort(U.begin(), U.end(), sweepCompNextX);
#ifdef DEBUG
            for (int i = 0; i < static_cast<int>(C.size()) - 1; ++i) {
                assert(sweepCompNextX(C[i + 1], C[i]));
                assert(!sweepCompCurrentX(C[i], C[i + 1]));
                assert(!sweepCompCurrentX(C[i + 1], C[i]));
            }
#endif
            vector<int> merged;
            merge(U.begin(), U.end(), C.rbegin(), C.rend(), back_inserter(merged), sweepCompNextX);
            const auto sr = hint;
            auto sl = sr;
            if (sl != M.begin()) {
                advance(sl, -1);
            } else {
                sl = M.end();
            }
            for (auto i = merged.begin(), e = merged.end(); i != e; ++i) {
                M.insert(hint, *i); // inserting from least to greatest according to nextX
                // all equal according to currentX
                // multiset (C++11) inserts in upper_bound position
            }
            auto findNewEvent = [&input, &currentX, &Q](int sl, int sr, Rational y) {
                if (segmentIntersect(input[sl].first, input[sl].second, input[sr].first, input[sr].second)) {
                    auto p = crossSegment(input[sl].first, input[sl].second, input[sr].first, input[sr].second);
                    if (p.x > currentX || (p.x == currentX && p.y > y)) {
//                        clog << "Wykryto nowe przecięcie w " << p.x << " " << p.y << " między " << sl << " a " << sr <<
//                        "\n";
                        Q[p]; //.intersections.push_back(sl);
//                        Q[p].intersections.push_back(sr);
                    } /*else {
                        clog << "Zignorowano przecięcie w " << p.x << " " << p.y << " między " << sl << " a " << sr <<
                        "\n";
                    }*/
                }/* else {
                    clog << "Nie ma przecięcia między " << sl << " a " << sr << "\n";
                }*/
            };
            if (merged.empty()) {
                if (sr != M.end() && sr != M.begin()) {
                    findNewEvent(*sl, *sr, point.y);
                }
            } else {
                if (sl != M.end()) {
                    auto sp = merged.front();
                    findNewEvent(*sl, sp, point.y);
                }
                if (sr != M.end()) {
                    auto spp = merged.back();
                    findNewEvent(spp, *sr, point.y);
                }
            }
        }
        for (auto &p : output) {
            cout << p.first << ": " << p.second << '\n';
        }
    }
    return 0;
}