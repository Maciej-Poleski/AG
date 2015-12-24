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
#include <unistd.h>
#include <fcntl.h>

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

static bool operator<(Point<int> lhs, Point<int> rhs)
{
    return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}

class Column
{
public:
    Column();

    // posortowane współrzędne
    template<class Iterator>
    Column(Iterator begin, Iterator end);

    Column(const Column &) = delete;

    Column(Column &&) = default;

    Column &operator=(Column &&other);

    ~Column();

    void add(int y, int diff);

    // y1 < y2
    int query(int y1, int y2);

private:
    struct Node
    {
        int begin, end;
        int sum;
    };

    Node *tree;
    size_t size;

    size_t find(size_t root, int y);

    int queryImpl(size_t node, int y1, int y2);
};

template<typename Numeric>
Numeric mpow2(Numeric n)
{
    Numeric i;
    for (i = 1; i < n; i *= 2);
    return i;

}

Column::Column() : tree(nullptr), size(0)
{ }

template<class Iterator>
Column::Column(Iterator begin, Iterator end)
{
    auto length = distance(begin, end);
    this->size = mpow2(length);
    this->tree = new Node[this->size * 2];
    auto i = this->size;
    for (auto j = begin; j != end; ++j, ++i) {
        this->tree[i].begin = this->tree[i].end = *j;
        this->tree[i].sum = 0;
    }
    int markValue = this->tree[i - 1].end + 1;
    for (; i < this->size * 2; ++i, ++markValue) {
        this->tree[i].begin = this->tree[i].end = markValue;
        this->tree[i].sum = 0;
    }
    for (i = this->size - 1; i >= 1; --i) {
        tree[i].begin = tree[i * 2].begin;
        tree[i].end = tree[i * 2 + 1].end;
        tree[i].sum = 0;
    }
}

Column &Column::operator=(Column &&other)
{
    delete tree;
    size = other.size;
    tree = other.tree;
    other.tree = nullptr;
    return *this;
}

Column::~Column()
{
    delete[] tree;
}

void Column::add(int y, int diff)
{
    int node = find(1, y);
    //cerr << "size: " << size << " node: " << node << '\n';
    assert(node >= size);
    assert(node < size * 2);
    for (; node > 0; node /= 2) {
        tree[node].sum += diff;
    }
}

size_t Column::find(size_t root, int y)
{
    if (tree[root].begin == tree[root].end) {
        assert(tree[root].begin == y);
        return root;
    } else {
        if (tree[root * 2].end >= y) {
            return find(root * 2, y);
        } else {
            return find(root * 2 + 1, y);
        }
    }
}

int Column::query(int y1, int y2)
{
    assert(y1 <= y2);
    if (tree == nullptr) {
        return 0;
    }
    return queryImpl(1, y1, y2);
}

int Column::queryImpl(size_t node, int y1, int y2)
{
    assert(node > 0);
    assert(node < size * 2);
    assert(y1 <= y2);
    if (y1 <= tree[node].begin && tree[node].end <= y2) {
        return tree[node].sum;
    } else if (tree[node].end < y1) {
        return 0;
    } else if (tree[node].begin > y2) {
        return 0;
    } else {
        return queryImpl(node * 2, y1, y2) + queryImpl(node * 2 + 1, y1, y2);
    }
}

class Table
{
    // (x,y)
public:
    // sorted by x
    template<class Iterator>
    Table(Iterator begin, Iterator end);

    ~Table();

    void add(int x, int y, int diff);

    int query(int x1, int y1, int x2, int y2);

private:
    struct Node
    {
        int begin = 0, end = 0;
        Column column;
    };

    Node *tree;
    size_t size;

    vector<int> buildTree(size_t node, vector<vector<Point<int>>>::const_iterator baseRangesBegin,
                          vector<vector<Point<int>>>::const_iterator baseRangesEnd);

    void setImpl(size_t node, int x, int y, int diff);

    int queryImpl(size_t node, int x1, int y1, int x2, int y2);
};

template<class Iterator>
Table::Table(Iterator begin, Iterator end)
{
    vector<vector<Point<int>>> baseRanges;
    baseRanges.emplace_back(1, *begin);
    for (++begin; begin != end; ++begin) {
        if (begin->x == baseRanges.back().front().x) {
            baseRanges.back().push_back(*begin);
        } else {
            baseRanges.emplace_back(1, *begin);
        }
    }
    auto width = baseRanges.size();
    this->size = mpow2(width);
    tree = new Node[size * 2];
    buildTree(1, baseRanges.begin(), baseRanges.end());
}

vector<int> Table::buildTree(size_t node, vector<vector<Point<int>>>::const_iterator baseRangesBegin,
                             vector<vector<Point<int>>>::const_iterator baseRangesEnd)
{
    auto size = distance(baseRangesBegin, baseRangesEnd);
    if (size == 0) {
        tree[node].begin = tree[node].end = tree[node - 1].end + 1;
        return {};
    } else if (size == 1) {
        auto x = (*baseRangesBegin)[0].x;
        vector<int> reducedCoord;
        for (auto p : *baseRangesBegin)
            reducedCoord.push_back(p.y);
#ifdef DEBUG
        for (int i = 1; i < reducedCoord.size(); ++i) {
            //cerr << reducedCoord[i - 1] << " " << reducedCoord[i] << '\n';
            assert(reducedCoord[i - 1] < reducedCoord[i]);
        }
#endif
        tree[node].begin = tree[node].end = x;
        tree[node].column = Column(reducedCoord.begin(), reducedCoord.end());
        return move(reducedCoord);
    } else {
        // general case
        auto partition = mpow2(size / 2);
        auto leftSideMerged = buildTree(node * 2, baseRangesBegin, baseRangesBegin + partition);
        auto rightSideMerged = buildTree(node * 2 + 1, baseRangesBegin + partition, baseRangesEnd);
        vector<int> merged;
        merge(leftSideMerged.begin(), leftSideMerged.end(), rightSideMerged.begin(), rightSideMerged.end(),
              back_inserter(merged));
        merged.erase(unique(merged.begin(), merged.end()), merged.end());
        tree[node].begin = tree[node * 2].begin;
        tree[node].end = tree[node * 2 + 1].end;
        tree[node].column = Column(merged.begin(), merged.end());
        return move(merged);
    }
}

Table::~Table()
{
    delete[] tree;
}

void Table::add(int x, int y, int diff)
{
    setImpl(1, x, y, diff);
}

void Table::setImpl(size_t node, int x, int y, int diff)
{
    assert(node > 0);
    if (node >= size * 2 || tree[node].end == 0) {
        return;
    }
    if (x < tree[node].begin || x > tree[node].end) {
        return;
    }
    tree[node].column.add(y, diff);
    setImpl(node * 2, x, y, diff);
    setImpl(node * 2 + 1, x, y, diff);
}

int Table::query(int x1, int y1, int x2, int y2)
{
    return queryImpl(1, x1, y1, x2, y2);
}

int Table::queryImpl(size_t node, int x1, int y1, int x2, int y2)
{
    //cout << "QUERY " << x1 << " " << y1 << " " << x2 << " " << y2 << '\n';
    if (x2 < tree[node].begin || tree[node].end < x1) {
        return 0;
    } else if (x1 <= tree[node].begin && tree[node].end <= x2) {
        //cout << x1 << " (" << tree[node].begin << " " << tree[node].end << ") " << x2 << "\n";
        return tree[node].column.query(y1, y2);
    }
    int result = 0;
    result += queryImpl(node * 2, x1, y1, x2, y2);
    result += queryImpl(node * 2 + 1, x1, y1, x2, y2);
    return result;
}


int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        int n;
        cin >> n;
        vector<Point<int>> input(n);
        for (auto &p : input) {
            cin >> p.x >> p.y;
        }
        sort(input.begin(), input.end());
        Table table(input.begin(), input.end());
        int q;
        cin >> q;
        unordered_map<Point<int>, int> MAP;
        while (q--) {
            string command;
            cin >> command;
            if (command[0] == 'S') {
                assert(command == "SET");
                int x, y, k;
                cin >> x >> y >> k;
                auto oldK = MAP[{x, y}];
                MAP[{x, y}] = k;
                table.add(x, y, k - oldK);
            } else {
                assert(command == "QUERY");
                int x1, y1, x2, y2;
                cin >> x1 >> y1 >> x2 >> y2;
                cout << table.query(x1, y1, x2, y2) << '\n';
            }
        }

    }
    return 0;
}