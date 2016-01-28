#define DEBUG
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

/******************* ROZWIĄZANIE *****************/

static bool operator<(const Point<long> &lhs, const Point<long> &rhs)
{
    return lhs.y > rhs.y || (lhs.y == rhs.y && lhs.x < rhs.x);
}

static bool isCCW(Point<long> p1, Point<long> p2, Point<long> p3)
{
    return turnsLeft<long>(p1, p2, p3) && turnsLeft<long>(p2, p3, p1) && turnsLeft<long>(p3, p1, p2);
}


/****************** WERYFIKATOR ***********************/

struct Edge : pair<int, int>
{
    Edge(int a, int b) : pair<int, int>(a, b)
    {
        if (a > b) {
            using std::swap;
            swap(first, second);
        }
    }
};

namespace std {
    template<>
    struct hash<Edge>
    {
        size_t operator()(const Edge &o) const noexcept
        {
            return static_cast<size_t>(o.first) << 32 | o.second;
        }
    };
}

static vector<array<int, 3>> output;

static void printCCW(const vector<Point<long>> &input, int p1, int p2, int p3)
{
    if (isCCW(input[p1], input[p2], input[p3])) {
#ifdef DEBUG
        output.push_back({p1, p2, p3});
#endif
        cout << p1 << ' ' << p2 << ' ' << p3 << '\n';
    } else {
        assert(isCCW(input[p1], input[p3], input[p2]));
#ifdef DEBUG
        output.push_back({p1, p3, p2});
#endif
        cout << p1 << ' ' << p3 << ' ' << p2 << '\n';
    }
}

#ifndef DEBUG
#define invokeVerification(x)
#else

static void invokeVerification(const vector<Point<long>> &input)
{
    // sprawdź ilość
    assert(input.size() == output.size() + 2);

    // sprawdź CCW
    for (const auto &t : output) {
        assert(isCCW(input[t[0]], input[t[1]], input[t[2]]));
    }

    unordered_map<Edge, vector<int>> edgesToPoints;
    for (size_t i = 0; i < output.size(); ++i) {
        edgesToPoints[{output[i][0], output[i][1]}].push_back(i);
        edgesToPoints[{output[i][1], output[i][2]}].push_back(i);
        edgesToPoints[{output[i][2], output[i][0]}].push_back(i);
    }

    size_t singleEdgesCount = 0;
    for (const auto &edge : edgesToPoints) {
        const int n = input.size();
        if (edge.second.size() == 1) {
            // oczekiwana krawędź (zewnętrzna) wielokąta
            assert((edge.first.first == (edge.first.second + 1) % n) ||
                   ((edge.first.first + 1) % n == edge.first.second));
            singleEdgesCount += 1;
        } else {
            // krawędź wewnątrzna ma zawsze 2 sąsiednie trójkąty
            assert(edge.second.size() == 2);
        }
    }
    // sprawdź czy krawędzi zewnętrznych jest tyle ile trzeba
    assert(singleEdgesCount == input.size());

    // reset
    output.clear();
}

#endif

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
        invokeVerification(input);
    }
    return 0;
}