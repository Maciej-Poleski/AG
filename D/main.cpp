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

struct Event
{
    int x;
    int yLower, yUpper;
    bool atEnd;
};

static bool operator<(const Event &lhs, const Event &rhs)
{
    if (lhs.x < rhs.x) {
        return true;
    }
    if (lhs.x == rhs.x) {
        return !lhs.atEnd && rhs.atEnd;
    }
    return false;
}

struct Node
{
    int coveredRange;
    int coverCount;
    int beginRange;
    int endRange;
    int leftChild, rightChild;
};

static Node *tree;
static int root;

static int leftChild(int node)
{
    return tree[node].leftChild;
}

static int rightChild(int node)
{
    return tree[node].rightChild;
}

static int separationPoint(int node)
{
    assert(node >= 0);
    assert(node <= root);
    assert(tree[leftChild(node)].endRange == tree[rightChild(node)].beginRange);
    return tree[leftChild(node)].endRange;
}

static void fixCover(int node)
{
    assert(node >= 0);
    assert(node <= root);
    assert(tree[node].coverCount == 0); // Przypadek trywialny powinien zawsze być obsługiwany na miejsu
    tree[node].coveredRange = tree[leftChild(node)].coveredRange + tree[rightChild(node)].coveredRange;
}

static void insertRange(int node, int begin, int end)
{
    assert(begin < end);
    assert(node >= 0);
    assert(node <= root);
    if (tree[node].beginRange == begin && tree[node].endRange == end) {
        tree[node].coverCount += 1;
        tree[node].coveredRange = tree[node].endRange - tree[node].beginRange;
    } else {
        int separation = separationPoint(node);
        if (separation <= begin) {
            insertRange(rightChild(node), begin, end);
        } else if (end <= separation) {
            insertRange(leftChild(node), begin, end);
        } else {
            insertRange(leftChild(node), begin, separation);
            insertRange(rightChild(node), separation, end);
        }
        fixCover(node);
    }
}

static void removeRange(int node, int begin, int end)
{
    assert(begin < end);
    assert(node >= 0);
    assert(node <= root);
    if (tree[node].beginRange == begin && tree[node].endRange == end) {
        tree[node].coverCount -= 1;
        if (tree[node].coverCount == 0) {
            if (leftChild(node) != -1) {
                fixCover(node);
            } else {
                tree[node].coveredRange = 0;
            }
        }
    } else {
        int separation = separationPoint(node);
        if (separation <= begin) {
            removeRange(rightChild(node), begin, end);
        } else if (end <= separation) {
            removeRange(leftChild(node), begin, end);
        } else {
            removeRange(leftChild(node), begin, separation);
            removeRange(rightChild(node), separation, end);
        }
        if (tree[node].coverCount == 0) {
            fixCover(node);
        }
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    int z;
    cin >> z;
    while (z--) {
        vector<int> verticalPoints;
        vector<Event> events;
        int n;
        cin >> n;
        for (int i = 0; i < n; ++i) {
            int x1, y1, x2, y2;
            cin >> x1 >> y1 >> x2 >> y2;
            verticalPoints.push_back(y1);
            verticalPoints.push_back(y2);
            assert(y1 < y2);
            assert(x1 < x2);
            events.push_back({x1, y1, y2, false});
            events.push_back({x2, y1, y2, true});
        }
        sort(verticalPoints.begin(), verticalPoints.end());
        verticalPoints.erase(unique(verticalPoints.begin(), verticalPoints.end()), verticalPoints.end());
        const auto baseRangesCount = verticalPoints.size() - 1;
        tree = new Node[baseRangesCount * 2];
        vector<int> roots;
        int nextNode = 0;
        for (int i = 0; i < verticalPoints.size() - 1; ++i) {
            const auto currentNode = nextNode++;
            tree[currentNode] = {0, 0, verticalPoints[i], verticalPoints[i + 1], -1, -1};
            roots.push_back(currentNode);
        }
        while (roots.size() > 1) {
            vector<int> oldRoots;
            swap(roots, oldRoots);
            int i = 0;
            for (; i < oldRoots.size() - 1; i += 2) {
                const auto currentNode = nextNode++;
                const auto leftChild = oldRoots[i];
                const auto rightChild = oldRoots[i + 1];
                tree[currentNode] = {0, 0, tree[leftChild].beginRange, tree[rightChild].endRange, leftChild,
                                     rightChild};
                roots.push_back(currentNode);
            }
            if (i < oldRoots.size()) {
                assert(i == oldRoots.size() - 1);
                roots.push_back(oldRoots[i]);
            }
        }
        root = roots[0];
        roots.clear();
        sort(events.begin(), events.end());
        long area = 0;
        int lastX = events[0].x;
        for (auto event : events) {
            const auto coveredRange = tree[root].coveredRange;
            const auto currentX = event.x;
            area += (currentX - lastX) * coveredRange;
            lastX = currentX;
            if (!event.atEnd) {
                insertRange(root, event.yLower, event.yUpper);
            } else {
                removeRange(root, event.yLower, event.yUpper);
            }
        }
        cout << area << '\n';
        delete[] tree;
    }
    return 0;
}