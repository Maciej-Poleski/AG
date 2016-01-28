#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <unordered_map>

using namespace std;

struct Point
{
    int x, y;
};

vector<Point> readInput(const char *inputFile)
{
    ifstream input(inputFile);
    int n;
    input >> n;
    if (n != 1) {
        cerr << "Zła ilość zestawów testowych: " << n << " (a oczekiwano 1)\n";
        exit(2);
    }
    input >> n;
    vector<Point> result(n);
    for (auto &p:result) {
        input >> p.x >> p.y;
    }
    return result;
}

vector<int[3]> readOutput(const char *outputFile)
{
    ifstream output(outputFile);
    int n;
    output >> n;
    vector<int[3]> result(n);
    for (auto &t:result) {
        output >> t[0] >> t[1] >> t[2];
    }
    return result;
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

static bool isCCW(Point p1, Point p2, Point p3)
{
    return turnsLeft(p1, p2, p3) && turnsLeft(p2, p3, p1) && turnsLeft(p3, p1, p2);
}


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

int main(int argc, char **argv)
{
    ios_base::sync_with_stdio(false);
    if (argc != 3) {
        cerr << argv[0] << " [input] [output]\n";
        return 1;
    }
    vector<Point> input = readInput(argv[1]);
    vector<int[3]> output = readOutput(argv[2]);

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
}
