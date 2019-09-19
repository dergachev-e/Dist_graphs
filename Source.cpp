#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <set>
using namespace std;
struct Node{
    double x;
    double y;
    Node(double x_ = 0, double y_ = 0) {
        x = x_;
        y = y_;
    }
};

double dist(Node u, Node v) {
    return pow(u.x - v.x, 2) + pow(u.y - v.y, 2);
}

double len(Node u) {
    return pow(u.x, 2) + pow(u.y, 2);
}

double angle(Node u, Node v) {
    return (len(u) + len(v) - dist(u, v)) / (2.0 * sqrt(len(u))*sqrt(len(v)));
}
double angle1(Node u, Node v){
    return (len(u) + len(v) - 1.0) / (2.0 * sqrt(len(u))*sqrt(len(v)));
}

Node turn(double cos_a, double sin_a, Node v) {
    double x_ = cos_a*v.x - sin_a*v.y;
    double y_ = sin_a*v.x + cos_a*v.y;
    return Node(x_, y_);
}

vector <Node> turn60(vector <Node> U) {
    vector <Node> new_U;
    for (int i = 0; i < U.size(); ++i) {
        for (int j = 0; j < 6; ++j) {
            Node new_v = turn(0.5, sqrt(3)*0.5, U[i]);
            double flag = 0;
            for (int k = 0; k < new_U.size(); ++k) {
                if (abs(new_U[k].x - new_v.x) < 0.0000001 && abs(new_U[k].y - new_v.y) < 0.0000001)
                    flag = 1;
            }
            if (flag == 0)
                new_U.push_back(new_v);
            U[i] = new_v;
        }
    }
    return new_U;
}

vector <Node> Sum_Mink(vector <Node> U, vector <Node> V) {
    vector <Node> result;
    for (int i = 0; i < U.size(); ++i) {
        for (int j = 0; j < V.size(); ++j) {
            Node v(U[i].x + V[j].x, U[i].y + V[j].y);
            double flag = 0;
            for (int k = 0; k < result.size(); ++k) {
                if (abs(result[k].x - v.x) < 0.0000001 && abs(result[k].y - v.y) < 0.0000001)
                    flag = 1;
            }
            if (flag == 0)
                result.push_back(v);
        }
    }
    return result;
}

void Print_Dist_Graph(vector <Node> G, int number) {
    vector <pair <int, int> > Edges;
    int count = 0;
    for (int i = 0; i < G.size(); ++i) {
        for (int j = i + 1; j < G.size(); ++j) {
            if (abs(dist(G[i], G[j]) - 1) < 0.0000001) {
                count = count + 1;
                Edges.push_back(make_pair(i, j));
            }
        }
    }
    string s = to_string(number);
    cout << s;
    ofstream fout;
    fout.open("Graph" + s + ".txt");
    fout << "p edges ";
    fout << G.size() << " " << Edges.size() << endl;
    for (int i = 0; i < Edges.size(); ++i) {
        fout << "e " << Edges[i].first + 1 << " " << Edges[i].second + 1 << endl;
    }
    fout.close();
}

int main() {
    vector <Node> U1;
    U1.push_back(Node(1.0, 0.0));
    U1.push_back(turn(5.0 / 6.0, sqrt(11.0) / 6.0, Node(1.0, 0.0)));
    U1.push_back(turn(5.0 / 6.0, -sqrt(11.0) / 6.0, Node(1.0, 0.0)));
    U1.push_back(turn(sqrt(11.0) / (2.0 * sqrt(3.0)), 1.0 / (2.0 * sqrt(3.0)), Node(1.0, 0.0)));
    U1.push_back(turn(sqrt(11.0) / (2.0 * sqrt(3.0)), -1.0 / (2.0 * sqrt(3.0)), Node(1.0, 0.0)));
    U1 = turn60(U1);
    U1.push_back(Node(0.0, 0.0));
    vector <Node> V = Sum_Mink(U1, U1);
    vector <Node> V0;
    for (int i = 0; i < V.size(); ++i) {
        if (len(V[i]) < 1.0000001)
            V0.push_back(V[i]);
    }
    V0 = Sum_Mink(V0, U1);
    map < pair <long long, long long>, double> d;
    for (int i = 0; i < V0.size(); ++i) {
        for (int j = 0; j < V0.size(); ++j) {
            double cos_a = angle(V0[i], V0[j]);
            double cos_y = angle1(V0[i], V0[j]);
            double cos_b, sin_b;
            if (pow(cos_a, 2) > 1 || pow(cos_y, 2) > 1)
                continue;
            double s = V0[i].x * V0[j].y - V0[i].y * V0[j].x;
            if (s >= 0) {
                cos_b = cos_a*cos_y - sqrt(1 - pow(cos_a, 2))*sqrt(1 - pow(cos_y, 2));
                sin_b = sqrt(1 - pow(cos_b, 2));
            }
            if (s < 0) {
                if (cos_y < cos_a) {
                    cos_b = cos_a*cos_y + sqrt(1 - pow(cos_a, 2))*sqrt(1 - pow(cos_y, 2));
                    sin_b = sqrt(1 - pow(cos_b, 2));
                }
                if (cos_y >= cos_a) {
                    cos_b = cos_a*cos_y + sqrt(1 - pow(cos_a, 2))*sqrt(1 - pow(cos_y, 2));
                    sin_b = -sqrt(1 - pow(cos_b, 2));
                }
            }
            if (cos_b < sqrt(3) / 2.0)
                continue;
            cos_b = cos_b * 1000000000;
            long long cos_b1 = cos_b;
            sin_b = sin_b * 1000000000;
            long long sin_b1 = sin_b;
            if (d.find(make_pair(cos_b1, sin_b1)) == d.end()) {
                d[make_pair(cos_b1, sin_b1)] = 1;
            }
            else {
                d[make_pair(cos_b1, sin_b1)]++;
            }
        }
    }
    cout << d.size();
    set < pair< int, pair <double, double > > > angles;
    for (map < pair <long long, long long>, double>::iterator it = d.begin(); it != d.end(); ++it) {
        angles.insert(make_pair(-it->second, make_pair(it->first.first, it->first.second)));
    }
    for (set < pair< int, pair <double, double > > >::iterator it = angles.begin(); it != angles.end(); ++it) {
        cout << it->first << " " << it->second.first << " " << it->second.second << endl;
    }
    int cnt = 0;
    const int n = 100;
    vector <pair <int, int> > grid1[4*n][4*n];
    for (int i = 0; i < 4 * n; ++i) {
        for (int j = 0; j < 4 * n; j++) {
            double x1, y1, x2, y2, x3, y3, x4, y4;
            double l = j / n;

        }
    }
    for (set < pair< int, pair <double, double > > >::iterator it = angles.begin(); it != angles.end(); ++it) {
        vector <Node> G;
        vector <Node> grid[4*n][4*n];
        for (int i = 0; i < V0.size(); ++i) {
            G.push_back(V0[i]);
            int index_i = (V0[i].x + 2) * n;
            int index_j = (V0[i].y + 2) * n;
            grid[index_i][index_j].push_back(V0[i]);
        }
        for (int i = 0; i < V0.size(); ++i) {
            if (abs(V0[i].x) > 0.000001 || abs(V0[i].y) > 0.000001) {
                Node h = turn(it->second.first / 1000000000, it->second.second / 1000000000, V0[i]);
                double d = 1.0;
                for (int j = 0; j < V0.size(); ++j) {
                    d = min(d, dist(V0[j], h));
                }
                if (d > 0.000000000001) {
                    G.push_back(h);
                    int index_i = (h.x + 2) * n;
                    int index_j = (h.y + 2) * n;
                    grid[index_i][index_j].push_back(V0[i]);
                }
            }
        }
        int count = 0, index_v = 0, index_u = 0;
    /*/    for (int i = 0; i < G.size(); ++i) {
            index_v = 0;
            index_u = index_u + 1;
            for (int j = 0; j < G.size(); ++j) {
                index_v = index_v + 1;
                if (index_v < index_u) {
                    if (abs(((G[j].x - G[i].x)*(G[j].x - G[i].x) + (G[j].y - G[i].y)*(G[j].y - G[i].y)) - 1) < 0.0000001)
                        count++;
                }
            }
        }/*/
        //cout << count << " " << G.size() << endl;
        if (count >= 26750) {
            cnt++;
            cout << cnt;
            Print_Dist_Graph(G, cnt);
            ofstream lout;
            lout.open("log.txt", ios::app);
            lout << cnt << " " << it->second.first << " " << it->second.second << endl;
            lout.close();
        }
    }
}