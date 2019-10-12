#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <math.h>
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
int cnt = 0;

void Print_Dist_Graph(vector <Node> G, int number, set <pair <int, int> > Edges) {
    string s = to_string(number);
    cout << s << endl;
    ofstream fout;
    fout.open("Graph" + s + ".txt");
    fout << "p edges ";
    fout << G.size() << " " << Edges.size() << endl;
    for (set<pair<int, int>>::iterator it = Edges.begin(); it != Edges.end(); ++it) {
        fout << "e " << it->first + 1 << " " << it->second + 1 << endl;
    }
    fout.close();
}
const int n = 10;
set < pair< int, int> > set_grid[6 * n + 1][6 * n + 1];
int grid[6 * n + 1][6 * n + 1][1000]; // Сделать 2mn, где m - сколько сумм минковского
int grid_size[6 * n + 1][6 * n + 1];
int _count = 0;
void find_edges(int x, int y, int index, vector<pair<pair<int, int>, double>> sdvig, set <pair <int, int> > &Edges, vector <Node> G, vector < pair <int, int> > G_grid, double flag) {
    vector<pair<int, int>> neighbors;
    int i = G_grid[index].first;
    int j = G_grid[index].second;
    for (int l = 0; l < sdvig.size(); ++l) {
        int k = sdvig[l].first.first;
        int ind = sdvig[l].first.second;
        double _ind = sdvig[l].second + flag;
        if (i - k > 0 && j - ind >= 0) {
            neighbors.push_back(make_pair(i - k, j - ind));
            neighbors.push_back(make_pair(i - k - 1, j - ind));
        }
        if (i + k < 6 * n && j - ind >= 0) {
            neighbors.push_back(make_pair(i + k, j - ind));
            neighbors.push_back(make_pair(i + k + 1, j - ind));
        }
        if (i - k > 0 && j + ind < 6 * n + 1) {
            neighbors.push_back(make_pair(i - k, j + ind));
            neighbors.push_back(make_pair(i - k - 1, j + ind));
        }
        if (i + k < 6 * n  && j + ind < 6 * n + 1) {
            neighbors.push_back(make_pair(i + k, j + ind));
            neighbors.push_back(make_pair(i + k + 1, j + ind));
        }
        if (abs(ind - _ind) < 0.000001) {
            find_edges(0, 0, i, sdvig, Edges, G, G_grid, 1);
            find_edges(0, 1, i, sdvig, Edges, G, G_grid, 1);
            find_edges(1, 0, i, sdvig, Edges, G, G_grid, 1);
            find_edges(1, 1, i, sdvig, Edges, G, G_grid, 1);
        }
    }
   // cout << i << " " << j << endl;
    for (int l = 0; l < neighbors.size(); ++l) {
        for (int h = 0; h < grid_size[neighbors[l].first][neighbors[l].second]; ++h) {
            if (abs(dist(G[index], G[grid[neighbors[l].first][neighbors[l].second][h]]) - 1) < 0.000000001) {
                Edges.insert({ index, grid[neighbors[l].first][neighbors[l].second][h] });
                //cout << abs(dist(G[index], G[grid[neighbors[l].first][neighbors[l].second][h]])) << " ";
            }
        }
    }
}
void _insert_in_set(int x, int y, int k, int i, int j, int ind, double _ind) {
    _count++;
    if (_count % 1000 == 0)
        cout << _count << " ";
    if (i - k > 0 && j - ind >= 0) {
        set_grid[i][j].insert(make_pair(i - k, j - ind));
        set_grid[i][j].insert(make_pair(i - k - 1, j - ind));
    }
    if (i + k < 6 * n && j - ind >= 0) {
        set_grid[i][j].insert(make_pair(i + k, j - ind));
        set_grid[i][j].insert(make_pair(i + k + 1, j - ind));
    }
    if (i - k > 0 && j + ind < 6 * n + 1) {
        set_grid[i][j].insert(make_pair(i - k, j + ind));
        set_grid[i][j].insert(make_pair(i - k - 1, j + ind));
    }
    if (i + k < 6 * n  && j + ind < 6 * n + 1) {
        set_grid[i][j].insert(make_pair(i + k, j + ind));
        set_grid[i][j].insert(make_pair(i + k + 1, j + ind));
    }
    if (abs(ind - _ind) < 0.000001) {
        _insert_in_set(x, y, k, i, j, ind + 1, _ind);
        _insert_in_set(x, y, k, i, j, ind - 1, _ind);
    }
}
void insert_in_set(int x, int y, int k, int ind, double _ind) {
    for (int i = 0; i < 6 * n + 1; i++) {
        for (int j = 0; j < 6 * n + 1; j++) {
            _insert_in_set(x, y, k, i, j, ind, _ind);
        }
    }
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
    cout << "grid" << endl;
    vector<pair<pair<int, int>, double>> sdvig;
    double h = 1.0 / n;
    double x1, y1, x2, y2;
    int i = 0, j = 0;
    x1 = i*h;
    y1 = j*h;
    x2 = x1 + h;
    y2 = y1 + h;
    for (int k = 0; k < n; ++k) {
        double _ind = sqrt(1 / (h*h) - double(k*k));
        int ind = floor(_ind);
        sdvig.push_back({ {k,ind},_ind });
    }
    cout << "graph" << endl;
    for (set < pair< int, pair <double, double > > >::iterator it = angles.begin(); it != angles.end(); ++it) {
        vector <Node> G;
        vector < pair <int, int> > G_grid;
        for (int i = 0; i < 6 * n + 1; ++i)
            for (int j = 0; j < 6 * n + 1; ++j)
                grid_size[i][j] = 0;
        for (int i = 0; i < V0.size(); ++i) {
              G.push_back(V0[i]);
              int index_i = (V0[i].x + 3) * n;
              int index_j = (V0[i].y + 3) * n;
              grid[index_i][index_j][grid_size[index_i][index_j]] = G.size() - 1;
              grid_size[index_i][index_j]++;
              G_grid.push_back({ index_i, index_j });
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
                    int index_i = (h.x + 3) * n;
                    int index_j = (h.y + 3) * n;
                    grid[index_i][index_j][grid_size[index_i][index_j]] = G.size() - 1;
                    grid_size[index_i][index_j]++;
                    G_grid.push_back({ index_i, index_j });
                }
            }
        }
        set <pair<int, int> > Edges;
        for (int i = 0; i < G.size(); ++i) {
            int ind_i = G_grid[i].first;
            int ind_j = G_grid[i].second;
            find_edges(0, 0, i, sdvig, Edges, G, G_grid, 0);
            find_edges(0, 1, i, sdvig, Edges, G, G_grid, 0);
            find_edges(1, 0, i, sdvig, Edges, G, G_grid, 0);
            find_edges(1, 1, i, sdvig, Edges, G, G_grid, 0);
        }
        cout << Edges.size() << " ";
        if (Edges.size() >= 26750) {
            cnt++;
            Print_Dist_Graph(G, cnt, Edges);
            ofstream lout;
            lout.open("log.txt", ios::app);
            lout << cnt << " " << it->second.first << " " << it->second.second << endl;
            lout.close();
        }
        Edges.clear();
    }
}