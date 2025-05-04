#include <bits/stdc++.h>
#include "util.hpp"

using namespace std;

vector<unordered_map<int, bool>> get_random_coloring_graph(const vector<vector<int>>& graph, int n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 1);

    vector<unordered_map<int, bool>> res(n);
    for (int i = 0; i < n; i++) {
        for (int node : graph[i]) {
            bool has_color = dist(gen);
            res[i][node] = has_color;
            res[node][i] = has_color;
        }
    }

    return res;
}

bool HPA2() {
    // cout<<"XD\n";
    return false;
}


bool HPA1(const vector<vector<int>>& graph, int n, int start, int finish) {
    assert(0 <= start && start < n);
    assert(0 <= finish && finish < n);
    auto red_graph = get_random_coloring_graph(graph, n);
    list<int> even{start};
    list<int> odd{finish};
    list<int> even_outback;
    list<int> odd_outback;
    for (int i = 0; i < n; i++) {
        if (i == start || i == finish) continue;

        if (i % 2 == 0) even_outback.push_back(i);
        else odd_outback.push_back(i);
    }

    const int LIM = sqrt(n);

    // stage 1
    while (!even_outback.empty()) {
        int last_even = even.back();
        bool extended = false;
        for (auto it = even_outback.begin(); it != even_outback.end(); ++it) {
            int node = *it;
            if (red_graph[last_even][node]) {
                even.push_back(node);
                even_outback.erase(it);
                extended = true;
                break;
            }
        }
        if (!extended) break;
    }
    if (int(even_outback.size()) >= LIM) return HPA2();

    while (!odd_outback.empty()) {
        int first_odd = odd.front();
        bool extended = false;
        for (auto it = odd_outback.begin(); it != odd_outback.end(); ++it) {
            int node = *it;
            if (red_graph[first_odd][node]) {
                odd.push_front(node);
                odd_outback.erase(it);
                extended = true;
                break;
            }
        }
        if (!extended) break;
    }
    if (int(odd_outback.size()) >= LIM) return HPA2();

    // stage 2
    vector<int> even_arr(even.begin(), even.end());
    vector<int> odd_arr(odd.begin(), odd.end());

    vector<tuple<int, int, int>> pairs;
    for (int i = 0; i < int(even_arr.size()); ++i) {
        for (int j = 0; j < int(odd_arr.size()); ++j) {
            pairs.emplace_back(i + j, i, j);  // sorting by (i + j), then (i, j)
        }
    }
    sort(pairs.begin(), pairs.end());
    
    list<int> P;
    for (auto [sum, i, j] : pairs) {
        if (sum >= LIM) {
            return HPA2();
        }
        int node1 = even_arr[even_arr.size() - 1 - i];
        int node2 = odd_arr[j];
        if (red_graph[node1][node2]) {
            for (int k = 0; k <= int(even_arr.size()) - 1 - i; k++) {
                P.push_back(even_arr[k]);
            }
            for (int k = int(even_arr.size()) - i; k < int(even_arr.size()); k++) {
                even_outback.push_back(even_arr[k]);
            }
            for (int k = 0; k < j; k++) {
                odd_outback.push_back(odd_arr[k]);
            }
            for (int k = j; k < int(odd_arr.size()); k++) {
                P.push_back(odd_arr[k]);
            }
            break;
        }
    }

    if (int(P.size()) == n) return true;

    // stage 3
    int idx = P.size() - 1;
    auto iter = P.end();
    iter--;
    int x = *iter;
    while (!even_outback.empty() && idx > 3) {
        if (iter == P.begin()) break;

        int v = even_outback.front();
        int PredX = *prev(iter);
        if (red_graph[PredX][x] && red_graph[v][x]) {
            auto prev_iter = prev(iter);
            P.insert(prev_iter, v);
            even_outback.pop_front();
            iter--;
            x = *iter;
            idx--;
        }
        else if (red_graph[PredX][v]) {
            iter--;
            x = *iter;
            idx--;
        }
        else {
            iter--;
            iter--;
            x = *iter;
            idx -= 2;
        }
    }
    if (!even_outback.empty()) return HPA2();

    idx = P.size() - 1;
    iter = P.end();
    iter--;
    x = *iter;
    while (!odd_outback.empty() && idx > 3) {
        int v = odd_outback.front();
        int PredX = *prev(iter);
        if (red_graph[PredX][x] && red_graph[v][x]) {
            auto prev_iter = prev(iter);
            P.insert(prev_iter, v);
            odd_outback.pop_front();
        }
        else if (red_graph[PredX][v]) {
            iter--;
            x = *iter;
            idx--;
        }
        else {
            iter--;
            iter--;
            x = *iter;
            idx -= 2;
        }
    }
    if (!odd_outback.empty()) return HPA2();

    return int(P.size()) == n;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string INPUT = argv[1];
    int n, m;
    vector<vector<int>> graph;

    if (!Util::get_adjList(INPUT, n, m, graph)) {
        return 1;
    }
    cout << "Number of Vertices: " << n << ", number of edges: " << m << '\n';


    
    auto start = chrono::steady_clock::now();

    bool found = HPA1(graph, n, 0, n - 1);

    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
    Util::display_time(duration, Util::micro);

    cout << (found ? "Yes" : "No") << endl;

    return 0;
}
