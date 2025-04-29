#include <bits/stdc++.h>
#include "util.hpp"

using namespace std;

int get_random_idx(int size) {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, size - 1);
    return dist(rng);
}

void extendPath(int start, const vector<vector<int>>& sortedAdj, 
                const vector<vector<int>>& graph, vector<bool>& visited, vector<int>& path) {
    int n = graph.size();
    vector<int> deg(n);
    for (int i = 0; i < n; i++) {
        deg[i] = graph[i].size();
    }

    int prev_path_size = 0;
    while (int(path.size()) < n && int(path.size()) > prev_path_size) {
        for (int neigh : sortedAdj[start]) {
            if (visited[neigh]) continue;

            int mutual_count = 0;
            for (int mutual : graph[neigh]) {
                if (visited[mutual]) continue;
                
                mutual_count++;
            }
            
            bool does_it_create_unreachable_vertex = false;
            for (int mutual : graph[neigh]) {
                if (visited[mutual]) continue;

                if (deg[mutual] == 1 && mutual_count > 1) {
                    does_it_create_unreachable_vertex = true;
                    break;
                }
            }

            if (!does_it_create_unreachable_vertex) {
                path.push_back(neigh);
                visited[neigh] = true;
                start = neigh;
                for (int mutual : graph[start]) {
                    deg[mutual]--;
                }
                break;
            }
        }
        prev_path_size = path.size();
    }
}

// Phase 1
vector<int> buildInitialPath(int start, const vector<vector<int>>& sortedAdj, 
                             const vector<vector<int>>& graph) {
        
    int n = graph.size();
    vector<bool> visited(n, false);
    vector<int> path;
    path.reserve(n);
    path.push_back(start);
    visited[start] = true;
    extendPath(start, sortedAdj, graph, visited, path);

    return path;
}


// Given a current path, try to find a pivot i such that there is an edge from path[i] to the last vertex.
// Then, the new path becomes: path[0...i] concatenated with reverse(path[i+1...end]).
// Returns true if a rotation was performed.
void rotatePath(vector<int>& path, const vector<vector<int>>& graph) {
    int n = path.size();
    if (n <= 2) return;

    vector<unordered_set<int>> fast_lookup_graph(graph.size());
    for (int i = 0; i < int(graph.size()); i++) {
        for (int neigh : graph[i]) {
            fast_lookup_graph[i].insert(neigh);
        }
    }

    vector<int> valid_rotation_points;

    // Try all possible pivot indices from the beginning (except the last two vertices).
    for (int i = 0; i < n - 2; i++) {
        int node = path[i];
        // Check if there is an edge between path[i] and the last vertex.
        if (fast_lookup_graph[node].count(path.back())) {
            valid_rotation_points.push_back(node);
        }
    }

    if (valid_rotation_points.empty()) return;

    // paper does not say which node to use as pivot
    // so I'm picking a random one
    int rand_idx = get_random_idx(valid_rotation_points.size());
    int pivot = 0;
    for (int i = 0; i < int(path.size()); i++) {
        if (path[i] == valid_rotation_points[rand_idx]) {
            pivot = i;
            break;
        }
    }

    // Apply rotation: newPath = path[0..i] + reverse(path[i+1..end]).
    vector<int> new_path;
    for (int j = 0; j <= pivot; j++) {
        new_path.push_back(path[j]);
    }
    for (int j = n - 1; j > pivot; j--) {
        new_path.push_back(path[j]);
    }

    path = new_path;
}


vector<int> extendToHamiltonianPath(vector<int> path, const vector<vector<int>>& sortedAdj, 
                                    const vector<vector<int>>& graph, const vector<int>& degree) {
    int n = graph.size();
    vector<bool> visited(n, false);
    for (int v : path) {
        visited[v] = true;
    }
    
    int rotationAttempts = 0;
    const int max_rotation_attempts = n * n;
    while (int(path.size()) < n && rotationAttempts < max_rotation_attempts) {
        vector<int> old_path = path;

        if (degree[path.front()] > degree[path.back()]) {
            reverse(path.begin(), path.end());
        }

        rotatePath(path, graph);
        
        fill(visited.begin(), visited.end(), false);
        for (int v : path)
            visited[v] = true;
            
        extendPath(path.back(), sortedAdj, graph, visited, path);

        ++rotationAttempts;

        if (path == old_path) break;
        // for (auto it:path)cout<<it<<" ";
        // cout<<endl;
    }
    return path;
}


int main(int argc, char* argv[]) {
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string INPUT = argv[1];
    int n, m;
    vector<vector<int>> graph;
    
    if (!Util::get_adjList(INPUT, n, m, graph)) {
        return 0;
    }

    cout << "Number of Vertices: " << n << ", number of edges: " << m << endl;

    auto start = chrono::steady_clock::now();


    vector<int> degree(n);
    for (int i = 0; i < n; ++i)
        degree[i] = graph[i].size();

    vector<vector<int>> sortedAdj(n);
    for (int i = 0; i < n; ++i) {
        sortedAdj[i] = graph[i];
        sort(sortedAdj[i].begin(), sortedAdj[i].end(),
             [&](int a, int b){ return degree[a] < degree[b]; });
    }

    vector<int> Vd(n);
    for (int i = 0; i < n; i++) Vd[i] = i;
    sort(Vd.begin(), Vd.end(), [&](int a, int b){ return degree[a] > degree[b]; });
    int max_deg = degree[Vd[0]];


    vector<int> best;
    for (int u : Vd) {
        if (degree[u] != max_deg) break;
        auto p = buildInitialPath(u, sortedAdj, graph);

        if (p.size() > best.size()) best = move(p);
        if ((int)best.size() == n) break;
    }

    if (best.empty()) {
        auto end = chrono::steady_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
        Util::display_time(duration, Util::micro);
        cout << "No" << endl;
        return 0;
    } else if (int(best.size()) == n) {
        auto end = chrono::steady_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
        Util::display_time(duration, Util::micro);
        cout << "Yes" << endl;
        return 0;
    }

    // for (auto it : best) cout<<it<<" ";
    // cout<<endl;

    // ------------ Phase 2: Convert initial path into Hamiltonian path ------------
    vector<int> HamPath = extendToHamiltonianPath(best, sortedAdj, graph, degree);
    if (int(HamPath.size()) != n) {
        auto end = chrono::steady_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
        Util::display_time(duration, Util::micro);
        cout << "No" << endl;
    } else {
        auto end = chrono::steady_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
        Util::display_time(duration, Util::micro);
        cout << "Yes" << endl;
    }
    
    return 0;
}