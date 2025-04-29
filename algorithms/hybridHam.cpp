#include <bits/stdc++.h>
#include "util.hpp"

using namespace std;

// Check connectivity in the induced subgraph of unvisited vertices.
bool checkConnectivity(const vector<vector<int>>& graph, const vector<bool>& visited) {
    int n = graph.size();
    int start = -1;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            start = i;
            break;
        }
    }
    // If all vertices are visited or no unvisited vertex exists, connectivity is fine.
    if (start == -1) return true;
    
    vector<bool> seen(n, false);
    queue<int> q;
    q.push(start);
    seen[start] = true;
    int count = 0;
    while (!q.empty()){
        int u = q.front();
        q.pop();
        count++;
        for (int v : graph[u]) {
            if (!visited[v] && !seen[v]) {
                seen[v] = true;
                q.push(v);
            }
        }
    }
    int totalUnvisited = 0;
    for (int i = 0; i < n; i++) {
        if (!visited[i])
            totalUnvisited++;
    }
    return (count == totalUnvisited);
}


void extendPath(int current, const vector<vector<int>>& sortedAdj, 
                const vector<vector<int>>& graph, vector<bool>& visited, vector<int>& path) {
    int n = graph.size();
    bool extended = true;

    while (extended && int(path.size()) < n) {
        extended = false;
        // iterate over neighbours in increasing order (Va order)
        for (int cand : sortedAdj[current]) {
            if (!visited[cand]) {
                visited[cand] = true;
                // Check that the induced subgraph on remaining vertices is connected.
                if (checkConnectivity(graph, visited)) {
                    path.push_back(cand);
                    current = cand;
                    extended = true;
                    break;
                }
                // Otherwise, backtrack the visit marking.
                visited[cand] = false;
            }
        }
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


// Pósa-style rotation.
// Given a current path, try to find a pivot i such that there is an edge from path[i] to the last vertex.
// Then, the new path becomes: path[0...i] concatenated with reverse(path[i+1...end]).
// Returns true if a rotation was performed.
bool rotatePath(vector<int>& path, const vector<vector<int>>& graph) {
    int n = path.size();
    if (n < 2) return false;

    // Try all possible pivot indices from the beginning (except the very last vertex).
    for (int i = 0; i < n - 1; i++) {
        int node = path[i];
        // Check if there is an edge between path[i] and the last vertex.
        if (find(graph[node].begin(), graph[node].end(), path.back()) != graph[node].end()) {
            // Apply rotation: newPath = path[0..i] + reverse(path[i+1..end]).
            vector<int> newPath;
            for (int j = 0; j <= i; j++) {
                newPath.push_back(path[j]);
            }
            for (int j = n - 1; j > i; j--) {
                newPath.push_back(path[j]);
            }
            if(newPath != path) {  // if a different path is generated
                path = newPath;
                return true;
            }
        }
    }
    return false;  // no valid rotation found
}



// Phase 2: Given an initial path (possibly not Hamiltonian), try to extend it to a Hamiltonian path
// using rotational transformations and greedy extension.
vector<int> extendToHamiltonianPath(vector<int> path, const vector<vector<int>>& sortedAdj, 
                                    const vector<vector<int>>& graph, const vector<int>& degree) {
    int n = graph.size();
    vector<bool> visited(n, false);
    for (int v : path) {
        visited[v] = true;
    }
    
    int rotationAttempts = 0;
    const int maxRotationAttempts = n * n;
    while (int(path.size()) < n && rotationAttempts < maxRotationAttempts) {
        // Choose the end for rotation transformation.
        // We select the end with the highest degree.
        if (degree[path.front()] > degree[path.back()]) {
            // reverse the path so that highest degree vertex is at the end.
            reverse(path.begin(), path.end());
        }
        // Attempt a rotational transformation.
        bool rotated = rotatePath(path, graph);
        if (!rotated) {
            // Could not rotate: the algorithm fails to find a Hamiltonian path.
            return vector<int>(); // empty indicates failure
        }
        // Update visited: mark new vertices in path (they remain marked, but now order may change)
        fill(visited.begin(), visited.end(), false);
        for (int v : path)
            visited[v] = true;
        // Try to extend the new path using the same greedy extension as in Phase 1.
        int current = path.back();
        extendPath(current, sortedAdj, graph, visited, path);

        ++rotationAttempts;
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
    vector<vector<int>> graph(n);
    
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
    int max_deg = Vd[0];
    iota(Vd.begin(), Vd.end(), 0);
    sort(Vd.begin(), Vd.end(), [&](int a, int b){ return degree[a] > degree[b]; });

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
    }

    // ------------ Phase 2: Convert initial path into Hamiltonian path ------------
    vector<int> HamPath = extendToHamiltonianPath(best, sortedAdj, graph, degree);
    if (HamPath.empty()) {
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