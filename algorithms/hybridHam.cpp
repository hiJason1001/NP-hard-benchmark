#include <bits/stdc++.h>
using namespace std;

// Check connectivity in the induced subgraph of unvisited vertices.
bool checkConnectivity(const vector<vector<int>> &graph, const vector<bool> &visited) {
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

// Try to extend the current path from current vertex using greedy DFS
// (using sorted neighbor lists and feasibility connectivity check).
void extendPath(int current, const vector<vector<int>> &sortedAdj, 
                const vector<vector<int>> &graph, vector<bool> &visited, vector<int> &path) {
    int n = graph.size();
    bool extended = true;
    while (extended && path.size() < n) {
        extended = false;
        // iterate over neighbours in increasing order (Va order)
        for (int cand : sortedAdj[current]) {
            if (!visited[cand]) {
                visited[cand] = true;
                // Check that the induced subgraph on remaining vertices is connected.
                if (checkConnectivity(graph, visited)) {
                    // Candidate is acceptable.
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

// Pósa-style rotation.
// Given a current path, try to find a pivot i such that there is an edge from path[i] to the last vertex.
// Then, the new path becomes: path[0...i] concatenated with reverse(path[i+1...end]).
// Returns true if a rotation was performed.
bool rotatePath(vector<int> &path, const vector<vector<int>> &graph) {
    int n = path.size();
    int last = path.back();
    // Try all possible pivot indices from the beginning (except the very last vertex).
    for (int i = 0; i < n - 1; i++) {
        // Check if there is an edge between path[i] and the last vertex.
        if (find(graph[path[i]].begin(), graph[path[i]].end(), last) != graph[path[i]].end()) {
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

// Attempt Phase 1: build an initial path starting from a given vertex.
vector<int> buildInitialPath(int start, const vector<vector<int>> &sortedAdj, 
                             const vector<vector<int>> &graph) {
    int n = graph.size();
    vector<bool> visited(n, false);
    vector<int> path;
    int current = start;
    visited[current] = true;
    path.push_back(current);
    extendPath(current, sortedAdj, graph, visited, path);
    return path;
}

// Phase 2: Given an initial path (possibly not Hamiltonian), try to extend it to a Hamiltonian path
// using rotational transformations and greedy extension.
vector<int> extendToHamiltonianPath(vector<int> path, const vector<vector<int>> &sortedAdj, 
                                    const vector<vector<int>> &graph, const vector<int>& degree) {
    int n = graph.size();
    vector<bool> visited(n, false);
    // Mark vertices already in the path.
    for (int v : path) {
        visited[v] = true;
    }
    
    int rotationAttempts = 0;
    const int maxRotationAttempts = n * n * n;
    while (path.size() < n && rotationAttempts < maxRotationAttempts) {
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
            cout << "Phase 2 failed: no valid rotational transformation found.\n";
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

int main(){
    bool oneIndexed = true;
    int n;
    cin >> n;
    vector<vector<int>> graph(n);
    int curr = 0;
    while (1) {
        int a, b;
        cin >> a >> b;
        if (a == -1) break;
        if (oneIndexed) {
            a--;
            b--;
        }
        graph[a].push_back(b);
        graph[b].push_back(a);
    }

    vector<int> degree(n);
    for (int i = 0; i < n; i++){
        degree[i] = graph[i].size();
    }

    // Build sorted neighbour list for each vertex according to increasing degree (Va array).
    vector<vector<int>> sortedAdj(n);
    for (int i = 0; i < n; i++){
        sortedAdj[i] = graph[i];
        sort(sortedAdj[i].begin(), sortedAdj[i].end(), [&](int a, int b) {
            return degree[a] < degree[b];
        });
    }

    // Build Vd: array of vertices sorted in decreasing order of degree.
    vector<int> Vd(n);
    for (int i = 0; i < n; i++){
        Vd[i] = i;
    }
    sort(Vd.begin(), Vd.end(), [&](int a, int b){
        return degree[a] > degree[b];
    });

    // ------------ Phase 1: Create an initial path ------------
    vector<int> bestPath;
    // Try each high-degree vertex as starting vertex.
    for (int start : Vd) {
        vector<int> currPath = buildInitialPath(start, sortedAdj, graph);
        if (currPath.size() > bestPath.size())
            bestPath = currPath;
        // If a full Hamiltonian path is already found, we can break early.
        if (bestPath.size() == n)
            break;
    }
    cout << "Phase 1: Initial path of length " << bestPath.size() << "\n";
    if(bestPath.empty()){
        cout << "No initial path could be constructed.\n";
        return 0;
    }

    // If initial path already covers all vertices, proceed to Phase 2.
    if (bestPath.size() < n) {
        cout << "Initial path does not cover all vertices.\n";
        // You might choose to terminate or continue with Phase 2 attempts.
    }
    // ------------ Phase 2: Convert initial path into Hamiltonian path ------------
    vector<int> HamPath = extendToHamiltonianPath(bestPath, sortedAdj, graph, degree);
    if (HamPath.empty()){
        cout << "Algorithm failed to find a Hamiltonian path.\n";
        return 0;
    } else {
        cout << "Phase 2: Hamiltonian path found.\n";
    }
    
    return 0;
}