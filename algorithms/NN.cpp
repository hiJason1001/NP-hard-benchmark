#include <bits/stdc++.h>
// #include "util.hpp"
#include "util.hpp"

using namespace std;

bool Hamiltonian_path(const vector<vector<int>>& adjMatrix, int n) {
    vector<int> best_path;

    for (int start = 0; start < n; ++start) {
        vector<bool> visited(n, false);
        vector<int> path;
        int current = start;

        visited[current] = true;
        path.push_back(current);

        for (int step = 1; step < n; ++step) {
            int next = -1;
            for (int j = 0; j < n; ++j) {
                if (!visited[j] && adjMatrix[current][j] == 1) {
                    next = j;
                    break;
                }
            }

            if (next == -1) break;
            visited[next] = true;
            path.push_back(next);
            current = next;
        }

        if (path.size() == n) {
            return true;
        }
    }

    return false;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string INPUT = argv[1];
    int n, m;
    vector<vector<int>> adjMatrix;

    if (!Util::get_adjMatrix(INPUT, n, m, adjMatrix)) {
        return 1;
    }
    cout << "Number of Vertices: " << n << ", number of edges: " << m << '\n';


    
    auto start = chrono::steady_clock::now();

    bool found = Hamiltonian_path(adjMatrix, n);

    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
    Util::display_time(duration, Util::micro);

    cout << (found ? "Yes" : "No") << endl;

    return 0;
}
