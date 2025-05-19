#include <bits/stdc++.h>
#include "util.hpp"

using namespace std;



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
