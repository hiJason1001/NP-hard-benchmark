#include <bits/stdc++.h>
#include "util.hpp"

using namespace std;

bool Hamiltonian_path(vector<vector<int> >& adj, int N)
{
    vector<vector<bool>> dp(N, vector<bool>(1 << N)); 
 
    for (int i = 0; i < N; i++)
        dp[i][(1 << i)] = true;

    for (int i = 0; i < (1 << N); i++) {
 
        for (int j = 0; j < N; j++) {
 
            if (i & (1 << j)) {
 
                for (int k = 0; k < N; k++) {
 
                    if (i & (1 << k)
                        && adj[k][j]
                        && j != k
                        && dp[k][i ^ (1 << j)]) {

                        dp[j][i] = true;
                        break;
                    }
                }
            }
        }   
    }
 
    for (int i = 0; i < N; i++) {
        if (dp[i][(1 << N) - 1])
            return true;
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
