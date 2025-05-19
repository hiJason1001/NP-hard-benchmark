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


bool HPA3() {
    return false;
}

static std::vector<bool> make_in_env(const vector<vector<int>>& F, int n) {
    std::vector<bool> inEnv(n, false);
    for (auto &path : F)
        for (int v : path)
            inEnv[v] = true;
    return inEnv;
}

static void peel_high_fan(std::unordered_set<int>& X,
                          std::deque<int>& List,
                          const vector<vector<int>>& graph)
{
    while (true) {
        int sz = X.size();
        vector<int> toRemove;
        for (int v : X) {
            if ((int)graph[v].size() >= 3 * sz) {
                toRemove.push_back(v);
            }
        }
        if (toRemove.empty()) break;
        for (int v : toRemove) {
            X.erase(v);
            List.push_front(v);
        }
    }
}


static bool extend_envelope_newpath(vector<vector<int>>& F,
                                    int v,
                                    const vector<vector<int>>& graph,
                                    const vector<bool>& inEnv)
{
    // collect all endpoints of current envelope
    std::vector<int> endpoints;
    for (auto &path : F) {
        endpoints.push_back(path.front());
        endpoints.push_back(path.back());
    }

    // try all distinct pairs (u,w) in graph[v] that lie in endpoints
    for (int u : graph[v]) if (inEnv[u]) {
        for (int w : graph[v]) if (w != u && inEnv[w]) {
            // insert a new path [u, v, w]
            F.push_back({u, v, w});
            return true;
        }
    }
    return false;
}

static bool extend_envelope_endpoint(vector<vector<int>>& F,
                                     int v,
                                     const vector<vector<int>>& graph,
                                     const vector<bool> &inEnv)
{
    for (auto &path : F) {
        bool atFront = (path.front() == v);
        bool atBack  = (path.back()  == v);
        if (!atFront && !atBack) continue;
        for (int u : graph[v]) {
            if (!inEnv[u]) {
                if (atFront) path.insert(path.begin(), u);
                else         path.push_back(u);
                return true;
            }
        }
    }
    return false;
}

bool envelope_contains(const vector<vector<int>>& F, int v, int start, int finish) {
    for (auto &path : F) {
        if (v == start && path.front() == start) return true;
        if (v == finish && path.back() == finish) return true;
        // internal vertex check
        for (size_t i = 1; i+1 < path.size(); ++i) {
            if (path[i] == v) return true;
        }
    }
    return false;
}


bool find_envelope(const unordered_set<int>& X,
                   const vector<vector<int>>& graph,
                   int start, int finish,
                   vector<vector<int>>& F) {
    // Convert X to vector for permutation
    vector<int> Xv(X.begin(), X.end());
    int m = Xv.size();
    if (m == 0) return true;
    if (m > 8) return false;
    // If both start and finish are in X, ensure distinct paths
    bool startIn = X.count(start), finishIn = X.count(finish);
    sort(Xv.begin(), Xv.end());
    do {
        // Enforce start/finish positions
        if (startIn && Xv[0] != start) continue;
        if (finishIn && Xv.back() != finish) continue;
        // Build envelope paths by splitting on non-edges
        vector<vector<int>> cand;
        cand.push_back({Xv[0]});
        for (int i = 1; i < m; ++i) {
            int u = Xv[i-1], v = Xv[i];
            if (find(graph[u].begin(), graph[u].end(), v) != graph[u].end()) {
                // same path
                cand.back().push_back(v);
            } else {
                // new path
                cand.push_back({v});
            }
        }
        // Verify envelope properties a–d
        // (a) covers all X by construction
        // (b) internal vertices are from X: holds
        // (c) start/finish endpoints
        if (startIn) {
            bool ok=false;
            for (auto &P:cand) if (P.front()==start) ok=true;
            if (!ok) continue;
        }
        if (finishIn) {
            bool ok=false;
            for (auto &P:cand) if (P.back()==finish) ok=true;
            if (!ok) continue;
        }
        if (startIn && finishIn && cand.size()>1) {
            // ensure start and finish on different paths
            int pi=-1, pj=-1;
            for (int i = 0; i < int(cand.size()); ++i) {
                if (cand[i].front()==start) pi=i;
                if (cand[i].back()==finish) pj=i;
            }
            if (pi==-1||pj==-1||pi==pj) continue;
        }
        // (d) for each edge in paths, at least one end in X (holds)
        // Accept this envelope
        F = cand;
        return true;
    } while (next_permutation(Xv.begin(), Xv.end()));
    return false;
}


bool stage2_build_envelope(int start, int finish,
                           const vector<vector<int>>& graph,
                           std::unordered_set<int> &T,
                           vector<vector<int>>& F)
{
    std::unordered_set<int> X = T;
    std::deque<int> List;
    peel_high_fan(X, List, graph);

    if (!find_envelope(X, graph, start, finish, F))
        return false;

    while (!List.empty()) {
        int v = List.front();
        List.pop_front();

        auto inEnv = make_in_env(F, graph.size());

        if (inEnv[v]) {
            X.insert(v);
            continue;
        }

        bool inserted = false;

        if (!inserted && extend_envelope_endpoint(F, v, graph, inEnv)) {
            inserted = true;
        }

        if (!inserted && extend_envelope_newpath(F, v, graph, inEnv)) {
            inserted = true;
        }

        if (!inserted) {
            return false;
        }

        auto afterEnv = make_in_env(F, graph.size());
        assert(afterEnv[v]);

        int countV = 0;
        for (auto &path : F)
            for (int x : path)
                if (x == v) ++countV;
        assert(countV == 1);
        X.insert(v);
    }


    return true;
}

static bool extend_chain_backtrack(
    vector<int>& P,
    const vector<vector<int>>& orange,
    vector<bool>& used,
    int rn,
    bool back)
{
    struct Frame { int v, depth, idx; };
    vector<Frame> stk;
    vector<int> trail;
    stk.push_back({ back ? P.back() : P.front(), 0, 0 });

    while (!stk.empty()) {
        auto [v, depth, idx] = stk.back(); stk.pop_back();

        auto &nbrs = orange[v];
        bool advanced = false;
        for (int i = idx; i < (int)nbrs.size(); ++i) {
            int y = nbrs[i];
            if (used[y]) continue;
            // choose y
            used[y] = true;
            trail.push_back(y);
            // if full length, commit:
            if (depth+1 == rn) {
                if (back) {
                    for (int z : trail) P.push_back(z);
                } else {
                    for (int j = (int)trail.size()-1; j >= 0; --j)
                        P.insert(P.begin(), trail[j]);
                }
                return true;
            }
            // push resume and descend
            stk.push_back({ v, depth, i+1 });
            stk.push_back({ y, depth+1, 0 });
            advanced = true;
            break;
        }
        if (!advanced) {
            // backtrack one choice
            if (!trail.empty()) {
                used[trail.back()] = false;
                trail.pop_back();
            }
        }
    }
    // clear any leftover trail marks
    for (int y : trail) used[y] = false;
    return false;
}

bool stage3_extend_orange(int n,
    const std::unordered_set<int>& T,
    std::vector<std::vector<int>>& F,
    const std::vector<std::vector<int>>& orange)
{
    int k = F.size();
    int t = T.size();
    int limit = n/12 - 3*t;
    int denom = 2*(k-1);
    int rn = (denom > 0 && limit > 0) ? (limit / denom) : 0;

    std::vector<bool> used(n,false);
    for (auto &P : F)
        for (int v : P)
            used[v] = true;

    for (int i = 0; i+1 < k; ++i) {
        // extend back of F[i]
        if (!extend_chain_backtrack(F[i], orange, used, rn, true))
            return false;

        // post-back extension check
        int szb = F[i].size();
        for (int d = 1; d <= rn; ++d) {
            int y = F[i][szb-d];
            assert(std::find(orange[F[i][szb-d-1]].begin(),
                             orange[F[i][szb-d-1]].end(), y)
                   != orange[F[i][szb-d-1]].end());
        }

        // extend front of F[i+1]
        if (!extend_chain_backtrack(F[i+1], orange, used, rn, false))
            return false;

        // post-front extension check
        for (int d = 0; d < rn; ++d) {
            int y = F[i+1][d];
            assert(std::find(orange[F[i+1][d+1]].begin(),
                             orange[F[i+1][d+1]].end(), y)
                   != orange[F[i+1][d+1]].end());
        }

        // global no-duplicate check
        {
          std::unordered_set<int> seen;
          for (auto &P : F)
            for (int v : P) {
              if (!used[v]) return false;
              if (!seen.insert(v).second) return false;
            }
        }
    }
    return true;
}


// Stage 4: Sew F-paths into one path Q0 using yellow edges
bool stage4_sew_yellow(
    vector<vector<int>>& F,
    vector<int>& Q0,
    const vector<vector<int>>& yellow)
{
    int k = F.size();
    Q0.clear();
    // for each i from 0..k-2, find xi in last 10% of F[i], yi in first 10% of F[i+1]
    // such that yellow[xi][yi], then splice
    for (int i = 0; i < k-1; ++i) {
        auto &A = F[i];
        auto &B = F[i+1];
        int a_start = (9*A.size())/10;
        int b_end   = B.size()/10;
        bool linked = false;
        for (int ia = a_start; ia < (int)A.size(); ++ia) {
            for (int ib = 0; ib < b_end; ++ib) {
                int x = A[ia], y = B[ib];
                // check yellow edge
                for (int nbr : yellow[x]) if (nbr==y) {
                    // splice: trim A after ia and B before ib, then join
                    Q0.insert(Q0.end(), A.begin(), A.begin()+ia+1);
                    Q0.insert(Q0.end(), B.begin()+ib, B.end());
                    linked = true;
                    break;
                }
            }
            if (linked) break;
        }
        if (!linked) return false;  // fail ⇒ HPA3
        // prepare for next iteration:
        A = Q0;  // merge into single growing path
        Q0.clear();
    }
    // the last merged result lives in A = F[k-1]:
    Q0 = F.back();
    return true;
}

bool stage5_partition(int n,
    const std::vector<int>& Q0,
    std::vector<std::vector<int>>& Qs,
    const std::vector<std::vector<int>>& yellow)
{
    // 1) Mark Q0
    std::vector<bool> inQ0(n,false);
    for (int v : Q0) inQ0[v] = true;

    // 2) used[v]=true if v in Q0 or already in some Qi
    std::vector<bool> used = inQ0;

    // threshold = |Q0|/12
    int thresh = Q0.size() / 12;

    Qs.clear();
    std::deque<int> X;  // use deque for O(1) front inserts

    for (int start = 0; start < n; ++start) {
        if (used[start]) continue;
        // begin new path
        X.clear();
        X.push_back(start);
        used[start] = true;

        // Forward extension
        while (true) {
            int x = X.back();
            bool ext = false;
            for (int y : yellow[x]) {
                if (!used[y]) {
                    X.push_back(y);
                    used[y] = true;
                    ext = true;
                    break;
                }
            }
            if (!ext) break;
        }

        // Snapshot used for rollback
        auto used_snapshot = used;

        // Backward extension
        while (true) {
            int x = X.front();
            bool ext = false;
            for (int y : yellow[x]) {
                if (!used[y]) {
                    X.push_front(y);
                    used[y] = true;
                    ext = true;
                    break;
                }
            }
            if (!ext) break;
        }

        // Check anchors into Q0
        auto count_conn = [&](int u) {
            int c = 0;
            for (int w : yellow[u])
                if (inQ0[w]) ++c;
            return c;
        };

        if (count_conn(X.front()) >= thresh &&
            count_conn(X.back())  >= thresh)
        {
            // accept
        }
        else {
            // rollback backward marks if needed
            // we require a circuit: find w in X.back()’s neighbors back to X.front()
            used = used_snapshot;

            int u = X.front();
            // build set of X for O(1) test
            std::unordered_set<int> Xset(X.begin(), X.end());

            int best_idx = -1;
            // scan neighbors of u once
            for (int w : yellow[u]) {
                if (Xset.count(w)) {
                    // position of w in X:
                    auto it = std::find(X.begin(), X.end(), w);
                    best_idx = std::distance(X.begin(), it);
                    break;
                }
            }
            if (best_idx < 0) return false;

            // truncate and close circuit by re-appending u
            X.resize(best_idx+1);
            X.push_back(u);
        }

        // Ensure no duplicate within X
        {
          std::unordered_set<int> seen;
          for (int v : X)
            assert(seen.insert(v).second);
        }

        // Move X into Qs
        Qs.emplace_back(X.begin(), X.end());
    }
    return true;
}


// Stage 6: Merge Q0, Q1..Qm into one Hamiltonian path using green and yellow edges
bool stage6_merge(int n,
    const unordered_set<int>& T,
    vector<int>& R,                     // accumulates the merged path (initially Q0)
    vector<vector<int>>& Qs,
    const vector<vector<int>>& green,
    const vector<vector<int>>& yellow) 
{
    auto is_green = [&](int u, int v) {
        return find(green[u].begin(), green[u].end(), v) != green[u].end();
    };
    auto is_yellow = [&](int u, int v) {
        return find(yellow[u].begin(), yellow[u].end(), v) != yellow[u].end();
    };

    // We will merge each Qi into R in turn
    for (auto &Q : Qs) {
        bool merged = false;

        // --- Case 3: Q is a circuit (first == last) ---
        if (!Q.empty() && Q.front() == Q.back()) {
            int L = (int)Q.size() - 1;  // drop duplicate last
            // Try every break point in circuit Q
            for (int j = 0; j < L && !merged; ++j) {
            // Build Qpath = Q[j..L-1] + Q[0..j]
                vector<int> Qpath;
                for (int x = j; x < L; ++x) Qpath.push_back(Q[x]);
                for (int x = 0; x < j; ++x) Qpath.push_back(Q[x]);

                // Now try to splice this Qpath into R
                int Rsz = R.size();
                for (int i = 0; i + 1 < Rsz && !merged; ++i) {
                    int u = R[i], up = R[i+1];
                    if (T.count(u) || T.count(up)) continue;
                    int v = Qpath.front(), vp = Qpath[1];
                    if (is_green(u, v) && is_green(up, vp)) {
                        // splice: [Start..u] + Qpath + [up..Finish]
                        vector<int> newR;
                        newR.insert(newR.end(), R.begin(), R.begin()+i+1);
                        newR.insert(newR.end(), Qpath.begin(), Qpath.end());
                        newR.insert(newR.end(), R.begin()+i+1, R.end());
                        R = move(newR);
                        merged = true;
                    }
                }
            }
            if (!merged) return false;
            continue;
        }

        // --- Case 4: Q is an open path ---
        // Build A = {w in R\T | successor w' in R is in T, and yellow edge (w, Q.front())}
        //       B = {w in R\T | predecessor w' in R is in T, and yellow edge (Q.back(), w)}
        vector<int> A, B;
        int Rsz = R.size();
        for (int i = 0; i + 1 < Rsz; ++i) {
            int w = R[i], wsucc = R[i+1];
            if (!T.count(w) && T.count(wsucc) && is_yellow(w, Q.front()))
                A.push_back(w);
        }
        for (int i = 1; i < Rsz; ++i) {
            int w = R[i], wpred = R[i-1];
            if (!T.count(w) && T.count(wpred) && is_yellow(Q.back(), w))
                B.push_back(w);
        }
        if (A.empty() || B.empty()) return false;

        // Find I: the shortest prefix of R that contains >= |A|/2 of A or >= |B|/2 of B
        int halfA = (A.size()+1)/2, halfB = (B.size()+1)/2;
        vector<bool> inA(Rsz,false), inB(Rsz,false);
        for (int i = 0; i < Rsz; ++i) {
            if (find(A.begin(), A.end(), R[i])!=A.end()) inA[i]=true;
            if (find(B.begin(), B.end(), R[i])!=B.end()) inB[i]=true;
        }
        int cntA=0, cntB=0, Iend=0;
        for (; Iend < Rsz; ++Iend) {
            if (inA[Iend]) cntA++;
            if (inB[Iend]) cntB++;
            if (cntA >= halfA || cntB >= halfB) break;
        }
        if (cntA < halfA) {
            // reverse Q and swap A/B sets
            reverse(Q.begin(), Q.end());
            swap(A,B);
            // recompute half sizes and index sets
            halfA = (A.size()+1)/2;
            inA.assign(Rsz,false);
            for (int i=0;i<Rsz;++i)
                if (find(A.begin(),A.end(),R[i])!=A.end()) 
                    inA[i]=true;
            // recompute Iend
            cntA=0;
            for (Iend=0;Iend<Rsz;++Iend) {
                if (inA[Iend]) cntA++;
                if (cntA >= halfA) break;
            }
        }

        // Now choose u = first in A ∩ [0..Iend], v = first in B ∩ (Iend..end)
        int u=-1, v=-1;
        for (int i = 0; i <= Iend; ++i) {
            if (inA[i]) { u = R[i]; break; }
        }
        for (int i = Iend+1; i < Rsz; ++i) {
            if (inB[i]) { v = R[i]; break; }
        }
        if (u<0 || v<0) return false;

        // Let u' be successor of u in R; let v' be predecessor of v in R
        int idxu = find(R.begin(), R.end(), u) - R.begin();
        int idxv = find(R.begin(), R.end(), v) - R.begin();
        int up = R[idxu+1];
        int vp = R[idxv-1];

        if (!is_green(up, vp)) return false;

        // splice: [Start..u] + Q + reverse [up..v] + [vp..Finish]
        vector<int> newR;
        newR.insert(newR.end(), R.begin(), R.begin()+idxu+1);
        newR.insert(newR.end(), Q.begin(), Q.end());
        // reverse segment from up to v
        vector<int> revseg;
        for (int i = idxu+1; i <= idxv; ++i)
            revseg.push_back(R[i]);
        reverse(revseg.begin(), revseg.end());
        newR.insert(newR.end(), revseg.begin(), revseg.end());
        newR.insert(newR.end(), R.begin()+idxv, R.end());

        R = move(newR);
        merged = true;
    }

    return R.size() == (size_t)n;
}


bool HPA2(const vector<vector<int>>& graph, int n, int start, int finish) {
    // generate all the other colouring graphs
    // k, q_green, q_yellow, q_orange are adjustable
    int k = ceil(10 * log2(n));
    double q_green = 0.9;
    double q_yellow = 0.6;
    double q_orange = 0.3;
    vector<vector<int>> count(n, vector<int>(n, 0));
    for (int trial = 0; trial < k; ++trial) {
        auto red = get_random_coloring_graph(graph, n);
        for (int src = 0; src < n; src++) {
            for (auto& [dest, col] : red[src]) {
                if (col) {
                    count[src][dest]++;
                    count[dest][src]++;
                }
            }
        }
    }

    vector<vector<int>> green_graph(n);
    vector<vector<int>> yellow_graph(n);
    vector<vector<int>> orange_graph(n);
    for (int src = 0; src < n; src++) {
        for (int dest : graph[src]) {
            double freq = double(count[src][dest]) / k;
            if (freq >= q_green) {
                green_graph[src].push_back(dest);
                green_graph[dest].push_back(src);
            }
            else if (freq >= q_yellow) {
                yellow_graph[src].push_back(dest);
                yellow_graph[dest].push_back(src);
            }
            else if (freq >= q_orange) {
                orange_graph[src].push_back(dest);
                orange_graph[dest].push_back(src);
            }
        }
    }

    auto M = [&](int v) -> int {
        return min({orange_graph[v].size(), yellow_graph[v].size(), green_graph[v].size()});
    };

    int nlminus2 = ceil(static_cast<double>(n) / (std::pow(std::log2(static_cast<double>(n)), 2.0)));
    nlminus2 = max(nlminus2, 1);

    unordered_set<int> T_prime;
    for (int v = 0; v < n; v++) {
        if (M(v) < nlminus2) {
            T_prime.insert(v);
        }
    }
    vector<int> VV;
    for (int i = 0; i < n; i++) {
        VV.push_back(i);
    }

    if (T_prime.size() >= ceil(3.0 / -std::log2(q_green))) {
        return HPA3();
    }
    
    T_prime.insert(start);
    T_prime.insert(finish);

    vector<vector<int>> F;
    if (!stage2_build_envelope(start, finish, graph, T_prime, F)) return HPA3();

    if (!stage3_extend_orange(n, T_prime, F,orange_graph)) return HPA3();
    vector<int> Q0;
    if (!stage4_sew_yellow(F, Q0, yellow_graph)) return HPA3();
    vector<vector<int>> Qs;
    if (!stage5_partition(n,Q0,Qs,yellow_graph)) return HPA3();
    vector<bool> inTarr(n,false);
    for(int v:T_prime) inTarr[v]=true;
    vector<int> R = Q0;
    if (!stage6_merge(n, T_prime, R, Qs, green_graph, yellow_graph)) return HPA3();

    return true;
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
    if (int(even_outback.size()) >= LIM) return HPA2(graph, n, start, finish);

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
    if (int(odd_outback.size()) >= LIM) return HPA2(graph, n, start, finish);

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
            return HPA2(graph, n, start, finish);
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
    if (!even_outback.empty()) return HPA2(graph, n, start, finish);

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
    if (!odd_outback.empty()) return HPA2(graph, n, start, finish);

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
    cout << "Number of Vertices: " << n << ", number of edges: " << m << endl;


    
    auto start = chrono::steady_clock::now();

    bool found = false;
    for (int src = 0; src < n; src++) {
        for (int dest = src + 1; dest < n; dest++) {
            found = HPA1(graph, n, src, dest);
            if (found) break;
        }
        if (found) break;
    }

    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
    Util::display_time(duration, Util::micro);

    cout << (found ? "Yes" : "No") << endl;

    return 0;
}
