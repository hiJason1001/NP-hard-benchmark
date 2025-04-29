#pragma once

#include <iostream>
#include <chrono>
#include <vector>

namespace Util {
    bool get_adjMatrix(const std::string& filename, int& n, int& m, std::vector<std::vector<int>>& adjMatrix) {
        std::ifstream file(filename);
    
        if (!file.is_open()) {
            std::cerr << "Failed to open the file!" << std::endl;
            return false;
        }

        file >> n >> m;
        adjMatrix.assign(n, std::vector<int>(n, 0));
        
        for (int i = 0; i < m; i++) {
            int u, v;
            file >> u >> v;
            adjMatrix[u][v] = 1;
            adjMatrix[v][u] = 1;
        }

        return true;
    }

    bool get_adjList(const std::string& filename, int& n, int& m, std::vector<std::vector<int>>& adjList) {
        std::ifstream file(filename);
    
        if (!file.is_open()) {
            std::cerr << "Failed to open the file!" << std::endl;
            return false;
        }

        file >> n >> m;
        adjList.assign(n, std::vector<int>());
        
        for (int i = 0; i < m; i++) {
            int u, v;
            file >> u >> v;
            adjList[u].push_back(v);
            adjList[v].push_back(u);
        }

        return true;
    }

    enum Second {
        milli,
        micro,
        nano,
    };

    void display_time(auto time, Second sec=milli) {
        if (sec == milli) {
            constexpr long long LIMIT_MS = 2LL * 3600 * 1000;
            if (time > LIMIT_MS) {
                std::cout << "Time Limit Exceeded " << std::endl;
            } else {
                int64_t hours = time / (3600 * 1000);
                int64_t minutes = (time % (3600 * 1000)) / (60 * 1000);
                int64_t seconds = (time % (60 * 1000)) / 1000;
                int64_t milliseconds = time % 1000;
                std::cout << hours << " hours " 
                          << minutes << " minutes " 
                          << seconds << " seconds " 
                          << milliseconds << " milliseconds " << std::endl;
            }
        }
        else if (sec == micro) {
            constexpr long long LIMIT_US = 2LL * 3600 * 1000 * 1000;
            if (time > LIMIT_US) {
                std::cout << "Time Limit Exceeded " << std::endl;
            } else {
                int64_t hours = time / (3600LL * 1000 * 1000);
                int64_t minutes = (time % (3600LL * 1000 * 1000)) / (60LL * 1000 * 1000);
                int64_t seconds = (time % (60LL * 1000 * 1000)) / (1000LL * 1000);
                int64_t milliseconds = (time % (1000LL * 1000)) / 1000LL;
                int64_t microseconds = time % 1000LL;

                std::cout << hours << " hours "
                        << minutes << " minutes "
                        << seconds << " seconds "
                        << milliseconds << " milliseconds "
                        << microseconds << " microseconds " << std::endl;
            }
        }
        else if (sec == nano) {
            constexpr long long LIMIT_NS = 2LL * 3600 * 1000 * 1000 * 1000;
            if (time > LIMIT_NS) {
                std::cout << "Time Limit Exceeded " << std::endl;
            } else {
                int64_t hours = time / (3600LL * 1000 * 1000 * 1000);
                int64_t minutes = (time % (3600LL * 1000 * 1000 * 1000)) / (60LL * 1000 * 1000 * 1000);
                int64_t seconds = (time % (60LL * 1000 * 1000 * 1000)) / (1000LL * 1000 * 1000);
                int64_t milliseconds = (time % (1000LL * 1000 * 1000)) / (1000LL * 1000);
                int64_t microseconds = (time % (1000LL * 1000)) / 1000LL;
                int64_t nanoseconds = time % 1000LL;
        
                std::cout << hours << " hours "
                          << minutes << " minutes "
                          << seconds << " seconds "
                          << milliseconds << " milliseconds "
                          << microseconds << " microseconds "
                          << nanoseconds << " nanoseconds " << std::endl;
            }
        }
    }
};