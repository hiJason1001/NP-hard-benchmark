#include <bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;

namespace Util {
    enum Second {
        milli,
        micro,
        nano,
    };

    std::string display_time(int64_t time, Second sec = milli) {
        std::ostringstream oss;

        if (sec == milli) {
            constexpr int64_t LIMIT_MS = 2LL * 3600 * 1000;
            if (time > LIMIT_MS) {
                return "Time Limit Exceeded\n";
            } else {
                int64_t hours = time / (3600 * 1000);
                int64_t minutes = (time % (3600 * 1000)) / (60 * 1000);
                int64_t seconds = (time % (60 * 1000)) / 1000;
                int64_t milliseconds = time % 1000;
                oss << hours << " hours "
                    << minutes << " minutes "
                    << seconds << " seconds "
                    << milliseconds << " milliseconds\n";
            }
        } else if (sec == micro) {
            constexpr int64_t LIMIT_US = 2LL * 3600 * 1000 * 1000;
            if (time > LIMIT_US) {
                return "Time Limit Exceeded\n";
            } else {
                int64_t hours = time / (3600LL * 1000 * 1000);
                int64_t minutes = (time % (3600LL * 1000 * 1000)) / (60LL * 1000 * 1000);
                int64_t seconds = (time % (60LL * 1000 * 1000)) / (1000LL * 1000);
                int64_t milliseconds = (time % (1000LL * 1000)) / 1000LL;
                int64_t microseconds = time % 1000LL;
                oss << hours << " hours "
                    << minutes << " minutes "
                    << seconds << " seconds "
                    << milliseconds << " milliseconds "
                    << microseconds << " microseconds\n";
            }
        } else if (sec == nano) {
            constexpr int64_t LIMIT_NS = 2LL * 3600 * 1000 * 1000 * 1000;
            if (time > LIMIT_NS) {
                return "Time Limit Exceeded\n";
            } else {
                int64_t hours = time / (3600LL * 1000 * 1000 * 1000);
                int64_t minutes = (time % (3600LL * 1000 * 1000 * 1000)) / (60LL * 1000 * 1000 * 1000);
                int64_t seconds = (time % (60LL * 1000 * 1000 * 1000)) / (1000LL * 1000 * 1000);
                int64_t milliseconds = (time % (1000LL * 1000 * 1000)) / (1000LL * 1000);
                int64_t microseconds = (time % (1000LL * 1000)) / 1000LL;
                int64_t nanoseconds = time % 1000LL;
                oss << hours << " hours "
                    << minutes << " minutes "
                    << seconds << " seconds "
                    << milliseconds << " milliseconds "
                    << microseconds << " microseconds "
                    << nanoseconds << " nanoseconds\n";
            }
        }

        return oss.str();
    }
};

bool extractGraphInfo(const std::string& tspFile, int& n, int& m) {
    std::ifstream fin(tspFile);
    if (!fin) return false;

    std::string line;
    bool inEdgeSection = false;
    m = 0;

    while (std::getline(fin, line)) {
        if (line.find("DIMENSION") != std::string::npos) {
            size_t pos = line.find(":");
            if (pos != std::string::npos)
                n = std::stoi(line.substr(pos + 1));
        } else if (line == "EDGE_DATA_SECTION") {
            inEdgeSection = true;
        } else if (line == "EOF") {
            break;
        } else if (inEdgeSection) {
            std::istringstream iss(line);
            int u, v;
            if (iss >> u >> v && u != -1 && v != -1)
                ++m;
        }
    }

    return true;
}

bool parseLKHOutput(const std::string& logFile, int numNodes) {
    std::ifstream fin(logFile);
    if (!fin) return false;

    std::string line;
    bool foundValidTour = false;

    std::regex timeRegex(R"((\d+ hours \d+ minutes \d+ seconds \d+ milliseconds \d+ microseconds))");
    std::regex costRegex(R"(Cost\s*=\s*(\d+))");

    while (std::getline(fin, line)) {
        std::smatch match;

        if (!foundValidTour && std::regex_search(line, match, costRegex)) {
            int cost = std::stoi(match[1]);
            if (cost == numNodes) {
                foundValidTour = true;
                break;
            }
        }
    }

    return foundValidTour;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <folder_with_par_tsp> <lkh_executable> <output_file>\n";
        return 1;
    }

    std::string folder = argv[1];
    std::string lkhExe = argv[2];
    std::string outputFile = argv[3];

    std::ofstream fout(outputFile);
    if (!fout) {
        std::cerr << "Failed to open output file\n";
        return 1;
    }

    std::vector<fs::directory_entry> files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".par") continue;
        files.push_back(entry);
    }

    std::sort(files.begin(), files.end(), [](const fs::directory_entry& a, const fs::directory_entry& b) {
        auto extractNumber = [](const std::string& name) {
            std::smatch match;
            if (std::regex_search(name, match, std::regex(R"((\d+)$)")))
                return std::stoi(match[1]);
            return 0;
        };
        return extractNumber(a.path().stem().string()) < extractNumber(b.path().stem().string());
    });

    for (const auto& entry : files) {
        std::string parPath = entry.path().string();
        std::string baseName = entry.path().stem().string();
        std::string tspPath = folder + "/" + baseName + ".tsp";
        std::string logPath = folder + "/" + baseName + ".log";

        std::string command = lkhExe + " " + parPath + " > " + logPath + " 2>&1";
        auto start = std::chrono::steady_clock::now();
        int ret = std::system(command.c_str());
        auto end = std::chrono::steady_clock::now();

        std::string timeStr;
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        long ms = duration / 1000;
        long us = duration % 1000;
        if (ret != 0) {
            ms = 0;
            us = 0;
        }

        timeStr = "0 hours 0 minutes 0 seconds " + std::to_string(ms) + " milliseconds " + std::to_string(us) + " microseconds";
        
        int numVertices = 0, numEdges = 0;
        if (!extractGraphInfo(tspPath, numVertices, numEdges)) {
            std::cerr << "Could not parse " << tspPath << "\n";
            continue;
        }
        
        bool found = parseLKHOutput(logPath, numVertices);
        if (ret != 0) {
            found = true;
        }
        fout << baseName << ".hcp\n";
        fout << "Number of Vertices: " << numVertices << ", number of edges: " << numEdges << "\n";
        fout << Util::display_time(duration, Util::micro);
        fout << (found ? "Yes" : "No") << "\n";
    }

    return 0;
}