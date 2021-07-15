#include <fiuncho/utils/Result.h>
#include <functional>
#include <unordered_map>
#include <vector>

bool has_repeated_elements(const std::vector<Result<int, float>> &results)
{
    // String function for a vector of integers
    auto vec_to_str = [](const std::vector<int> &v) {
        std::string s;
        for (auto i = v.begin(); i < v.end(); ++i) {
            s += std::to_string(*i);
        }
        return s;
    };
    std::unordered_map<std::string, int> count;
    // Count occurrences
    for (auto r = results.begin(); r != results.end(); ++r) {
        auto s = vec_to_str(r->combination);
        ++count[s];
        if (count[s] > 1) {
            return true;
        }
    }
    // Otherwise return false
    return false;
}

bool ascending_combinations(const std::vector<Result<int, float>> &v)
{
    // For each result
    for (auto r = v.begin(); r != v.end(); ++r) {
        // Check if every index is strictly superior than its previous
        for (size_t i = 1; i < r->combination.size(); ++i) {
            if (r->combination[i - 1] >= r->combination[i]) {
                // If that's the case, return false
                return false;
            }
        }
    }
    // Otherwise, return true
    return true;
}

bool matches_mpi3snp_output(const std::vector<Result<int, float>> &v)
{
    std::vector<std::vector<int>> mpi3snp_output{
        {0, 8, 9}, {0, 5, 9}, {5, 6, 9}, {0, 7, 9}, {1, 3, 9}, {0, 1, 9},
        {0, 3, 9}, {1, 5, 7}, {0, 6, 9}, {0, 3, 4}, {3, 8, 9}, {0, 4, 5},
        {0, 5, 6}, {0, 4, 9}, {5, 7, 9}, {4, 5, 7}, {5, 7, 8}, {0, 5, 7},
        {3, 5, 6}, {3, 5, 7}, {3, 6, 9}, {0, 5, 8}, {0, 1, 4}, {4, 5, 9},
        {1, 6, 9}, {0, 4, 7}, {0, 2, 5}, {0, 2, 9}, {0, 3, 8}, {4, 5, 6},
        {5, 6, 7}, {1, 8, 9}, {2, 6, 7}, {0, 3, 7}, {0, 2, 3}, {0, 1, 8},
        {0, 4, 8}, {0, 2, 8}, {1, 7, 9}, {6, 7, 8}, {3, 7, 8}, {0, 2, 6},
        {2, 3, 7}, {0, 1, 3}, {1, 3, 7}, {2, 3, 9}, {0, 2, 4}, {1, 3, 8},
        {3, 7, 9}, {2, 4, 6}, {3, 5, 9}, {0, 3, 5}, {2, 3, 5}, {4, 5, 8},
        {1, 7, 8}, {1, 5, 9}, {0, 7, 8}, {4, 8, 9}, {2, 8, 9}, {3, 5, 8},
        {0, 1, 5}, {0, 3, 6}, {7, 8, 9}, {5, 8, 9}, {2, 4, 5}, {1, 3, 5},
        {3, 4, 9}, {2, 5, 7}, {2, 3, 6}, {0, 4, 6}, {2, 6, 8}, {3, 4, 5},
        {2, 7, 8}, {6, 8, 9}, {4, 6, 9}, {2, 3, 8}, {0, 2, 7}, {2, 6, 9},
        {0, 6, 8}, {6, 7, 9}, {2, 3, 4}, {0, 1, 2}, {1, 4, 9}, {0, 1, 6},
        {1, 5, 8}, {3, 4, 8}, {1, 3, 4}, {0, 1, 7}, {1, 6, 8}, {1, 2, 6},
        {2, 5, 6}, {1, 6, 7}, {3, 4, 7}, {4, 7, 9}, {3, 4, 6}, {1, 2, 8},
        {0, 6, 7}, {2, 5, 9}, {2, 7, 9}, {1, 4, 7}
    };

    // For each result
    auto r1 = v.begin();
    auto r2 = mpi3snp_output.begin();
    for (; r1 != v.end() && r2 != mpi3snp_output.end(); ++r1, ++r2) {
        // Check if its the same as in mpi3snp
        if (r1->combination != *r2) {
            // If that's the case, return false
            return false;
        }
    }
    // Otherwise, return true
    return true;
}