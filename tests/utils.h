#include <fiuncho/utils/Result.h>
#include <cmath>
#include <unordered_map>
#include <vector>

template <typename T>
std::string vec_to_str(const std::vector<T> &v) {
    std::string s;
    for (auto i = v.begin(); i < v.end(); ++i) {
        s += std::to_string(*i);
        s += ",";
    }
    return s;
};

bool has_repeated_elements(const std::vector<Result<int, float>> &results)
{
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
    std::vector<std::pair<std::vector<int>, float>> mpi3snp_output{
        {{4, 8, 9}, 0.198944}, {{6, 8, 9}, 0.198866}, {{3, 8, 9}, 0.198622},
        {{7, 8, 9}, 0.198564}, {{5, 8, 9}, 0.198526}, {{0, 8, 9}, 0.198519},
        {{1, 8, 9}, 0.198300}, {{2, 8, 9}, 0.198216}, {{5, 6, 7}, 0.001916},
        {{0, 5, 7}, 0.001585}, {{0, 4, 5}, 0.001559}, {{1, 6, 7}, 0.001548},
        {{1, 2, 7}, 0.001522}, {{3, 5, 8}, 0.001495}, {{4, 6, 7}, 0.001482},
        {{6, 7, 8}, 0.001478}, {{1, 4, 5}, 0.001474}, {{0, 1, 9}, 0.001451},
        {{0, 4, 9}, 0.001429}, {{3, 5, 7}, 0.001391}, {{4, 7, 9}, 0.001374},
        {{0, 2, 7}, 0.001346}, {{4, 5, 9}, 0.001341}, {{0, 1, 2}, 0.001338},
        {{2, 3, 5}, 0.001320}, {{1, 7, 9}, 0.001315}, {{1, 5, 7}, 0.001315},
        {{2, 6, 8}, 0.001313}, {{1, 7, 8}, 0.001309}, {{0, 4, 8}, 0.001266},
        {{0, 3, 5}, 0.001260}, {{1, 4, 9}, 0.001258}, {{3, 7, 9}, 0.001254},
        {{0, 1, 7}, 0.001249}, {{0, 3, 7}, 0.001245}, {{2, 6, 7}, 0.001240},
        {{0, 2, 4}, 0.001238}, {{4, 5, 7}, 0.001235}, {{5, 7, 8}, 0.001220},
        {{0, 1, 5}, 0.001203}, {{4, 5, 8}, 0.001201}, {{0, 1, 3}, 0.001183},
        {{2, 7, 8}, 0.001179}, {{2, 4, 6}, 0.001165}, {{2, 5, 7}, 0.001165},
        {{0, 1, 4}, 0.001163}, {{3, 6, 7}, 0.001153}, {{0, 2, 3}, 0.001152},
        {{5, 6, 8}, 0.001148}, {{4, 5, 6}, 0.001136}, {{3, 7, 8}, 0.001119},
        {{0, 4, 7}, 0.001115}, {{6, 7, 9}, 0.001113}, {{3, 4, 5}, 0.001111},
        {{2, 3, 6}, 0.001110}, {{1, 5, 8}, 0.001107}, {{2, 3, 7}, 0.001107},
        {{0, 6, 7}, 0.001070}, {{0, 5, 8}, 0.001070}, {{0, 1, 6}, 0.001068},
        {{0, 6, 8}, 0.001067}, {{0, 2, 5}, 0.001058}, {{2, 3, 8}, 0.001057},
        {{5, 7, 9}, 0.001055}, {{2, 4, 7}, 0.001054}, {{0, 2, 8}, 0.001054},
        {{1, 3, 7}, 0.001050}, {{2, 5, 8}, 0.001047}, {{3, 4, 6}, 0.001047},
        {{0, 3, 8}, 0.001038}, {{0, 3, 6}, 0.001034}, {{2, 4, 9}, 0.001031},
        {{1, 3, 9}, 0.001030}, {{3, 4, 9}, 0.001029}, {{0, 3, 9}, 0.001028},
        {{1, 2, 9}, 0.001026}, {{2, 7, 9}, 0.001016}, {{1, 2, 8}, 0.001008},
        {{1, 4, 7}, 0.001007}, {{1, 3, 8}, 0.000986}, {{4, 6, 9}, 0.000982},
        {{1, 5, 9}, 0.000981}, {{1, 3, 5}, 0.000981}, {{0, 3, 4}, 0.000963},
        {{0, 4, 6}, 0.000960}, {{3, 6, 8}, 0.000932}, {{2, 5, 6}, 0.000931},
        {{0, 2, 6}, 0.000926}, {{0, 7, 9}, 0.000922}, {{0, 1, 8}, 0.000917},
        {{0, 5, 6}, 0.000917}, {{2, 4, 5}, 0.000910}, {{1, 2, 5}, 0.000899},
        {{3, 5, 6}, 0.000896}, {{1, 5, 6}, 0.000883}, {{0, 7, 8}, 0.000881},
        {{4, 6, 8}, 0.000872}, {{2, 6, 9}, 0.000854}, {{1, 6, 8}, 0.000806},
        {{0, 2, 9}, 0.000801},
    };
    // Create a map in which the MI of each combination from MPI3SNP is
    // normalized
    std::unordered_map<std::string, float> norm_mi;
    const float mpi3snp_max_mi = mpi3snp_output[0].second;
    for (size_t i = 0; i < mpi3snp_output.size(); ++i){
        auto p = mpi3snp_output[i];
        auto s = vec_to_str(p.first);
        norm_mi[s] = p.second / mpi3snp_max_mi;
    }
    // Compare the normalized MI from MPI3SNP with the same value from fiuncho
    const float fiuncho_max_mi = v[0].val;
    for (auto r : v){
        auto s = vec_to_str(r.combination);
        if (std::abs(r.val / fiuncho_max_mi - norm_mi[s]) > 0.0001){
            return false;
        }
    }
    return true;
}
