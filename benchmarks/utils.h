#include <string>
#include <vector>

std::vector<int> split_into_ints(const std::string &s, const char sep)
{
    std::vector<int> ints;
    size_t pos = s.find(sep, 0), prev = 0;
    while (pos != std::string::npos) {
        ints.push_back(atoi(s.substr(prev, pos - prev).c_str()));
        prev = pos + 1;
        pos = s.find(sep, pos + 1);
    }
    ints.push_back(atoi(s.substr(prev).c_str()));
    return ints;
}
