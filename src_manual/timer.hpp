#pragma once
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>


namespace pi {

class Timer {
public:
  Timer() : start_time(std::chrono::high_resolution_clock::now()) {}

  void reset() { start_time = std::chrono::high_resolution_clock::now(); }

  std::string elapsed_str() const {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time)
            .count();

    // Ensure non-negative duration
    if (duration < 0)
      duration = 0;

    long long ms = duration % 1000;
    long long s = (duration / 1000) % 60;
    long long m = (duration / (1000 * 60)) % 60;
    long long h = duration / (1000 * 60 * 60);

    std::stringstream ss;
    ss << h << " h " << m << " m " << s << " s " << ms << " ms";
    return ss.str();
  }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

} // namespace pi
