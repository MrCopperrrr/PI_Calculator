#include <iostream>
#include <string>
#include <chrono>
#include "../Core/core_api.h" // Include the API from the Core project

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <threads> <digits> [outfile]\n";
        std::cerr << "Example: " << argv[0] << " 8 100000 pi_output.txt\n";
        return 1;
    }

    int threads = 8;
    long long digits = 100000;
    std::string outfile = "pi_output.txt";

    try {
        if (argc >= 2) threads = std::stoi(argv[1]);
        if (argc >= 3) digits = std::stoll(argv[2]);
        if (argc >= 4) outfile = std::string(argv[3]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Starting Pi calculation..." << std::endl;
    std::cout << "Threads: " << threads << std::endl;
    std::cout << "Digits: " << digits << std::endl;
    std::cout << "Output File: " << outfile << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Call the function from Core.dll
    compute_pi(digits, threads, outfile.c_str());

    auto t1 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

    long long ms = ms_int.count();
    long long hours = ms / (1000LL * 60 * 60);
    ms %= (1000LL * 60 * 60);
    long long minutes = ms / (1000LL * 60);
    ms %= (1000LL * 60);
    long long seconds = ms / 1000LL;
    long long milliseconds = ms % 1000LL;

    std::cout << "\nCalculation finished!" << std::endl;
    std::cout << "Execution time: "
        << hours << "h "
        << minutes << "m "
        << seconds << "s "
        << milliseconds << "ms" << std::endl;

    return 0;
}