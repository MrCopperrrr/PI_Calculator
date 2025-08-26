#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>
#include "..\Core\Core.h"
#include "..\GPU\gpu_multiplier.h" // Include the GPU header for the initialization function

void write_to_file(const std::string& filename, const std::string& content) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    outFile << content;
    outFile.close();
    std::cout << "Successfully wrote Pi to " << filename << std::endl;
}

std::string format_duration(std::chrono::milliseconds ms) {
    if (ms.count() == 0) return "0 millisecond(s)";
    long long total_ms = ms.count();
    constexpr long long ms_in_second = 1000;
    constexpr long long ms_in_minute = ms_in_second * 60;
    constexpr long long ms_in_hour = ms_in_minute * 60;
    constexpr long long ms_in_day = ms_in_hour * 24;
    long long days = total_ms / ms_in_day; total_ms %= ms_in_day;
    long long hours = total_ms / ms_in_hour; total_ms %= ms_in_hour;
    long long minutes = total_ms / ms_in_minute; total_ms %= ms_in_minute;
    long long seconds = total_ms / ms_in_second; total_ms %= ms_in_second;
    long long milliseconds = total_ms;
    std::stringstream ss;
    if (days > 0) ss << days << " day(s) ";
    if (hours > 0) ss << hours << " hour(s) ";
    if (minutes > 0) ss << minutes << " minute(s) ";
    if (seconds > 0) ss << seconds << " second(s) ";
    if (milliseconds > 0 || ss.str().empty()) ss << milliseconds << " millisecond(s)";
    std::string result = ss.str();
    if (!result.empty() && result.back() == ' ') result.pop_back();
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <digits_after_decimal> <num_threads> [output_file]" << std::endl;
        return 1;
    }

    // --- GPU INITIALIZATION CHECK ---
    std::cout << "Initializing GPU context..." << std::endl;
    if (!initialize_gpu_context()) {
        std::cerr << "Aborting due to GPU initialization failure. Please check your NVIDIA drivers and hardware." << std::endl;
        // Pause to allow user to see the error message before the console closes
        std::cout << "Press Enter to exit...";
        std::cin.get();
        return 1;
    }
    std::cout << "GPU context initialized successfully." << std::endl;
    // --- END OF CHECK ---

    try {
        int digits = std::stoi(argv[1]);
        int num_threads = std::stoi(argv[2]);
        if (digits <= 0 || num_threads <= 0) {
            std::cerr << "Error: Number of digits and threads must be positive." << std::endl;
            return 1;
        }

        std::cout << "Calculating " << digits << " digits of Pi after the decimal point using " << num_threads << " threads..." << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now();
        std::string pi = PiCalculator::Chudnovsky::calculate(digits, num_threads);
        auto end_time = std::chrono::high_resolution_clock::now();

        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Calculation finished in " << format_duration(elapsed_ms) << "." << std::endl;

        if (argc == 4) {
            std::string filename = argv[3];
            write_to_file(filename, pi);
        }
        else {
            int preview_digits = (digits > 100) ? 100 : digits;
            std::cout << "Pi: " << pi.substr(0, preview_digits + 2) << "..." << std::endl;
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid number provided for digits or threads." << std::endl;
        return 1;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Error: Number is too large. Please enter a reasonable number of digits." << std::endl;
        return 1;
    }

    return 0;
}