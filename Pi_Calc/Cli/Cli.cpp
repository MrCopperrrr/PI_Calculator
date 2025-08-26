#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream> // Required for std::stringstream
#include "..\Core\Core.h"

void write_to_file(const std::string& filename, const std::string& content)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    outFile << content;
    outFile.close();
    std::cout << "Successfully wrote Pi to " << filename << std::endl;
}

/**
 * @brief Formats a duration in milliseconds into a human-readable string.
 *
 * @param ms The duration to format, as a std::chrono::milliseconds object.
 * @return A string in the format: X year(s) X month(s) X day(s)...
 */
std::string format_duration(std::chrono::milliseconds ms)
{
    if (ms.count() == 0) {
        return "0 millisecond(s)";
    }

    // Use long long to avoid overflow for very long durations
    long long total_ms = ms.count();

    // Define durations in milliseconds (using average values for month/year)
    constexpr long long ms_in_second = 1000;
    constexpr long long ms_in_minute = ms_in_second * 60;
    constexpr long long ms_in_hour = ms_in_minute * 60;
    constexpr long long ms_in_day = ms_in_hour * 24;
    constexpr long long ms_in_month = static_cast<long long>(ms_in_day * 30.4375); // Average month
    constexpr long long ms_in_year = static_cast<long long>(ms_in_day * 365.25);  // Average year

    // Calculate each unit
    long long years = total_ms / ms_in_year;
    total_ms %= ms_in_year;

    long long months = total_ms / ms_in_month;
    total_ms %= ms_in_month;

    long long days = total_ms / ms_in_day;
    total_ms %= ms_in_day;

    long long hours = total_ms / ms_in_hour;
    total_ms %= ms_in_hour;

    long long minutes = total_ms / ms_in_minute;
    total_ms %= ms_in_minute;

    long long seconds = total_ms / ms_in_second;
    total_ms %= ms_in_second;

    long long milliseconds = total_ms;

    // Build the output string
    std::stringstream ss;
    if (years > 0)   ss << years << " year(s) ";
    if (months > 0)  ss << months << " month(s) ";
    if (days > 0)    ss << days << " day(s) ";
    if (hours > 0)   ss << hours << " hour(s) ";
    if (minutes > 0) ss << minutes << " minute(s) ";
    if (seconds > 0) ss << seconds << " second(s) ";

    // Always show milliseconds, especially if the duration is less than a second
    if (milliseconds > 0 || ss.str().empty()) {
        ss << milliseconds << " millisecond(s)";
    }

    std::string result = ss.str();
    // Remove potential trailing space
    if (!result.empty() && result.back() == ' ') {
        result.pop_back();
    }

    return result;
}


int main(int argc, char* argv[])
{
    if (argc < 3 || argc > 4)
    {
        std::cerr << "Usage: " << argv[0] << " <digits_after_decimal> <num_threads> [output_file]" << std::endl;
        return 1;
    }

    try
    {
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

        // Calculate the duration and cast it to milliseconds
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        // Use the new formatting function
        std::cout << "Calculation finished in " << format_duration(elapsed_ms) << "." << std::endl;

        if (argc == 4)
        {
            std::string filename = argv[3];
            write_to_file(filename, pi);
        }
        else
        {
            int preview_digits = (digits > 100) ? 100 : digits;
            std::cout << "Pi: " << pi.substr(0, preview_digits + 2) << "..." << std::endl;
        }
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "Error: Invalid number provided for digits or threads." << std::endl;
        return 1;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Error: Number is too large. Please enter a reasonable number of digits." << std::endl;
        return 1;
    }

    return 0;
}