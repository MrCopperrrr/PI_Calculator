#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
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
        std::chrono::duration<double> elapsed = end_time - start_time;

        std::cout << "Calculation finished in " << elapsed.count() << " seconds." << std::endl;

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