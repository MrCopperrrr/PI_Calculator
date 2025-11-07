#include "src/compute_pi.h"
using namespace std;
using namespace std::chrono;

void preview_pi_from_file(const std::string& filename, int num_digits_after_decimal) {
    std::cout << "\n--- Previewing Pi from output file ---\n";

    // Mở file để đọc
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " to preview Pi.\n";
        return;
    }

    // Tổng số ký tự cần đọc: "3." (2 ký tự) + số chữ số thập phân
    int chars_to_read = 2 + num_digits_after_decimal;

    // Tạo một buffer để chứa các ký tự đọc được
    std::vector<char> buffer(chars_to_read);

    // Đọc một khối ký tự từ file vào buffer
    file.read(buffer.data(), chars_to_read);

    // file.gcount() trả về số ký tự đã thực sự đọc được,
    // phòng trường hợp file ngắn hơn yêu cầu.
    std::streamsize chars_actually_read = file.gcount();

    // In ra các ký tự đã đọc được
    std::cout << "First " << num_digits_after_decimal << " decimal digits of Pi: ";
    std::cout.write(buffer.data(), chars_actually_read);
    std::cout << "...\n" << std::endl;

    // file sẽ tự động được đóng khi biến 'file' ra khỏi scope
}

void view_last_pi_digits(const std::string& filename, int num_digits = 20) {

    // Mở file để đọc
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " to read Pi digits.\n";
        return;
    }

    // Di chuyển con trỏ đến cuối file để xác định kích thước
    file.seekg(0, std::ios::end);
    std::streampos file_size = file.tellg();
    
    if (file_size <= 2) { // File quá ngắn (chỉ có "3." hoặc ít hơn)
        std::cerr << "Error: File too short to contain Pi digits.\n";
        return;
    }

    // Tính toán vị trí bắt đầu đọc (file_size - num_digits)
    // Vì file bắt đầu bằng "3.", ta bỏ qua 2 ký tự đầu
    std::streampos start_pos = std::max((std::streampos)2, file_size - (std::streamoff)num_digits);
    
    // Di chuyển con trỏ đến vị trí bắt đầu đọc
    file.seekg(start_pos);

    // Đọc các ký tự cuối
    std::string last_digits;
    std::getline(file, last_digits);
    
    // Lấy đúng num_digits ký tự cuối (bỏ qua ký tự xuống dòng nếu có)
    if (last_digits.length() > num_digits) {
        last_digits = last_digits.substr(last_digits.length() - num_digits);
    }

    std::cout << "Last " << num_digits << " decimal digits of Pi: ..." << last_digits << std::endl;
    std::cout << std::endl;
}



long long parse_digits(const string &digits_str) {
    char suffix = toupper(digits_str.back());
    long long multiplier = 1;

    if (suffix == 'K') multiplier = 1000;
    else if (suffix == 'M') multiplier = 1000000;
    else if (suffix == 'B') multiplier = 1000000000;

    if (isdigit(suffix)) {
        return atoll(digits_str.c_str());
    } else {
        return atoll(digits_str.substr(0, digits_str.size() - 1).c_str()) * multiplier;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <threads> <digits>\n";
        return 1;
    }

    int threads = 16;
    long long digits = 100000;
    string outfile = "output.txt";

    if (argc >= 2) threads = atoi(argv[1]);
    if (argc >= 3) digits = parse_digits(argv[2]);
    if (argc >= 4) outfile = string(argv[3]);

    auto t0 = high_resolution_clock::now();
    compute_pi_bs_parallel(digits, threads, outfile);
    auto t1 = high_resolution_clock::now();

    preview_pi_from_file(outfile, 20);
    view_last_pi_digits(outfile, 20);

    long long ms = duration_cast<milliseconds>(t1 - t0).count();
    long long hours = ms / (1000LL*60*60);
    ms %= (1000LL*60*60);
    long long minutes = ms / (1000LL*60);
    ms %= (1000LL*60);
    long long seconds = ms / 1000LL;
    long long milliseconds = ms % 1000LL;

    cout << "Calculation time: "
         << hours << " hour "
         << minutes << " minute "
         << seconds << " second "
         << milliseconds << " millisecond" << endl;


    return 0;
}
