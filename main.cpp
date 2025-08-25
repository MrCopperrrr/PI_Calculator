#include "compute_pi.h"
using namespace std;
using namespace std::chrono;


int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <threads> <digits>\n";
        return 1;
    }

    int threads = 8;
    long long digits = 100000;
    string outfile = "pi_bs_par_stream_output.txt";

    if (argc >= 2) threads = atoi(argv[1]);
    if (argc >= 3) digits = atoll(argv[2]);
    if (argc >= 4) outfile = string(argv[3]);

    auto t0 = high_resolution_clock::now();
    compute_pi_bs_parallel(digits, threads, outfile);
    auto t1 = high_resolution_clock::now();

    long long ms = duration_cast<milliseconds>(t1 - t0).count();
    long long hours = ms / (1000LL*60*60);
    ms %= (1000LL*60*60);
    long long minutes = ms / (1000LL*60);
    ms %= (1000LL*60);
    long long seconds = ms / 1000LL;
    long long milliseconds = ms % 1000LL;

    cout << "Execution time: "
         << hours << " hour "
         << minutes << " minute "
         << seconds << " second "
         << milliseconds << " millisecond" << endl;

    return 0;
}
