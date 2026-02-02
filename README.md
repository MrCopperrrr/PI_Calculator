# High-Performance PI Calculator

Dự án tính toán số Pi với hiệu năng cao sử dụng thuật toán **Chudnovsky**, kỹ thuật **Binary Splitting** và trình nhân số lớn dựa trên **NTT (Number Theoretic Transform)**. Mục tiêu là đạt tốc độ xử lý tiệm cận các công cụ hàng đầu hiện nay như `y-cruncher`.

## 🚀 Tính năng hiện tại
- **Chudnovsky Algorithm**: Fast convergence (14 digits per iteration).
- **GMP Integration**: Uses GNU Multi-Precision library for assembly-optimized math.
- **Binary Splitting**: Optimized sequence calculation, parallelized with OpenMP.
- **High-Resolution Timer**: Format `xx h xx m xx s xx ms`.

## 🚀 Performance (Benchmarks)
- **1M digits**: ~0.17s
- **10M digits**: ~2.7s
- **100M digits**: ~41s (Total including I/O), ~24s (Calculation only)

## 🛠 Yêu cầu hệ thống
- **Compiler**: GCC/G++ (C++17+)
- **Build System**: CMake (3.10+)
- **Library**: **GMP** (libgmp-dev)
- **Parallelism**: OpenMP

## 🔨 Hướng dẫn cài đặt và Chạy

### 1. Chuẩn bị thư mục Build
Mở terminal tại thư mục gốc của dự án và chạy:
```bash
mkdir build
cd build
```

### 2. Cấu hình với CMake
Sử dụng các cờ tối ưu hóa cao nhất cho phần cứng của bạn:
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### 3. Biên dịch
```bash
cmake --build . --config Release
```

### 4. Chạy chương trình
Sau khi biên dịch thành công, chạy file thực thi:
- **Windows**: `.\Release\pi_calc.exe` hoặc `.\pi_calc.exe` (tùy bộ sinh mã).
- **Linux/macOS**: `./pi_calc`

### 5. Build and run
```bash
cmake --build . --config Release ; .\pi_calc.exe 100k
```

## 📈 Lộ trình phát triển (Roadmap)
- [x] Nền tảng BigInt & Binary Splitting.
- [x] Hệ thống đo thời gian an toàn.
- [x] Bộ nhân NTT cơ bản.
- [ ] Phép chia số lớn (Long Division) & Căn bậc hai (Newton-Raphson).
- [ ] Đa luồng (Multi-threading) với OpenMP.
- [ ] Tối ưu hóa bộ nhớ Swap (Disk I/O) cho hàng tỷ chữ số.

## ⚠️ Lưu ý về hiệu năng
Để đạt tốc độ cao nhất, chương trình được cấu hình mặc định để tận dụng tối đa kiến trúc CPU hiện tại qua cờ `-mavx2 -march=native`. Nếu bạn gặp lỗi khi chạy trên các máy quá cũ, hãy điều chỉnh trong `CMakeLists.txt`.

---
*Dự án đang trong quá trình phát triển để chinh phục các mốc giới hạn mới của số Pi.*
