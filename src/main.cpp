#include <cpuid.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

size_t max_size = 1 << 25;
int alignment = 1 << 12;  // 4k
char* memory = (char*)aligned_alloc(alignment, max_size);
void reset_cache() {
    static volatile char* buffer = (char*)aligned_alloc(alignment, max_size);
    for (int i = 0; i < max_size; i += 16)
        buffer[i] = buffer[i - 1];
}

void make_accesses(char* cur_addr) {
    while (cur_addr != nullptr) {
        cur_addr = *(char**)cur_addr;
    }
}

void create_accesses_array(char* array, std::vector<size_t>& pos) {
    if (pos.size() > 1) {
        auto rnd = std::default_random_engine{};
        std::shuffle(pos.begin() + 1, pos.end(), rnd);
    }
    for (int i = 0; i < pos.size() - 1; i++)
        *reinterpret_cast<char**>(array + pos[i]) = array + pos[i + 1];
    *reinterpret_cast<char**>(array + pos[pos.size() - 1]) = nullptr;
}

std::vector<size_t> create_pos(size_t start, size_t count, size_t stride) {
    std::vector<size_t> pos(count);
    for (int i = 0; i < count; i++)
        pos[i] = start + i * stride;

    return pos;
}

int measure(char* array, int accesses_cnt) {
        make_accesses(array);
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < accesses_cnt; i++) {
        make_accesses(array);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds ms_double = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return ms_double.count();
}

void print_cpu_version() {
    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0, 0, 0, 0};

    __cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    unsigned int nExIds = CPUInfo[0];

    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    for (unsigned int i = 0x80000000; i <= nExIds; ++i) {
        __cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);

        if (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
    }

    std::cout << "CPU Version: " << CPUBrandString << std::endl;
}

int measure_size(int cur_size) {
    int sum = 0;
    for (int bench = 0; bench < 5; bench++) {
        auto pos = create_pos(0, cur_size / 16, 16);
        create_accesses_array(memory, pos);
        sum += measure(memory, 10000);
    }
    return sum;
}

int measure_assoc(int cache_size, int cur_assoc) {
    int sum = 0;
    for (int bench = 0; bench < 25; bench++) {
        auto pos = create_pos(0, cur_assoc, cache_size);
        create_accesses_array(memory, pos);
        sum += measure(memory, 1000000);
    }
    return sum;
}

int measure_block(int cache_size, int assoc, int cur_size) {
    int sum = 0;
    for (int bench = 0; bench < 100; bench++) {
        int way_size = cache_size / assoc;
        auto pos1 = create_pos(0, assoc / 2, way_size);
        auto pos2 = create_pos(way_size * (assoc / 2) + cur_size, assoc / 2 + 1, way_size);
        pos1.insert(pos1.end(), pos2.begin(), pos2.end());
        create_accesses_array(memory, pos1);
        sum += measure(memory, 1000000);
    }
    return sum;
}

// double measure_ws_and_assoc() {

//     auto pos = create_pos(0, cur_assoc, cache_size);
//     create_accesses_array(memory, pos);
//     return measure(memory);
// }

class JumpDetector {
   public:
    JumpDetector(double jump_coef) : jump_coef(jump_coef) {}
    void add_measure(int measure) {
        if (prev_measure != 0) {
            int cur_delta = measure - prev_measure;
            if (delta_sum != 0) {
                int mean_delta = delta_sum / measures_count;
                int diff = std::abs(cur_delta - mean_delta);
                std::cout << " jump?: " << (double)diff / std::abs(mean_delta) << " mean: " << mean_delta << " cur: " << cur_delta << "\n";
                if (diff > std::abs(mean_delta) * jump_coef)
                    is_jump_ = true;
            }
            delta_sum += cur_delta;
            measures_count++;
        }
        prev_measure = measure;
    }

    bool is_jump() const { return is_jump_; }

   private:
    bool is_jump_ = false;
    double jump_coef;
    int measures_count = 0;
    int delta_sum = 0;
    int prev_measure = 0;
};

int assoc() {
    double prev_time = 0;
    std::ofstream fout("measure_assoc");
    int res_assoc = 0;
    JumpDetector jump_detector(5);
    for (int cur_assoc = 1; cur_assoc <= 16; cur_assoc++) {
        auto cur_time = measure_assoc(32 * 1024, cur_assoc);
        std::cout << "cur_assoc: " << cur_assoc << ", time: " << (cur_time) << ", delta: " << (cur_time - prev_time) << "\n";
        fout << cur_assoc << ' ' << cur_time << '\n';
        jump_detector.add_measure(cur_time);
        if (jump_detector.is_jump()) {
            res_assoc = cur_assoc - 1;
            break;
        }
        prev_time = cur_time;
    }
    fout.close();
    return res_assoc;
}

int size() {
    int res_size = 0;
    double prev_time = 0;
    std::ofstream fout("measure_size");
    JumpDetector jump_detector(2);
    for (int cur_size = 1; cur_size <= 64; cur_size++) {
        auto cur_time = measure_size(cur_size * 1024);
        //std::cout << "cur_size: " << cur_size << ", time: " << (cur_time) << ", delta: " << (cur_time - prev_time) << "\n";
        fout << cur_size << ' ' << cur_time << '\n';

        jump_detector.add_measure(cur_time);
        if (jump_detector.is_jump()) {
            res_size = cur_size - 1;
            break;
        }
        prev_time = cur_time;
    }
    fout.close();
    return res_size;
}

int block() {
    int res_block = 0;
    double prev_time = 0;
    std::ofstream fout("measure_block");
    JumpDetector jump_detector(10);
    for (int cur_block = 16; cur_block <= 128; cur_block += 8) {
        auto cur_time = measure_block(32 * 1024, 8, cur_block);
        std::cout << "cur_block: " << cur_block << ", time: " << (cur_time) << ", delta: " << (cur_time - prev_time) << "\n";
        fout << cur_block << ' ' << cur_time << '\n';
        jump_detector.add_measure(cur_time);
        if (jump_detector.is_jump()) {
            res_block = cur_block;
            break;
        }
        prev_time = cur_time;
    }
    fout.close();
    return res_block;
}

int main() {
    print_cpu_version();
    int res_assoc = assoc();
    std::cout << "Assoc: " << res_assoc << "\n";
    int res_size = size();
    std::cout << "Size: " << res_size << "\n";
    int res_block = block();
    std::cout << "Block size: " << res_block << "\n";
}
