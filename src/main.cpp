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

#include "reg.h"

size_t max_size = 1 << 25;
int alignment = 1 << 12;  // 4k
char* memory = (char*)aligned_alloc(alignment, max_size);

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

float measure(char* array, int accesses_cnt) {
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<float, std::milli> duration;

    clock::time_point start = clock::now();
    for (int i = 0; i < accesses_cnt; i++) {
        make_accesses(array);
    }
    duration elapsed = clock::now() - start;
    return elapsed.count();
}

template <class Func>
double bench(int bench_cnt, Func one_measure_f) {
    double sum = 0;
    for (int bench = 0; bench < bench_cnt; bench++) {
        sum += one_measure_f();
    }
    return sum;
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

double measure_size(int cur_size) {
    return bench(25, [&]() {
        auto pos = create_pos(0, cur_size / 16, 16);
        create_accesses_array(memory, pos);
        return measure(memory, 100000);
    });
}

double measure_assoc(int cache_size, int cur_assoc) {
    return bench(25, [&]() {
        auto pos = create_pos(0, cur_assoc, cache_size);
        create_accesses_array(memory, pos);
        return measure(memory, 50000000);
    });
}

double measure_block(int cache_size, int assoc, int cur_size) {
    double sum = 0;
    int way_size = cache_size / assoc;
    return bench(25, [&]() {
        auto pos1 = create_pos(0, assoc / 2, way_size);
        auto pos2 = create_pos(way_size * (assoc / 2) + cur_size, assoc / 2 + 1, way_size);
        pos1.insert(pos1.end(), pos2.begin(), pos2.end());
        create_accesses_array(memory, pos1);
        return measure(memory, 50000000);
    });
}

// double measure_ws_and_assoc() {

//     auto pos = create_pos(0, cur_assoc, cache_size);
//     create_accesses_array(memory, pos);
//     return measure(memory);
// }

class JumpDetector {
   public:
    explicit JumpDetector(double jump_coef) : jump_coef(jump_coef) {}
    void add_measure(int x, double measure) {
        float predicted = r.predict(x);
        std::cout << "    predict: " << predicted << " actual: " << measure << " diff: " << measure - predicted << "\n";
        r.addValue(x, measure);

        if (measures_count >= 2) {
            double error_coef = (measure - predicted) / predicted;
            std::cout << "    error coef: " << error_coef << "\n";
            if (jump_coef > 0 && error_coef > jump_coef || jump_coef < 0 && error_coef < jump_coef)
                is_jump_ = true;
        }
        measures_count++;
        prev_measure = measure;
    }

    bool is_jump() const { return is_jump_; }

   private:
    regression r;
    bool is_jump_ = false;
    double jump_coef;
    int measures_count = 0;
    double delta_sum = 0;
    double prev_measure = 0;
};

int assoc(int cache_size) {
    std::ofstream fout("measure_assoc");
    int res_assoc = 0;
    JumpDetector jump_detector(0.5);
    for (int cur_assoc = 1; cur_assoc <= 32; cur_assoc++) {
        auto cur_time = measure_assoc(cache_size, cur_assoc);
        std::cout << "cur_assoc: " << cur_assoc << "\n";
        fout << cur_assoc << ' ' << cur_time << '\n';
        jump_detector.add_measure(cur_assoc, cur_time);
        if (jump_detector.is_jump()) {
            res_assoc = cur_assoc - 1;
            break;
        }
    }
    fout.close();
    return res_assoc;
}

int size() {
    int res_size = 0;
    std::ofstream fout("measure_size");
    JumpDetector jump_detector(0.5);
    for (int cur_size = 1; cur_size <= 256; cur_size *= 2) {
        auto cur_time = measure_size(cur_size * 1024);
        std::cout << "cur_size: " << cur_size << "\n";
        fout << cur_size << ' ' << cur_time << '\n';

        jump_detector.add_measure(cur_size, cur_time);
        if (jump_detector.is_jump()) {
            res_size = cur_size / 2;
            break;
        }
    }
    fout.close();
    return res_size * 1024;
}

int block(int cache_size, int assoc) {
    int res_block = 0;
    std::ofstream fout("measure_block");
    JumpDetector jump_detector(0.5);
    for (int cur_block = 2 * 1024; cur_block >= 2; cur_block /= 2) {
        auto cur_time = measure_block(cache_size, assoc, cur_block);
        std::cout << "cur_block: " << cur_block << "\n";
        fout << cur_block << ' ' << cur_time << '\n';
        jump_detector.add_measure(cur_block, cur_time);
        if (jump_detector.is_jump()) {
            res_block = cur_block * 2;
            break;
        }
    }
    fout.close();
    return res_block;
}

int main() {
    print_cpu_version();
    int res_size = size();
    std::cout << "Size: " << res_size << "\n";
    int res_assoc = assoc(res_size);
    std::cout << "Assoc: " << res_assoc << "\n";
    int res_block = block(res_size, res_assoc);
    std::cout << "Block size: " << res_block << "\n";

    std::cout << "Size: " << res_size << "B, Assoc: " << res_assoc << ", Block size: " << res_block << "B\n";
}
