#include <cpuid.h>

#include <chrono>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "timer.h"

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

double measure(char* array, int bench_cnt, int accesses_cnt) {
    timer timer;

    // double time = 0;
    for (int bench = 0; bench < bench_cnt; bench++) {
        // reset_cache();
        timer.restart();
        // auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < accesses_cnt; i++) {
            make_accesses(array);
        }
        // auto end = std::chrono::high_resolution_clock::now();
        // time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        timer.nextLap();
    }
    return timer.lapAvg();
    // return time / bench_cnt;
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
    auto pos = create_pos(0, cur_size / 16, 16);
    create_accesses_array(memory, pos);
    return measure(memory, 10, 10000);
}

double measure_assoc(int cache_size, int cur_assoc) {
    auto pos = create_pos(0, cur_assoc, cache_size);
    create_accesses_array(memory, pos);
    return measure(memory, 10, 10000000);
}

double measure_block(int cache_size, int assoc, int cur_size) {
    int way_size = cache_size / assoc;

    auto pos1 = create_pos(0, assoc/2, way_size);
    auto pos2 = create_pos(way_size*(assoc/2) + cur_size, assoc/2 + 1, way_size);
    pos1.insert(pos1.end(), pos2.begin(), pos2.end());
    create_accesses_array(memory, pos1);
    return measure(memory, 50, 1000000);
}

// double measure_ws_and_assoc() {

//     auto pos = create_pos(0, cur_assoc, cache_size);
//     create_accesses_array(memory, pos);
//     return measure(memory);
// }

int main() {
    double prev_time = 0, prev_delta = 0;
    std::ofstream fout("measure_assoc");

    double max_delta = INT32_MIN;
    // int res_assoc = 0;
    // for (int cur_assoc = 1; cur_assoc <= 16; cur_assoc++) {
    //     auto cur_time = measure_assoc(32 * 1024, cur_assoc);
    //     if (prev_time == 0)
    //         prev_time = cur_time;
    //     double delta = cur_time - prev_time;
    //     if (prev_delta == 0)
    //         prev_delta = delta;
    //     if (delta - prev_delta > max_delta) {
    //         max_delta = delta - prev_delta;
    //         res_assoc = cur_assoc - 1;
    //     }
    //     std::cout << "cur_assoc: " << cur_assoc << ", time: " << (cur_time) << ", delta: " << delta << "\n";
    //     fout << cur_assoc << ' ' << cur_time << '\n';
    //     prev_time = cur_time;
    //     prev_delta = delta;
    // }
    // fout.close();

    // max_delta = INT32_MIN;
    // int res_size = 0;
    // prev_delta = 0;
    // prev_time = 0;
    // fout = std::ofstream("measure_size");
    // for (int cur_size = 1; cur_size <= 64; cur_size++) {
    //     auto cur_time = measure_size(cur_size * 1024);
    //     if (prev_time == 0)
    //         prev_time = cur_time;
    //     double delta = cur_time - prev_time;
    //     if (prev_delta == 0)
    //         prev_delta = delta;
    //     if (delta - prev_delta > max_delta) {
    //         max_delta = delta - prev_delta;
    //         res_size = cur_size - 1;
    //     }
    //     std::cout << "cur_size: " << cur_size << ", time: " << (cur_time) << ", delta: " << delta << "\n";
    //     fout << cur_size << ' ' << cur_time << '\n';
    //     prev_time = cur_time;
    //     prev_delta = delta;
    // }
    // fout.close();

    max_delta = INT32_MIN;
    int res_block = 0;
    prev_delta = 0;
    prev_time = 0;
    fout = std::ofstream("measure_block");
    for (int cur_block = 16; cur_block <= 128; cur_block+=8) {
        auto cur_time = measure_block(32 * 1024, 8, cur_block);
        if (prev_time == 0)
            prev_time = cur_time;
        double delta = cur_time - prev_time;
        if (prev_delta == 0)
            prev_delta = delta;
        if (delta - prev_delta > max_delta) {
            max_delta = delta - prev_delta;
            res_block = cur_block - 1;
        }
        std::cout << "cur_block: " << cur_block << ", time: " << (cur_time) << ", delta: " << delta << "\n";
        fout << cur_block << ' ' << cur_time << '\n';
        prev_time = cur_time;
        prev_delta = delta;
    }
    fout.close();

    // std::cout << "Res assoc: " << res_assoc << ", res size: " << res_size;
}
