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

class regression {
    std::vector<double> x;
    std::vector<double> y;
    double coeff{0};
    double constTerm{0};
    double sum_xy{0};
    double sum_x{0};
    double sum_y{0};
    double sum_x_square{0};
    double sum_y_square{0};

   public:
    void calculateCoefficient() {
        double N = x.size();
        double numerator = (N * sum_xy - sum_x * sum_y);
        double denominator = (N * sum_x_square - sum_x * sum_x);
        coeff = numerator / denominator;
    }

    void calculateConstantTerm() {
        double N = x.size();
        double numerator = (sum_y * sum_x_square - sum_x * sum_xy);
        double denominator = (N * sum_x_square - sum_x * sum_x);
        constTerm = numerator / denominator;
    }

    void addValue(double xi, double yi) {
        sum_xy += xi * yi;
        sum_x += xi;
        sum_y += yi;
        sum_x_square += xi * xi;
        sum_y_square += yi * yi;
        x.push_back(xi);
        y.push_back(yi);
    }

    double predict(double x) {
        calculateCoefficient();
        calculateConstantTerm();
        return coeff * x + constTerm;
    }
};

size_t max_size = 1 << 25;
int alignment = 1 << 12;  // 4k
char* memory = (char*)aligned_alloc(alignment, max_size);
auto rnd = std::default_random_engine{};
std::ofstream log_file("log.txt");

inline void make_accesses(char* cur_addr) {
    while (cur_addr != nullptr) {
        cur_addr = *(char**)cur_addr;
    }
}

void create_accesses_array(char* array, std::vector<size_t>& pos) {
    if (pos.size() > 1) {
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

double measure_accesses(char* array, int accesses_cnt) {
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, std::micro> duration;

        make_accesses(array);
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
    return sum / bench_cnt;
}


double measure_size(int cur_size) {
    return bench(10000, [&]() {
        auto pos = create_pos(0, cur_size / 16, 16);
        create_accesses_array(memory, pos);
        return measure_accesses(memory, 100);
    });
}

double measure_assoc(int cache_size, int cur_assoc) {
    return bench(1000000, [&]() {
        auto pos = create_pos(0, cur_assoc, cache_size);
        create_accesses_array(memory, pos);
        return measure_accesses(memory, 500);
    });
}

double measure_block(int cache_size, int assoc, int cur_size) {
    int way_size = cache_size / assoc;
    return bench(1000000, [&]() {
        auto pos1 = create_pos(0, assoc, way_size);
        auto pos2 = create_pos(way_size * assoc + cur_size, assoc, way_size);
        pos1.insert(pos1.end(), pos2.begin(), pos2.end());
        create_accesses_array(memory, pos1);
        return measure_accesses(memory, 100);
    });
}

class JumpDetector {
   public:
    explicit JumpDetector(double jump_coef) : jump_coef(jump_coef) {}
    void add_measure(int x, double measure) {
        double predicted = r.predict(x);
        log_file << "    predict: " << predicted << " actual: " << measure << " diff: " << measure - predicted << std::endl;
        r.addValue(x, measure);

        if (measures_count >= 2) {
            double error_coef = (measure - predicted) / predicted;
            log_file << "    error coef: " << error_coef << std::endl;
            if (jump_coef > 0 && error_coef > jump_coef || jump_coef < 0 && error_coef < jump_coef)
                is_jump_ = true;
        }
        measures_count++;
    }

    bool is_jump() const { return is_jump_; }

   private:
    regression r;
    bool is_jump_ = false;
    double jump_coef;
    int measures_count = 0;
};

std::pair<int, int> get_size_and_assoc() {
    int ws = 8 * 1024;
    int assoc_old = 0;
    while (true) {
        log_file << "ws: " << ws << std::endl;
        JumpDetector jump_detector(0.1);
        for (int assoc = 1;; assoc++) {
            log_file << "assoc: " << assoc << std::endl;
            double time = bench(1000000, [&]() {
                auto pos = create_pos(0, assoc, ws);
                create_accesses_array(memory, pos);
                return measure_accesses(memory, 100);
            });
            jump_detector.add_measure(assoc, time);
            if (jump_detector.is_jump() && assoc > 3) {
                if (assoc - 1 == 2 * assoc_old)
                    return {ws * 2 * assoc_old, assoc_old};
                assoc_old = assoc - 1;
                break;
            }
        }
        ws /= 2;
    }
}

int get_assoc(int cache_size) {
    std::ofstream fout("measure_assoc");
    int res_assoc = 0;
    JumpDetector jump_detector(0.2);
    for (int cur_assoc = 1; cur_assoc <= 32; cur_assoc++) {
        auto cur_time = measure_assoc(cache_size, cur_assoc);
        log_file << "cur_assoc: " << cur_assoc << std::endl;
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

int get_size() {
    int res_size = 0;
    std::ofstream fout("measure_size");
    JumpDetector jump_detector(0.3);
    for (int cur_size = 1; cur_size <= 256; cur_size *= 2) {
        auto cur_time = measure_size(cur_size * 1024);
        log_file << "cur_size: " << cur_size << std::endl;
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

int get_block(int cache_size, int assoc) {
    int res_block = 0;
    std::ofstream fout("measure_block");
    JumpDetector jump_detector(0.4);
    for (int cur_block = 2 * 1024; cur_block >= 2; cur_block /= 2) {
        auto cur_time = measure_block(cache_size, assoc, cur_block);
        log_file << "cur_block: " << cur_block << std::endl;
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
    std::cout << "Calculating...\n";
    // auto [s, a] = get_size_and_assoc();
    // std::cout << s << ' ' << a;
    // return 0;
    int res_size = get_size();
    if (res_size == 0) {
        std::cout << "Can't find out cache size...";
        return 1;
    }
    std::cout << "Size: " << res_size << "B\n";

    int res_assoc = get_assoc(res_size);
    if (res_assoc == 0) {
        std::cout << "Can't find out cache associativity...";
        return 1;
    }
    std::cout << "Assoc: " << res_assoc << std::endl;
    
    int res_block = get_block(res_size, res_assoc);
    if (res_block == 0) {
        std::cout << "Can't find out cache block size...";
        return 1;
    }
    std::cout << "Block size: " << res_block << "B\n";
}
