#pragma once
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
