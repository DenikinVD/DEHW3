#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <algorithm>

struct Phase_dot {
    Phase_dot(double _x, double _y) : x(_x), y(_y) {}
    Phase_dot() : x(0), y(0) {}
    double x;
    double y;
};

void print_sol(std::ofstream &out, const std::vector<Phase_dot>& sol) {
    out << sol.size() << "\n";
    for (size_t j = 0; j != sol.size(); ++j) {
        out << std::fixed << std::setprecision(6) << sol[j].x << " "
        << std::fixed << std::setprecision(6) << sol[j].y << "\n";
    }
}

Phase_dot unlinsys(double x, double y) {
    double dx = sin(x) + exp(y) - 1;
    double dy = sin(x - y);
    return Phase_dot(dx, dy);
}

template<typename FUNC>
std::vector<Phase_dot> eyler_with_borders(Phase_dot init_pos,
                                          std::pair<double, double> x_borders,
                                          std::pair<double, double> y_borders,
                                          double step, FUNC ode) {
    std::vector<Phase_dot> res;
    Phase_dot pos = init_pos;
    int step_counter = 0;
    while (x_borders.first <= pos.x &&
           x_borders.second >= pos.x &&
           y_borders.first <= pos.y &&
           y_borders.second >= pos.y && step_counter < 10000) {
        res.push_back(pos);
        Phase_dot dv = ode(pos.x, pos.y);
        pos.x += dv.x * step;
        pos.y += dv.y * step;
        ++step_counter;
        //std::cout << pos.x << " " << pos.y << "\n";
    }

    return res;
}


std::vector<std::vector<Phase_dot>> separatrisses(std::pair<double, double> x_borders,
                                                  std::pair<double, double> y_borders, double step) {
    std::vector<std::vector<Phase_dot>> res(4);
    std::vector<Phase_dot> poses(4);
    std::vector<double> shift{-1 - sqrt(2), -1 + sqrt(2), -1 + sqrt(2), -1 - sqrt(2)};
    for (int i = 0; i != 4; ++i) {
        int counter = 0;
        poses[i] = Phase_dot(0, 0);
        while (x_borders.first <= poses[i].x &&
               x_borders.second >= poses[i].x &&
               y_borders.first <= poses[i].y &&
               y_borders.second >= poses[i].y && counter < 10000) {
            res[i].push_back(poses[i]);
            poses[i].y += shift[i] * step;
            poses[i].x += step;
            ++counter;
            shift[i] = sin(poses[i].x - poses[i].y) / (sin(poses[i].x) + exp(poses[i].y) - 1);
            /*if (fabs(sin(poses[i].x) + exp(poses[i].y) - 1) > 0.00001) {

            } else {
                break;
            }*/
        }
        step = -step;
    }
    return res;
}

int main() {
    std::ofstream out;
    out.open("HW3.txt");
    std::vector<std::vector<Phase_dot>> portraits;
    std::pair<double, double> x_borders(-5, 5);
    std::pair<double, double> y_borders(-5, 5);
    double len_step_x = (x_borders.second - x_borders.first) / 11;
    double len_step_y = (y_borders.second - y_borders.first) / 11;
    for (int i = 1; i != 11; ++i) {
        Phase_dot init_pos_1 (x_borders.first + static_cast<double>(i) * len_step_x,
                              y_borders.second - 0.00005);
        Phase_dot init_pos_2 (x_borders.first + static_cast<double>(i) * len_step_x,
                              y_borders.first + 0.00005);
        Phase_dot init_pos_3 (x_borders.first + 0.00005,
                              y_borders.first + static_cast<double>(i) * len_step_y);
        Phase_dot init_pos_4 (x_borders.second - 0.00005,
                              y_borders.first + static_cast<double>(i) * len_step_y);
        //init_pos_4.x << " " << init_pos_4.y << "\n";
        auto res1 = eyler_with_borders(init_pos_1, x_borders, y_borders, 0.005, unlinsys);
        auto res2 = eyler_with_borders(init_pos_2, x_borders, y_borders, 0.005, unlinsys);
        auto res3 = eyler_with_borders(init_pos_3, x_borders, y_borders, 0.005, unlinsys);
        auto res4 = eyler_with_borders(init_pos_4, x_borders, y_borders, 0.005, unlinsys);
        if (res1.size() > 100) {
            portraits.push_back(res1);
        }
        if (res2.size() > 100) {
            portraits.push_back(res2);
        }
        if (res3.size() > 100) {
            portraits.push_back(res3);
        }
        if (res4.size() > 100) {
            portraits.push_back(res4);
        }

    }
    Phase_dot last_dot(-5 + 0.00005, 3);
    auto res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(-5 + 0.00005, 3.5);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(-5 + 0.00005, 4.5);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(0, -1);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(0, -1);
    res = eyler_with_borders(last_dot, x_borders, y_borders, -0.0005, unlinsys);
    reverse(res.begin(), res.end());
    portraits.push_back(res);
    last_dot = Phase_dot(0, -2);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(0, -2);
    res = eyler_with_borders(last_dot, x_borders, y_borders, -0.0005, unlinsys);
    reverse(res.begin(), res.end());
    portraits.push_back(res);
    last_dot = Phase_dot(0, -3);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(0, -3);
    res = eyler_with_borders(last_dot, x_borders, y_borders, -0.0005, unlinsys);
    reverse(res.begin(), res.end());
    portraits.push_back(res);
    last_dot = Phase_dot(0, -4);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(0, -4);
    res = eyler_with_borders(last_dot, x_borders, y_borders, -0.0005, unlinsys);
    reverse(res.begin(), res.end());
    portraits.push_back(res);
    last_dot = Phase_dot(2, 0);
    res = eyler_with_borders(last_dot, x_borders, y_borders, 0.0005, unlinsys);
    portraits.push_back(res);
    last_dot = Phase_dot(2, 0);
    res = eyler_with_borders(last_dot, x_borders, y_borders, -0.0005, unlinsys);
    reverse(res.begin(), res.end());
    portraits.push_back(res);
    out << portraits.size() << "\n";
    for (auto& port : portraits) {
        print_sol(out, port);
    }
    out.close();
    auto separ = separatrisses(x_borders, y_borders, 0.0005);
    out.open("HW3_separat.txt");
    out << "4\n";
    for (auto& sep : separ) {
        print_sol(out, sep);
    }
    out.close();
}
