#pragma once
#include <iostream>
#include <methods/explicit_rgkutta.hpp>
namespace methods {
template <int Dimensions>
class ExplicitAdams final : private ExplicitBase<Dimensions> {
  using Vec = Eigen::Vector<double, Dimensions>;
  using System = std::function<Vec(const double, Vec)>;

  using ExplicitBase<Dimensions>::t_n;
  using ExplicitBase<Dimensions>::y_n;
  using ExplicitBase<Dimensions>::delta;
  using ExplicitBase<Dimensions>::f_diffs;
  using ExplicitBase<Dimensions>::system;
  using ExplicitBase<Dimensions>::solution;

  void compute_f_diffs() {
    f_diffs.pop_front();
    f_diffs.push_back(system(t_n, y_n));
  }

  void compute_next_point() override {
    compute_f_diffs();
    y_n += delta * ((55 / 24) * f_diffs[0] - (59 / 24) * f_diffs[1] +
                    (37 / 24) * f_diffs[2] - (9 / 24) * f_diffs[3]);
    t_n += delta;
  }

public:
  [[nodiscard]] std::vector<Vec> compute(const size_t n) override {

    for (size_t i = 0; i < n - 4; ++i) {
      compute_next_point();
      solution.push_back(y_n);
    }

    return solution;
  }

  explicit ExplicitAdams<Dimensions>(
      const double delta, const Vec y0,
      const std::function<Vec(const double, const Vec)> &system)
      : ExplicitBase<Dimensions>(delta, y0, system) {

    ExplicitRGKutta<Dimensions> rgkutta_solver(delta, y0, system);

    auto sol = rgkutta_solver.compute(4);
    f_diffs = rgkutta_solver.get_f_diffs();

    for (auto &&x : sol) {
      std::cout << x << std::endl;
    }

    y_n = sol.back();

    t_n = delta * 4;
  }

  ~ExplicitAdams<Dimensions>() override = default;
};
}; // namespace methods
