#pragma once

#include <Eigen/Dense>
#include <methods/explicit_base.hpp>
#include <vector>

namespace methods {
template <int Dimensions>
class ExplicitRGKutta final : private ExplicitBase<Dimensions> {
  using Vec = Eigen::Vector<double, Dimensions>;
  using System = std::function<Vec(const double, Vec)>;

  using ExplicitBase<Dimensions>::t_n;
  using ExplicitBase<Dimensions>::y_n;
  using ExplicitBase<Dimensions>::delta;
  using ExplicitBase<Dimensions>::f_diffs;
  using ExplicitBase<Dimensions>::system;
  using ExplicitBase<Dimensions>::solution;

  void compute_next_point() override {
    Vec sum = Vec::Zero();

    sum += delta * k1(t_n, y_n);
    sum += delta * k2(t_n, y_n);
    sum += delta * k3(t_n, y_n);
    sum += delta * k4(t_n, y_n);

    y_n += sum;

    t_n += delta;
  }

  Vec k1(const double t, const Vec y) { return system(t, y) / 6; }
  Vec k2(const double t, const Vec y) {
    return system(t + delta / 2, y + delta * k1(t, y) / 2) / 3;
  }
  Vec k3(const double t, const Vec y) {
    return system(t + delta / 2, y + delta * k2(t, y) / 2) / 3;
  }
  Vec k4(const double t, const Vec y) {
    return system(t + delta, y + delta * k3(t, y)) / 6;
  }

public:
  [[nodiscard]] std::deque<Vec> get_f_diffs() { return f_diffs; }

  [[nodiscard]] std::vector<Vec> compute(const size_t n) override {

    for (size_t i = 0; i < n; ++i) {
      f_diffs.push_back(system(t_n, y_n));
      compute_next_point();
      solution.push_back(y_n);
    }

    return solution;
  }

  explicit ExplicitRGKutta<Dimensions>(
      const double delta, const Vec y0,
      const std::function<Vec(const double, const Vec)> &system)
      : ExplicitBase<Dimensions>{delta, y0, system} {}

  ~ExplicitRGKutta<Dimensions>() override = default;
};

}; // namespace methods
