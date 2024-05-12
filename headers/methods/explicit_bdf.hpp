#pragma once
#include <methods/explicit_rgkutta.hpp>
namespace methods {
template <int Dimensions>
class ExplicitBDF final : public ExplicitBase<Dimensions> {
private:
  using Vec = Eigen::Vector<double, Dimensions>;
  using System = std::function<Vec(const double, Vec)>;
  std::deque<Vec> y_ns;

  using ExplicitBase<Dimensions>::t_n;
  using ExplicitBase<Dimensions>::y_n;
  using ExplicitBase<Dimensions>::delta;
  using ExplicitBase<Dimensions>::f_diffs;
  using ExplicitBase<Dimensions>::system;
  using ExplicitBase<Dimensions>::solution;

  void compute_next_point() override {
    Vec y_old;
    y_n = y_ns[3];

    for (size_t i = 0; i < 100; ++i) {
      y_old = y_n;
      y_n = (48.0 / 25) * y_ns[3] - (36.0 / 25) * y_ns[2] +
            (16.0 / 25) * y_ns[1] - (3.0 / 25) * y_ns[0] +
            (12.0 / 25) * delta * system(t_n, y_n);

      Vec diff = y_n - y_old;
      if (std::abs(diff.norm()) < 1e-8)
        break;
    }

    t_n += delta;
  }

public:
  [[nodiscard]] std::vector<Vec> compute(const size_t n) override {

    for (size_t i = 0; i < n - 4; ++i) {
      compute_next_point();
      solution.push_back(y_n);
      y_ns.pop_front();
      y_ns.push_back(y_n);
    }

    return solution;
  }

  explicit ExplicitBDF<Dimensions>(
      const double delta, const Vec y0,
      const std::function<Vec(const double, const Vec)> &system)
      : ExplicitBase<Dimensions>(delta, y0, system) {

    ExplicitRGKutta<Dimensions> rgkutta_solver(delta, y0, system);

    solution = rgkutta_solver.compute(4);
    for (auto &&x : solution) {
      y_ns.push_back(x);
    }
    y_n = solution.back();

    t_n = delta * 4;
  }

  ~ExplicitBDF<Dimensions>() override = default;
};
}; // namespace methods
