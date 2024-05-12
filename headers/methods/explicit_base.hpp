#pragma once
#include <Eigen/Dense>
#include <deque>
#include <vector>

namespace methods {
template <int Dimensions> class ExplicitBase {

public:
  using Vec = Eigen::Vector<double, Dimensions>;
  using System = std::function<Vec(const double, Vec)>;
  Vec y_n;
  double delta = 0;
  std::deque<Vec> f_diffs;

  std::function<Vec(const double, const Vec)> system;
  std::vector<Vec> solution;

  double t_n;
  Vec _y0;

  virtual void compute_next_point() = 0;

  [[nodiscard]] virtual std::vector<Vec> compute(const size_t n) = 0;

  explicit ExplicitBase<Dimensions>(
      const double delta, const Vec y0,
      const std::function<Vec(const double, const Vec)> &system)
      : y_n{y0}, delta{delta}, _y0{y0}, system{system} {}

  virtual ~ExplicitBase<Dimensions>() {}
};
}; // namespace methods
