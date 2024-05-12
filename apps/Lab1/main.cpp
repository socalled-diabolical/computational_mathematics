#include "SFML/Graphics/CircleShape.hpp"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <limits>
#include <methods/explicit_adams.hpp>
#include <methods/explicit_bdf.hpp>
#include <methods/explicit_rgkutta.hpp>
#define XSIZE 1920
#define YSIZE 1080

using Vec2 = Eigen::Vector<double, 2>;
Vec2 order1odu(const double t, const Vec2 vec);
std::vector<Vec2> center_and_resize(const std::vector<Vec2> &vec);

int main() {
  Vec2 x0{2, 0};
  std::function<Vec2(const double, const Vec2)> system = order1odu;

  methods::ExplicitBDF<2> method(1e-4, x0, system);
  auto sol = method.compute(100 / 1e-4);

  sol = center_and_resize(sol);
  sf::RenderWindow window(sf::VideoMode(XSIZE, YSIZE), "AAAAAAAAAAAAAAAA");
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed)
        window.close();
    }
    for (auto &&x : sol) {
      sf::CircleShape shape(5);
      shape.setPosition(x[0], x[1]);

      window.draw(shape);
    }
    window.display();
  }

  return 0;
}

Vec2 order1odu(const double t, const Vec2 vec) {
  double param = 1;
  return Vec2(vec[1], param * (1 - vec[0] * vec[0]) * vec[1] - vec[0]);
}

std::vector<Vec2> center_and_resize(const std::vector<Vec2> &vec) {
  std::vector<Vec2> out_vec;

  double x_max = std::numeric_limits<double>::min();
  double y_max = std::numeric_limits<double>::min();

  for (auto &&x : vec) {
    if (x_max < x[0])
      x_max = x[0];

    if (y_max < x[1])
      y_max = x[1];
  }

  double x_coef = XSIZE / (4 * x_max);
  double y_coef = YSIZE / (4 * y_max);

  Vec2 temp = Vec2::Zero();

  for (auto &&x : vec) {
    temp[0] = x_coef * x[0] + static_cast<double>(XSIZE) / 2;
    temp[1] = y_coef * x[1] + static_cast<double>(YSIZE) / 2;
    out_vec.push_back(temp);
  }
  return out_vec;
}
