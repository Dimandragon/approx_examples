#include <icecream.hpp>

#include <cmath>
#include <makima.hpp>
#include <optional>
#include <vector>
#include <matplot/matplot.h>

struct SquarePolynome {
  double a;
  double b;
  double c;

  void solve(auto x1, auto x2, auto x3, auto y1, auto y2, auto y3) {
    b = ((y1 - y3) * (x1 * x1 - x2 * x2) - (y1 - y2) * (x1 * x1 - x3 * x3)) /
        ((x1 - x2) * (x1 - x3) * (x2 - x3));
    a = (y1 - y2 - b * (x1 - x2)) / (x1 * x1 - x2 * x2);
    c = y3 - a * x3 * x3 - b * x3;
  }

  double compute(auto x) { return a * x * x + b * x + c; }

  double derive(auto x) { return b + 2.0 * a; }
};

struct Linear {
  double k;
  double b;

  void solve(auto x1, auto x2, auto y1, auto y2) {
    auto dx = x2 - x1;
    auto dy = y2 - y1;

    k = dy / dx;
    b = x1 - k * x1;
  }

  double compute(auto x) { return k * x + b; }
  double derive(auto x) { return k; }
};

// Интерполяция модифицированным сплайном Акимы
// see https://www.mathworks.com/help/matlab/ref/makima.html
template <typename T> struct ModifiedAkimaBasedWithNoTrain {
  // using boost::math::interpolators::makima;
  std::optional<boost::math::interpolators::makima<std::vector<double>>>
      spline = {};
  std::optional<SquarePolynome> square_polynom = {};
  std::optional<Linear> linear = {};

  double min_bound = 0.0;
  double max_bound = 0.0;

  void loadData(const T &x, const T &y) {
    std::vector<double> x_(x.size());
    std::vector<double> y_(x.size());
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = x[i];
      x_[i] = el;
    }
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    double min = 9999999999999999999999.0;  // todo
    double max = -9999999999999999999999.0; // todo
    for (auto const &el : x_) {
      if (el < min) {
        min = el;
      }
      if (el > max) {
        max = el;
      }
    }
    min_bound = min;
    max_bound = max;

    if (x_.size() > 3) {
      spline = boost::math::interpolators::makima<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  void loadData(const T &y) {
    // IC(y.size());
    std::vector<double> x_(y.size());
    std::vector<double> y_(y.size());
    for (int i = 0; i < y.size(); i++) {
      x_[i] = (double)i;
    }
    for (int i = 0; i < y.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    // IC(y_.size());
    min_bound = 0.0;
    max_bound = y.size() - 1;

    if (x_.size() > 3) {
      spline = boost::math::interpolators::makima<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  template <typename IdxT> double compute(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return (*spline)(idx);
      } else if (square_polynom) {
        return square_polynom->compute(idx);
      } else if (linear) {
        return linear->compute(idx);
      } else {
        return 0.0;
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size);
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size);
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
      }
    }
  }

  template <typename IdxT> double computeDerive(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return spline->prime(idx);
      } else if (square_polynom) {
        return square_polynom->derive(idx);
      } else if (linear) {
        return linear->derive(idx);
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }

  void show(double step, double first, double last, bool plot_original = false) {
    using SampleType = double;
    std::vector<SampleType> plotting_data = {};
    for (auto i = first; i < last; i += step) {
      plotting_data.push_back(compute<double>(i));
    }
    if (plot_original) {
      std::vector<SampleType> plotting_data_original = {};
      for (auto i = first; i < last; i += step) {
        plotting_data_original.push_back(compute<int>(i));
      }
      matplot::plot(plotting_data);
      matplot::hold(true);
      matplot::plot(plotting_data_original);
      matplot::hold(false);
      matplot::show();
      return;
    }
    matplot::plot(plotting_data);
    matplot::show();
  }
};


int main() {
  ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;
  auto size = 100;

  std::vector<double> y_data(size);
  std::vector<double> x_data(size);
  for (int i = 0; i < size; i++) {
    y_data[i] = std::rand();
    x_data[i] = i;
  }
  auto avg = 0.0;
  for (int i = 0; i < y_data.size(); i++) {
    avg += y_data[i];
  }
  for (int i = 0; i < y_data.size(); i++) {
    y_data[i] /= avg;
  }

  approximator.loadData(x_data, y_data);

  approximator.show(1.0, 0, size);
  approximator.show(0.1, 0, size, true);

  return 0;
}