#include "pocketfft_hdronly.h"
#include <cmath>
#include <complex>
#include <matplot/matplot.h>
#include <vector>

#include <iostream>

template <typename T> void vec_sin(std::vector<std::complex<T>> &v) {
  for (int i = 0; i < v.size(); i++)
    v[i] = std::complex<T>(std::sin((float)i) - 1.0, 0.f);
}

template <typename T>
void fftc2c(std::vector<std::complex<T>> const &data_in,
            std::vector<std::complex<T>> &data_out) {
  auto len = data_in.size();
  pocketfft::shape_t shape{len};
  pocketfft::stride_t stridef(shape.size());
  size_t tmpf = sizeof(std::complex<T>);
  for (int i = shape.size() - 1; i >= 0; --i) {
    stridef[i] = tmpf;
    tmpf *= shape[i];
  }
  pocketfft::shape_t axes;
  for (size_t i = 0; i < shape.size(); ++i) {
    axes.push_back(i);
  }
  pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::FORWARD,
                 data_in.data(), data_out.data(), (T)1.0 / (T)len);
}

template <typename T>
void ifftc2c(std::vector<std::complex<T>> const &data_in,
             std::vector<std::complex<T>> &data_out) {
  auto len = data_in.size();
  pocketfft::shape_t shape{len};
  pocketfft::stride_t stridef(shape.size());
  size_t tmpf = sizeof(std::complex<T>);
  for (int i = shape.size() - 1; i >= 0; --i) {
    stridef[i] = tmpf;
    tmpf *= shape[i];
  }
  pocketfft::shape_t axes;
  for (size_t i = 0; i < shape.size(); ++i) {
    axes.push_back(i);
  }
  pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::BACKWARD,
                 data_in.data(), data_out.data(), static_cast<T>(1));
}

int main() {
  size_t len = 101;

  std::vector<std::complex<float>> data(len);
  vec_sin(data);

  std::vector<std::complex<float>> res = data;

  std::vector<float> plotting_vector_first_real;
  std::vector<float> plotting_vector_first_imag;
  for (auto i = 0; i < 100; i++) {
    plotting_vector_first_real.push_back(data[i].real());
    plotting_vector_first_imag.push_back(data[i].imag());
  }

  matplot::subplot(2, 2, 0);
  matplot::plot(plotting_vector_first_real);
  matplot::hold(matplot::on);
  matplot::plot(plotting_vector_first_imag);

  fftc2c(data, res);

  std::vector<float> plotting_vector_second_real{};
  std::vector<float> plotting_vector_second_imag{};
  for (auto i = 0; i < 100; i++) {
    plotting_vector_second_real.push_back(0.);
    plotting_vector_second_imag.push_back(0.);
  }

  for (auto i = 0; i < res.size(); i++) {
    plotting_vector_second_real[i] = res[i].real();
    plotting_vector_second_imag[i] = res[i].imag();
  }
  matplot::subplot(2, 2, 1);
  matplot::plot(plotting_vector_second_real);
  matplot::hold(matplot::on);
  matplot::plot(plotting_vector_second_imag);
  matplot::hold(matplot::off);

  ifftc2c(res, data);

  std::vector<float> plotting_vector_last_real;
  std::vector<float> plotting_vector_last_imag;

  for (auto i = 0; i < 100; i++) {
    plotting_vector_last_real.push_back(data[i].real());
    plotting_vector_last_imag.push_back(data[i].imag());
  }
  matplot::subplot(2, 2, 2);
  matplot::plot(plotting_vector_last_real);
  matplot::hold(matplot::on);
  matplot::plot(plotting_vector_last_imag);
  matplot::hold(matplot::off);
  matplot::show();
  return 0;
}