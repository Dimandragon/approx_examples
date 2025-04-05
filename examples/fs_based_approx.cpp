#include <icecream.hpp>

#include "pocketfft_hdronly.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <icecream.hpp>
#include <iostream>
#include <matplot/matplot.h>
#include <string>
#include <vector>

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

template <typename T>
void fftc2c(std::vector<std::complex<T>> const &data_in,
            std::vector<std::complex<T>> &data_out, size_t len, size_t pad) {
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
                 data_in.data() + pad, data_out.data() + pad,
                 ((T)1.0) / ((T)len));
}

template <typename T>
void ifftc2c(std::vector<std::complex<T>> const &data_in,
             std::vector<std::complex<T>> &data_out, size_t len, size_t pad) {
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
                 data_in.data() + pad, data_out.data() + pad,
                 static_cast<T>(1));
}

// Аппроксимация рядом фурье
struct FourierSeriesBasedApprox {
  // static constexpr bool is_signal_approximator = true;
  size_t signal_size;
  bool is_actual = false;
  int polynoms_count;

  int polynoms_count_on_tile;
  int tile_size;

  std::vector<std::complex<double>> fourier_series;
  std::vector<std::complex<double>> approximated_data;

  template <typename SignalType>
  FourierSeriesBasedApprox(const SignalType &signal_in) {
    signal_size = signal_in.size();
    fourier_series = std::vector<std::complex<double>>(signal_in.size());
    approximated_data = std::vector<std::complex<double>>(signal_in.size());
    polynoms_count = signal_in.size();
    polynoms_count_on_tile = signal_in.size();
    tile_size = signal_in.size();
  }

  template <typename SignalType> void computeFSFromData(SignalType &signal_in) {
    for (size_t i = 0; i < signal_in.size(); i++) {
      approximated_data[i] = {signal_in[i], 0.0};
    }
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      fftc2c(approximated_data, fourier_series, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    fftc2c(approximated_data, fourier_series, approximated_data.size() - pad,
           pad);
  }

  size_t mirrorIdx(size_t const i) const {
    return i / tile_size * tile_size + i / tile_size * tile_size + tile_size -
           i;
  }

  void applyMirror(size_t const i) {
    if (i % tile_size != 0) {
      auto const mirror_idx = mirrorIdx(i);
      fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0,
                                    -fourier_series[i].imag() / 2.0};
      fourier_series[i] = fourier_series[i] / 2.0;
    }
  }

  void computeData() {
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      ifftc2c(fourier_series, approximated_data, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad,
            pad);
  }

  void computeTile(size_t const idx) {
    auto const i = idx / tile_size;
    auto const pad = i * tile_size;
    if (approximated_data.size() > pad + tile_size) {
      ifftc2c(fourier_series, approximated_data, tile_size, pad);
    } else {
      ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad,
              pad);
    }
  }

  void setApproxOrderRatio(double const ratio) {
    polynoms_count_on_tile =
        static_cast<int>(static_cast<double>(tile_size) * 0.5 * ratio);
  }

  template <typename IdxT> std::complex<double> computeSample(IdxT idx) {
    using SampleType = double;
    std::complex<double> accum = {0.0, 0.0};

    if (idx <= approximated_data.size() && idx >= 0) {
      auto tile_first_idx = static_cast<int>(idx) / tile_size * tile_size;
      if (tile_first_idx + tile_size < approximated_data.size()) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }

      } else {
        SampleType w =
            std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
            static_cast<SampleType>(approximated_data.size() - tile_first_idx);
        for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
      std::cout << accum << std::endl;
    } else if (idx > approximated_data.size()) {
      auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
      SampleType w =
          std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
          static_cast<SampleType>(approximated_data.size() - tile_first_idx);
      for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
        accum += fourier_series[i + tile_first_idx] *
                 std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
      }
    } else {
      if (approximated_data.size() > tile_size) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      } else {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(approximated_data.size());
        for (auto i = 0; i < approximated_data.size(); i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
    }

    return accum;
  }

  template <typename SampleType, typename IdxType>
  SampleType compute(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx].real();
      } else {
        computeData();
        return approximated_data[idx].real();
      }
    } else {
      return computeSample(idx).real();
    }
  }

  template <typename SampleType, typename IdxType>
  std::complex<SampleType> computeComplex(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx];
      } else {
        computeData();
        return approximated_data[idx];
      }
    } else {
      return computeSample(idx);
    }
  }

  void show(double step, bool plot_original = false) {
    using SampleType = double;
    is_actual = false;
    std::vector<SampleType> plotting_data = {};
    for (auto i = 0.0; i < signal_size; i += step) {
      plotting_data.push_back(computeSample<double>(i).real());
    }
    if (plot_original) {
      std::vector<SampleType> plotting_data_original = {};
      for (auto i = 0.0; i < signal_size; i += step) {
        plotting_data_original.push_back(computeSample<int>(i).real());
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

void createFill(auto &signal) {
  for (auto i = 0; i < signal.size(); i++) {
    signal[i] = std::rand();
  }
}

void plotComplex(auto &signal) {
  std::vector<double> real;
  std::vector<double> imag;
  for (auto i = 0; i < signal.size(); i++) {
    real.push_back(signal[i].real());
    imag.push_back(signal[i].imag());
  }
  matplot::plot(real);
  matplot::hold(matplot::on);
  matplot::plot(imag, "--");
  matplot::hold(matplot::off);
  matplot::show();
}

int main() {
  int size = 102;
  std::vector<double> signal1;
  using SignalT = decltype(signal1);
  // if (len == 100) fs[51].re == fs[49].re, fs[51].im = -fs[49].im; fs[50]is
  // uniq if len == 101 its 50 and 51 0 and last are different if 10 then 4 and
  // 6 if 11 then 5 and 6
  for (auto i = 0; i < size; i++) {
    signal1.push_back(0.);
  }
  createFill(signal1);

  matplot::plot(signal1);
  matplot::show();

  SignalT signal2;

  FourierSeriesBasedApprox approximator(signal1);

  approximator.setApproxOrderRatio(1.);
  approximator.tile_size = size;
  approximator.computeFSFromData(signal1);

  plotComplex(approximator.fourier_series);

  approximator.show(1.0);
  approximator.show(0.1, true);

  return 0;
}