#include "icecream.hpp"
#include <cmath>
#include <vector>
#include "interpolation.h"
#include <matplot/matplot.h>

enum class RBFKind { TPS, Gaussian, Bell, Multiquadric, MultiquadricAuto };
enum class LinTermKind { None, Const, Linear };
struct RBFBasedWithNoTrain {
  alglib::rbfmodel model = alglib::rbfmodel();
  alglib::rbfreport report = alglib::rbfreport();

  LinTermKind linterm_kind = LinTermKind::Linear;

  int n_dims_x = 1;
  int n_dims_y = 1;

  double lambda_v = 0.0;
  /*
  Only for TPS, Multiquadric or MultiquadricAuto

  lambda_v -   smoothing parameter, LambdaV>=0, defaults to 0.0:
          * LambdaV=0 means that no smoothing is applied,  i.e.  the
            spline tries to pass through all dataset points exactly
          * LambdaV>0 means that a smoothing thin  plate  spline  is
            built, with larger LambdaV corresponding to models  with
            less nonlinearities. Smoothing spline reproduces  target
            values at nodes with small error; from the  other  side,
            it is much more stable.
            Recommended values:
            * 1.0E-6 for minimal stability improving smoothing
            * 1.0E-3 a good value to start experiments; first results
              are visible
            * 1.0 for strong smoothing
  */

  double r_base = 20;
  double n_layers = 20;
  double lambda_n_s = 0.0;
  double search_r = 0.4;

  /*************************************************************************
  This function sets support radius parameter  of  hierarchical  (version 2)
  RBF constructor.

  Hierarchical RBF model achieves great speed-up  by removing from the model
  excessive (too dense) nodes. Say, if you have RBF radius equal to 1 meter,
  and two nodes are just 1 millimeter apart, you  may  remove  one  of  them
  without reducing model quality.

  Support radius parameter is used to justify which points need removal, and =
  alglib::rbfreport() which do not. If two points are less than
  SUPPORT_R*CUR_RADIUS  units  of distance apart, one of them is removed from
  the model. The larger  support radius  is, the faster model  construction  AND
  evaluation are.  However, too large values result in "bumpy" models.

  search_r       -   support radius coefficient, >=0.
              Recommended values are [0.1,0.4] range, with 0.1 being
              default value.

  *************************************************************************/

  /*
  Only for Gaussian or Bell
  S       -   RBF model, initialized by rbfcreate() call
  RBase   -   RBase parameter, RBase>0
  NLayers -   NLayers parameter, NLayers>0, recommended value  to  start
              with - about 5.
  LambdaNS-   >=0, nonlinearity penalty coefficient, negative values are
              not allowed. This parameter adds controllable smoothing to
              the problem, which may reduce noise. Specification of non-
              zero lambda means that in addition to fitting error solver
              will  also  minimize   LambdaNS*|S''(x)|^2  (appropriately
              generalized to multiple dimensions.

              Specification of exactly zero value means that no  penalty
              is added  (we  do  not  even  evaluate  matrix  of  second
              derivatives which is necessary for smoothing).

              Calculation of nonlinearity penalty is costly - it results
              in  several-fold  increase  of  model  construction  time.
              Evaluation time remains the same.

              Optimal  lambda  is  problem-dependent and requires  trial
              and  error.  Good  value to  start  from  is  1e-5...1e-6,
              which corresponds to slightly noticeable smoothing  of the
              function.  Value  1e-2  usually  means  that  quite  heavy
              smoothing is applied.

  TUNING ALGORITHM

  In order to use this algorithm you have to choose three parameters:
  * initial radius RBase
  * number of layers in the model NLayers
  * penalty coefficient LambdaNS

  Initial radius is easy to choose - you can pick any number  several  times
  larger  than  the  average  distance between points. Algorithm won't break
  down if you choose radius which is too large (model construction time will
  increase, but model will be built correctly).

  Choose such number of layers that RLast=RBase/2^(NLayers-1)  (radius  used
  by  the  last  layer)  will  be  smaller than the typical distance between
  points.  In  case  model  error  is  too large, you can increase number of
  layers.  Having  more  layers  will make model construction and evaluation
  proportionally slower, but it will allow you to have model which precisely
  fits your data. From the other side, if you want to  suppress  noise,  you
  can DECREASE number of layers to make your model less flexible (or specify
  non-zero LambdaNS).

  TYPICAL ERRORS

  1. Using too small number of layers - RBF models with large radius are not
     flexible enough to reproduce small variations in the  target  function.
     You  need  many  layers  with  different radii, from large to small, in
     order to have good model.

  2. Using  initial  radius  which  is  too  small.  You will get model with
     "holes" in the areas which are too far away from interpolation centers.
     However, algorithm will work correctly (and quickly) in this case.

  */

  double alpha = 10.0;
  /*
  for Multiquadric only
  f(r)=sqrt(r^2+Alpha^2) - rbf
  */

  bool v3tol = true;
  double tol = 0.0001;
  /*
  As of ALGLIB 3.20.0, version 3 models include biharmonic RBFs, thin  plate
  splines, multiquadrics.

  Version 3 models are fit  with  specialized  domain  decomposition  method
  which splits problem into smaller  chunks.  Models  with  size  less  than
  the DDM chunk size are computed nearly exactly in one step. Larger  models
  are built with an iterative linear solver. This function controls accuracy
  of the solver.

  desired precision:
          * must be non-negative
          * should be somewhere between 0.001 and 0.000001
          * values higher than 0.001 make little sense   -  you  may
            lose a lot of precision with no performance gains.
          * values below 1E-6 usually require too much time to converge,
            so they are silenly replaced by a 1E-6 cutoff value. Thus,
            zero can be used to denote 'maximum precision'.
  */

  RBFKind kind = RBFKind::TPS;

  ~RBFBasedWithNoTrain() {}

  template <typename T> void loadData(const T &x, const T &y) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::rbfcreate(1, 1, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    int N = x.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = x[i];
      data(i, 1) = y[i];
    }

    alglib::rbfsetpoints(model, data, N);
    alglib::rbfbuildmodel(model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  template <typename T> void loadData(const T &y) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::rbfcreate(1, 1, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    int N = y.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = i;
      data(i, 1) = y[i];
    }

    alglib::rbfsetpoints(model, data, N);
    alglib::rbfbuildmodel(model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  // U is 2d array
  // x and y are n_elems*n_dims_x and n_elems*n_dims_y 2d arrays
  template <typename U>
  void loadNDData(const U &x, const U &y, int n_dims_x, int n_dims_y,
                  int n_elems) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::real_2d_array data;
    int n_dims = n_dims_x + n_dims_y;
    data.setlength(n_elems, n_dims);

    for (int i = 0; i < n_elems; i++) {
      for (int j = 0; j < n_dims_x; j++) {
        data(i, j) = x[i][j];
      }
      for (int j = 0; j < n_dims_y; j++) {
        data(i, j + n_dims_x) = y[i][j];
      }
    }

    alglib::rbfcreate(n_dims_x, n_dims_y, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    alglib::rbfsetpoints(model, data, n_elems);
    alglib::rbfbuildmodel(model, report);

    this->n_dims_x = n_dims_x;
    this->n_dims_y = n_dims_y;
  }

  template <typename IdxT> double compute(IdxT idx) {
    alglib::real_1d_array idx_;
    idx_.setlength(1);
    idx_[0] = idx;
    alglib::real_1d_array value;
    value.setlength(1);
    alglib::rbfcalc(model, idx_, value);
    return value[0];
  }

  // idx and val are arrays of n_dims_x and n_dims_y sizes
  template <typename IdxT, typename ValueT>
  void compute(const IdxT &idx, ValueT &val) {
    alglib::real_1d_array idx_;
    idx_.setlength(n_dims_x);
    alglib::real_1d_array val_;
    val_.setlength(n_dims_y);
    for (int i = 0; i < n_dims_x; i++) {
      idx_[i] = idx[i];
    }
    alglib::rbfcalc(model, idx_, val_);
    for (int i = 0; i < n_dims_y; i++) {
      val[i] = val_[i];
    }
  }

  template <typename IdxT> double computeDerive(const IdxT &idx) {
    // todo
    alglib::real_1d_array idx_;
    idx_.setlength(1);
    idx_[0] = idx;
    alglib::real_1d_array value;
    value.setlength(1);
    alglib::real_1d_array derivative;
    derivative.setlength(1);
    alglib::rbfdiff(model, idx_, value, derivative);
    return derivative[0];
  }

  // idx and val are arrays of n_dims_x and n_dims_y sizes
  template <typename IdxT, typename DerivativeT>
  void computeDerive(const IdxT &idx, DerivativeT &der) {
    alglib::real_1d_array idx_;
    idx_.setlength(n_dims_x);
    alglib::real_1d_array val_, der_;
    val_.setlength(n_dims_y);
    der_.setlength(n_dims_y * n_dims_x);
    for (int i = 0; i < n_dims_x; i++) {
      idx_[i] = idx[i];
    }
    alglib::rbfdiff(model, idx_, val_, der_);
    for (int i = 0; i < n_dims_y * n_dims_x; i++) {
      der[i] = der_[i];
    }
  }

  void show1d(double step, double first, double last, bool plot_original = false, std::vector<double> * original_data = nullptr) {
    using SampleType = double;
    std::vector<SampleType> plotting_data = {};
    for (auto i = first; i < last; i += step) {
      plotting_data.push_back(compute<double>(i));
    }
    if (plot_original) {
      std::vector<SampleType> plotting_data_original = {};
      for (auto i = first; i < last; i += step) {
        plotting_data_original.push_back((*original_data)[static_cast<int>(i)]);
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
  int N = 102;
  RBFBasedWithNoTrain approximator;

  approximator.lambda_v = 1000.0; //smoothing parameter

  std::vector<double> data;

  for (int i = 0; i < N; i++) {
    data.push_back(std::rand());
  }

  approximator.loadData(data);

  
  approximator.show1d(1, 0, N, true, &data);
  approximator.show1d(0.1, 0, N, true, &data);

  N = 10;

  std::vector<std::vector<double>> data2d;
  for (int i = 0; i < N; i++) {
    data2d.push_back({});
    for (int j = 0; j < N; j++) {
      double r =
          sqrt((i - N / 2.0) * (i - N / 2.0) + (j - N / 2.0) * (j - N / 2.0));

      data2d[i].push_back(std::cos(r) / sqrt(r + 1));
      // IC(i, j, r, data2d[i][j]);
    }
  }

  std::vector<std::vector<double>> x_for_load;
  std::vector<std::vector<double>> y_for_load;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      x_for_load.push_back({});
      x_for_load[i * N + j].push_back(i);
      x_for_load[i * N + j].push_back(j);
      y_for_load.push_back({});
      y_for_load[i * N + j].push_back(data2d[i][j]);
    }
  }

  RBFBasedWithNoTrain approximator2d;
  approximator2d.kind = RBFKind::TPS;
  approximator2d.lambda_v = 10.0;
  approximator2d.loadNDData<decltype(x_for_load)>(x_for_load, y_for_load, 2, 1,
                                                  N * N);

  using namespace matplot;
  auto [X, Y] = meshgrid(iota(0, 0.1, N));
  auto Z = transform(X, Y, [&](double x, double y) {
    std::vector<double> idx = {x, y};
    std::vector<double> result = {0.0};
    approximator2d.compute<decltype(idx), decltype(result)>(idx, result);
    return result[0];
  });
  mesh(X, Y, Z);

  show();

  return 0;
}