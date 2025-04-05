# about
This is repo with some numerical approximations examples
# dependencies
 - git
 - c++ compiler (like gcc or clang)
 - xmake build system
 - cmake build system
 - gnuplot
# get started 
 1) clone repo via (you must have github ssh key)
```
git clone git@github.com:Dimandragon/approx_examples.git
cd approx_examples
```
 2) clone submodules via
 ```
 git submodule update --init --recursive
 ```
 3) build matplotplusplus
 ```
 ./matplot_build.sh
 ```
 4) build projects and create ide files
 ```
 xmake
 xmake project -k compile_commands
 ```
 5) run example via
 ```
 xmake r example_name
 ```
 # example names
 - fft-example
 - fs_based_approx_example
 - makima_approximator_example
 - rbf_spline_approx