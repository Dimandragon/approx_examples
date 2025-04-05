add_rules("mode.debug", "mode.release")

set_languages("c++23")

--xmake f --cxx=g++ --cc=gcc -m release --debugger=lldb
--xmake project -k compile_commands


target("icecream")
    set_kind("headeronly")
    --add_includedirs("$(projectdir)/icecream-cpp", {public = true})
    add_includedirs("$(projectdir)/my_icecream", {public = true})

target("boost_interpolate")
    set_kind("headeronly")
    add_includedirs("interpolate_from_boost", {public = true})

target("pocketfft")
    set_kind("headeronly")
    --add_headerfiles("pocketfft/pocketfft_hdronly.h", {public = true})
    add_includedirs("pocketfft", {public = true})

target("alglib")
    set_kind("static")
    add_includedirs("alglib-cpp/src", {public = true})
    add_files("alglib-cpp/src/*.cpp")


target("matplot++_external")
    set_kind("headeronly")
    --on_load(function(target)
        --#include("core.base.task")
        --task.run("build matplot++")
    --end)
    --add_headerfiles("$(projectdir)/matplotplusplus/source/matplot/matplot.h", {public = true})
    --add_headerfiles("$(projectdir)/matplotplusplus/source/3rd_party/nodesoup/include/nodesoup.h", {public = true})
    add_includedirs("matplotplusplus/build/local/source/matplot", {public = true})
    add_includedirs("matplotplusplus/source", {public = true})
    add_includedirs("matplotplusplus/source/3rd_party/nodesoup", {public = true})

    --add_headerfiles("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/include/glad/glad.h", {public = true})
    --add_headerfiles("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/include/KHR/khrplatform.h", {public = true})

    add_links("$(projectdir)/matplotplusplus/build/local/source/matplot/libmatplot.a", {public = true})
    --add_links("$(projectdir)/matplotplusplus/build/local/source/matplot/libmatplot_opengl.a", {public = true})
    add_links("$(projectdir)/matplotplusplus/build/local/source/3rd_party/libnodesoup.a", {public = true})
    --add_links("$(projectdir)/matplotplusplus/build/local/_deps/glad-build/libglad.a", {public = true})


target("fft-example")
    set_kind("binary")
    add_files("examples/fft.cpp")
    add_deps("pocketfft")
    add_deps("matplot++_external")
    
target("fs_based_approx_example")
    set_kind("binary")
    add_files("examples/fs_based_approx.cpp")
    add_deps("matplot++_external")
    add_deps("pocketfft")
    add_deps("icecream")

target("makima_approximator_example")
    set_kind("binary")
    add_files("examples/makima_approx.cpp")
    add_deps("matplot++_external")
    add_deps("boost_interpolate")
    add_deps("icecream")

target("rbf_spline_approx")
    set_kind("binary")
    add_files("examples/rbf_approx.cpp")
    add_deps("matplot++_external")
    add_deps("alglib")
    add_deps("icecream")
