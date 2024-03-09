add_library(Newton Newton.h Newton.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Newton PUBLIC Eigen3::Eigen)