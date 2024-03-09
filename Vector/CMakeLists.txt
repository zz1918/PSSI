add_library(Vector Vector.h Vector.cpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Vector PUBLIC Eigen3::Eigen)