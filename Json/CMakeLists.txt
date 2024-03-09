add_library(Json json.h json.cpp json.hpp)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Json PUBLIC Eigen3::Eigen)