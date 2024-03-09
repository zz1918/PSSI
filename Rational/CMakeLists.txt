add_subdirectory(Expression)
add_library(Rational rational.h rational.cpp)
target_link_libraries(Rational PUBLIC Expression)
target_include_directories(Rational PUBLIC Expression)