#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <cmath>
#include "complex_num.h"
#include <sstream>

using namespace complex;

TEST_CASE("EQ_Test") {
  CartCompNum x1 (4, 66), x2 (4, 66);
  REQUIRE(x1 == x1);
  REQUIRE(x1 == x2);
  REQUIRE(x1 == CartCompNum(4, 66));
  x1 = CartCompNum(2, 5);
  REQUIRE(x1 == CartCompNum(2, 5));
}

TEST_CASE("Sums Test") {
  CartCompNum x1 (15, 32), x2(12, 8), x3(-12, -8);
}