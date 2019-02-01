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
TEST_CASE("EqualTest") {
  CartCompNum x1(2,33), x2(2, 33);
  REQUIRE(x1 == x1);
  REQUIRE(x1 == x2);
  REQUIRE(x1 == CartCompNum(2, 33));
  x1 = CartCompNum(3, 10);
  REQUIRE(x1 == CartCompNum(3, 10));
}

TEST_CASE("SumTest") {
  CartCompNum x1(20, 32), x2(12, 8), x3(-12, -8);
  REQUIRE((x1 + x2) == CartCompNum(32, 40));
  REQUIRE((x1 + x3) == CartCompNum(8, 24));
  x1 += x1;
  REQUIRE(x1 == CartCompNum(40, 64));
}

TEST_CASE("MinusTest") {
  CartCompNum x1(20, 32), x2(12, 8), x3(-12, -8);
  REQUIRE((x1 - x2) == CartCompNum(8, 24));
  REQUIRE((x2 - x1) == CartCompNum(-8, -24));
  REQUIRE((x2 - x3) == CartCompNum(24, 16));
}

TEST_CASE("MultTest") {
  CartCompNum x1(7, 3), x2(5, -8), x3(-10, -4), x4(8, 0);
  REQUIRE(x1 * x2 == CartCompNum (59, -41));
  REQUIRE(x3 * x4 == CartCompNum (-80, -32));
}

TEST_CASE("DivTest") {
  CartCompNum x1(2, 5), x2(3, -2);
  REQUIRE((x1 / x2) == CartCompNum(-4.0 / 13.0, 19.0 / 13));
}

TEST_CASE("GET_SET_ARG_ABS") {
  CartCompNum x1(2, 5), x2(-3, -8), x3(0, 8), x4(0, -8);
  x1.SetIm(9);
  REQUIRE(x1 == CartCompNum(2, 9));
  x1.SetRe(-8);
  REQUIRE(x1 == CartCompNum(-8, 9));
  REQUIRE(x1.GetIm() == 9);
  REQUIRE(x1.GetRe() == -8);
  REQUIRE(std::abs(x1.arg() - 2.29) < 0.008);
  REQUIRE(std::abs(x1.abs() - 12.04159) < 0.002);
  x1.SetIm(-5);
  REQUIRE(std::abs(x1.arg()) > 2.58);
  REQUIRE(x3.arg() == (M_PI / 2));
  REQUIRE(x4.arg() == -(M_PI / 2));
  REQUIRE(x2 == PolarToCart(CartToPolar(x2)));
}

//polar numbers tests

TEST_CASE("STRING_IN"){
  PolarCompNum x2(12, 12);
  CartCompNum x1(10, 10);
  std::stringstream ss;
  ss << x2;
  REQUIRE(ss.str() == "(r=12,f=12)");
  ss.str("");
  ss << x1;
  REQUIRE(ss.str() == "(10,i*10)");
}

TEST_CASE("STRINGTEST_OUT"){
  PolarCompNum x2;
  std::stringstream ss;
  std::stringstream ss1;
  ss << "1 2";
  CartCompNum x1;
  ss >> x1;
  REQUIRE(x1 == CartCompNum(1,2));
  ss1 << ("2 3");
  ss1 >> x2;
  REQUIRE(x2 == PolarCompNum(2, 3));
}

TEST_CASE("EqualsTestPolar") {
  PolarCompNum x1(20, 32), x2(20, 32);
  REQUIRE(x1 == x1);
  REQUIRE(x1 == x2);
  REQUIRE(x1 == PolarCompNum(20, 32));
}

TEST_CASE("SumTestPolar") {
  PolarCompNum x1(20, 32), x2(12, 8), x3(-12, -8), x4(10, 10), x5(12, 8);
  REQUIRE((x1 + x2) == CartToPolar(PolarToCart(x1) + PolarToCart(x2)));
  REQUIRE((x1 + x3) == CartToPolar(PolarToCart(x1) + PolarToCart(x3)));
  x5 += x4;
  REQUIRE((x4 + x2) == x5);
}

TEST_CASE("MinusTestPolar") {
  PolarCompNum x1(20, 32), x2(12, 8), x3(-12, 0), x4(-12, M_PI), x5(10, 0);
  REQUIRE((x1 - x2) == CartToPolar(PolarToCart(x1) - PolarToCart(x2)));
  REQUIRE((x1 - x3) == CartToPolar(PolarToCart(x1) - PolarToCart(x3)));
  x1 -= x5;
  REQUIRE(x1 == CartToPolar(PolarToCart(PolarCompNum(20, 32)) - PolarToCart(PolarCompNum(10, 0))));
  REQUIRE((x1 - PolarCompNum(10)) == CartToPolar(PolarToCart(x1) - PolarToCart(x5)));
  x1 = PolarCompNum(10, 10);
  REQUIRE(x1 == PolarCompNum(10, 10));
}


TEST_CASE("MultTestPolar") {
  PolarCompNum x1(7, 3), x2(3, -8), x3(-10, -4), x4(8, 0);
  REQUIRE((x1 * x2) == PolarCompNum(21, -5));
  REQUIRE((x3 * PolarCompNum(8)) == (x3 * x4));
}

TEST_CASE("DivTestPolar") {
  PolarCompNum x1(2, 5), x2(3, -2);
  REQUIRE((x1 / x2) == PolarCompNum(2.0 / 3, 7));
}

TEST_CASE("GET_SET_POLAR") {
  PolarCompNum x1(9, 9);
  x1.SetF(8);
  REQUIRE(x1 == PolarCompNum(9, 8));
  x1.SetR(8);
  REQUIRE(x1 == PolarCompNum(8, 8));
}

TEST_CASE("OTHERTESTS") {
  CartCompNum x1, x2(10), x3(CartCompNum(10, 19));
  REQUIRE(x1 == CartCompNum());
  REQUIRE(x2 == CartCompNum(10));
  REQUIRE(x3 == CartCompNum(10, 19));
  CartCompNum x4 = x3;
  REQUIRE(x4 == CartCompNum(10, 19));
  x4 -= x3;
  REQUIRE(x4 == CartCompNum(0, 0));
  PolarCompNum y1, y2(10);
  REQUIRE(y2 == PolarCompNum(y2));
  REQUIRE(y1 == PolarCompNum());
  REQUIRE(y2.abs() == 10);
  REQUIRE(y2.arg() == 0);
}