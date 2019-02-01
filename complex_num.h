#pragma once
#include <iosfwd>
namespace complex {
class CartCompNum {
 private:
  double re_, im_;

 public:
  CartCompNum();
  //explicit CartCompNum(double);
  CartCompNum(double, double);
  CartCompNum(const CartCompNum &);

  double GetRe();
  double GetIm();

  void SetRe(double);
  void SetIm(double);

  CartCompNum operator=(const CartCompNum &);

  CartCompNum operator+();
  CartCompNum operator+(const CartCompNum &);
  CartCompNum operator+(const double &);
  CartCompNum operator+=(const CartCompNum &);
  friend CartCompNum operator+(const double &,
                               const CartCompNum &);

  CartCompNum operator-();
  CartCompNum operator-(const CartCompNum &);
  CartCompNum operator-(const double &);
  CartCompNum operator-=(const CartCompNum &);
  friend CartCompNum operator-(const double &,
                               const CartCompNum &);

  CartCompNum operator*(const CartCompNum &);
  CartCompNum operator*(const double &);
  CartCompNum operator*=(const CartCompNum &);
  friend CartCompNum operator*(const double &,
                               const CartCompNum &);

  CartCompNum operator/(const CartCompNum &);
  CartCompNum operator/(const double &);
  CartCompNum operator/=(const CartCompNum &);

  double abs();
  double arg();

  friend std::ostream &operator<<(std::ostream &, CartCompNum &);
  friend std::istream &operator>>(std::istream &, CartCompNum &);
  friend bool operator==(const CartCompNum &, const CartCompNum &);
};

class PolarCompNum {
 private:
  double r_, f_;
 public:
  PolarCompNum();
  explicit PolarCompNum(double r);
  PolarCompNum(double, double);
  PolarCompNum(const PolarCompNum &);

  double GetR();
  double GetF();
  void SetR(double);
  void SetF(double);

  PolarCompNum operator=(const PolarCompNum &);

  PolarCompNum operator+();
  PolarCompNum operator+(const PolarCompNum &);
  PolarCompNum operator+=(const PolarCompNum &);

  PolarCompNum operator-();
  PolarCompNum operator-(const PolarCompNum &);
  PolarCompNum operator-=(const PolarCompNum &);

  PolarCompNum operator*(const PolarCompNum &);
  double arg() const;
  double abs() const;

  PolarCompNum operator/(const PolarCompNum &);

  friend std::ostream &operator<<(std::ostream &, PolarCompNum &);
  friend std::istream &operator>>(std::istream &, PolarCompNum &);

  friend bool operator==(const PolarCompNum &, const PolarCompNum &);
  friend PolarCompNum bpow(const PolarCompNum& value, int n);
  friend CartCompNum PolarToCart(const PolarCompNum &);

};
CartCompNum PolarToCart(const PolarCompNum &);
PolarCompNum CartToPolar(CartCompNum &);
CartCompNum bpow(const CartCompNum &value, int n);
}
