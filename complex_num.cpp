#include "complex_num.h"
#include "cmath"
#include "iostream"

using namespace complex;

CartCompNum::CartCompNum() {};

CartCompNum::CartCompNum(double value) {
  re_ = value;
  im_ = 0;
}

CartCompNum::CartCompNum(double re, double im) {
  re_ = re;
  im_ = im;
}

CartCompNum::CartCompNum(const CartCompNum &prev) = default;

double CartCompNum::GetRe() {
  return re_;
}
double CartCompNum::GetIm() {
  return im_;
}

void CartCompNum::SetRe(double value) {
  re_ = value;
}

void CartCompNum::SetIm(double value) {
  im_ = value;
}

CartCompNum &CartCompNum::operator=(const CartCompNum &value) {
  re_ = value.re_;
  im_ = value.im_;
  return *this;
}

CartCompNum CartCompNum::operator+() {
  return CartCompNum(re_, im_);
}

CartCompNum CartCompNum::operator+(const double &value) {
  return CartCompNum(re_ + value, im_);
}

CartCompNum CartCompNum::operator+=(const CartCompNum &value) {
  re_ += value.re_;
  im_ += value.im_;
  return CartCompNum(re_, im_);
}

CartCompNum CartCompNum::operator+(const CartCompNum &value) {
  return CartCompNum(re_ + value.re_, im_ + value.im_);
}

CartCompNum complex::operator+(const double &first,
                               const CartCompNum &second) {
  return CartCompNum(first + second.re_, second.im_);
}

CartCompNum CartCompNum::operator-() {
  return CartCompNum(re_ * (-1), im_ * (-1));
}

CartCompNum CartCompNum::operator-(const double &value) {
  return CartCompNum(re_ - value, im_);
}

CartCompNum CartCompNum::operator-(const CartCompNum &value) {
  return CartCompNum(re_ - value.re_, im_ - value.im_);
}

CartCompNum CartCompNum::operator-=(const CartCompNum &value) {
  re_ -= value.re_;
  im_ -= value.im_;
  return CartCompNum(re_, im_);
}

CartCompNum complex::operator-(const double &first,
                               const CartCompNum &second) {
  return CartCompNum(first - second.re_, second.im_);
}

CartCompNum CartCompNum::operator*(const CartCompNum &value) {
  return CartCompNum(re_ * value.re_ - im_ * value.im_, re_ * value.im_ + im_ * value.re_);
}
CartCompNum CartCompNum::operator*(const double &value) {
  return CartCompNum(value * re_, value * im_);
}
CartCompNum CartCompNum::operator*=(const CartCompNum &value) {
  auto re_buff = re_;
  re_ = re_ * value.re_ - im_ * value.im_;
  im_ = re_buff * value.im_ + im_ * value.re_;
  return CartCompNum(re_, im_);
}
CartCompNum complex::operator*(const double &first,
                               const CartCompNum &second) {
  return CartCompNum(first * second.re_, first * second.im_);
}

CartCompNum CartCompNum::operator/(const double &value) {
  return CartCompNum(re_ / value, im_ / value);
}

CartCompNum CartCompNum::operator/(const CartCompNum &value) {
  CartCompNum buffer;
  double mult = value.re_ * value.re_ + value.im_ * value.im_;
  buffer.re_ = (re_ * value.re_ + im_ * value.im_) / mult;
  buffer.im_ = (im_ * value.re_ - re_ * value.im_) / mult;
  return buffer;
}

CartCompNum CartCompNum::operator/=(const CartCompNum &value) {
  auto buff = *this;
  re_ = (buff / value).re_;
  im_ = (buff / value).im_;
  return *this;
}

double CartCompNum::abs() const {
  return (sqrt(re_ * re_ + im_ * im_));
}

double CartCompNum::arg() const {
  return atan2(im_, re_);
}

CartCompNum complex::bpow(const CartCompNum &value, int n) {
  CartCompNum count(1, 0);
  CartCompNum x(1.0);
  CartCompNum buff_val = value;
  int buff_degree_ = n;
  if (!n) return CartCompNum(1);
  while (n) {
    if (n % 2 == 0) {
      n /= 2;
      buff_val = buff_val * buff_val;
    } else {
      n--;
      count = count * buff_val;
    }
  }
  if (buff_degree_ > 0) {
    return count;
  } else {
    return x / count;
  }
}

std::ostream &complex::operator<<(std::ostream &out, CartCompNum &value) {
  out << "(" << value.re_ << ",i*" << value.im_ << ")";
  return out;
}
std::istream &complex::operator>>(std::istream &in, CartCompNum &value) {
  in >> value.re_ >> value.im_;
  return in;
}

bool complex::operator==(const CartCompNum &left,
                         const CartCompNum &right) {
  return (std::abs(left.re_ - right.re_) < 0.00001 &&
      std::abs(left.im_ - right.im_) < 0.00001);
}
//Polar Complex Numbers
PolarCompNum::PolarCompNum() : r_(0), f_(0) {};
PolarCompNum::PolarCompNum(double value) : r_(value), f_(0) {}
PolarCompNum::PolarCompNum(double r, double f = 0) : r_(std::abs(r)), f_(f) {}
PolarCompNum::PolarCompNum(const PolarCompNum &value) : r_(value.r_), f_(value.f_) {};

double PolarCompNum::GetR() {
  return r_;
}

double PolarCompNum::GetF() {
  return f_;
}

void PolarCompNum::SetR(double value) {
  r_ = value;
}

void PolarCompNum::SetF(double value) {
  f_ = value;
}

PolarCompNum PolarCompNum::operator=(const PolarCompNum &value) {
  r_ = value.r_;
  f_ = value.f_;
  return PolarCompNum(r_, f_);
}

PolarCompNum PolarCompNum::operator+() {
  return PolarCompNum(r_, f_);
}

PolarCompNum PolarCompNum::operator+(const PolarCompNum &value) {
  CartCompNum buffer = PolarToCart(value) + PolarToCart(*this);
  return CartToPolar(buffer);
}

PolarCompNum PolarCompNum::operator+=(const PolarCompNum &value) {
  CartCompNum buff = PolarToCart(*this) + PolarToCart(value);
  r_ = CartToPolar(buff).r_;
  f_ = CartToPolar(buff).f_;
  return *this;
}

PolarCompNum PolarCompNum::operator-() {
  return PolarCompNum(-r_, f_);
}

PolarCompNum PolarCompNum::operator-(const PolarCompNum &value) {
  CartCompNum buffer = PolarToCart(*this) - PolarToCart(value);
  return CartToPolar(buffer);
}

PolarCompNum PolarCompNum::operator-=(const PolarCompNum &value) {
  CartCompNum buff = PolarToCart(*this) - PolarToCart(value);
  r_ = CartToPolar(buff).r_;
  f_ = CartToPolar(buff).f_;
  return *this;
}

PolarCompNum PolarCompNum::operator*(const PolarCompNum &value) {
  double r_buff_ = r_ * value.r_;
  double f_buff_ = f_ + value.f_;
  return PolarCompNum(r_buff_, f_buff_);
}

PolarCompNum PolarCompNum::operator*(double value) {
  return PolarCompNum(r_ * value, f_);
}

PolarCompNum PolarCompNum::operator*=(const PolarCompNum &value) {
  double r_buff_ = r_ * value.r_;
  double f_buff_ = f_ + value.f_;
  r_ = r_buff_;
  f_ = f_buff_;
  return *this;
}

PolarCompNum PolarCompNum::operator*=(double value) {
  r_ *= value;
  return *this;
}

PolarCompNum PolarCompNum::operator/(const PolarCompNum &value) {
  double r_buff_ = r_ / value.r_;
  double f_buff_ = f_ - value.f_;
  return PolarCompNum(r_buff_, f_buff_);
}

PolarCompNum PolarCompNum::operator/(double value){
  return PolarCompNum(r_ / value, f_);
}

PolarCompNum PolarCompNum::operator/=(double value){
  r_ /= value;
  return *this;
}

PolarCompNum PolarCompNum::operator/=(const PolarCompNum &value){
  double r_buff_ = r_ / value.r_;
  double f_buff_ = f_ - value.f_;
  r_ = r_buff_;
  f_ = f_buff_;
  return *this;
}

PolarCompNum complex::bpow(const PolarCompNum &value, int n) {
  return PolarCompNum(pow(value.r_, n), value.f_ * n);
}

double PolarCompNum::arg() const {
  return atan2(this->f_, this->r_) * 180 / M_PI;
}

double PolarCompNum::abs() const {
  return sqrt(this->r_ * this->r_ + this->f_ * this->f_);
}

std::ostream &complex::operator<<(std::ostream &out, PolarCompNum &c) {
  out << "(r=" << c.r_ << ",f=" << c.f_ << ")";
  return out;
}

std::istream &complex::operator>>(std::istream &in, PolarCompNum &c) {
  in >> c.r_ >> c.f_;
  return in;
}

bool complex::operator==(const PolarCompNum &first, const PolarCompNum &second) {
  return (first.r_ == second.r_ && first.f_ == second.f_);
}

CartCompNum complex::PolarToCart(const PolarCompNum &value) {
  return CartCompNum(value.r_ * cos(value.f_), value.r_ * sin(value.f_));
}

PolarCompNum complex::CartToPolar(const CartCompNum &value) {
  return PolarCompNum(value.abs(), value.arg());
}


