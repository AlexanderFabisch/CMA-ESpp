/**
 * @file utils.h
 * Contains some utility functions.
 */

#pragma once

#include <algorithm>
#include <string>

template<typename T> T square(T d) { return d*d; }
template<typename T>
T maxElement(const T* rgd, int len) { return *std::max_element(rgd, rgd + len); }
template<typename T>
T minElement(const T* rgd, int len) { return *std::min_element(rgd, rgd + len); }
template<typename T>
int maxIndex(const T* rgd, int len) { return std::max_element(rgd, rgd + len) - rgd; }
template<typename T>
int minIndex(const T* rgd, int len) { return std::min_element(rgd, rgd + len) - rgd; }
/** sqrt(a^2 + b^2) numerically stable. */
double myhypot(double a, double b);
