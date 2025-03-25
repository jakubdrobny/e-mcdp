#include "Timer.hpp"
#include <chrono>

Timer::Timer() : start_(std::chrono::high_resolution_clock::now()) {}

void Timer::reset() { start_ = std::chrono::high_resolution_clock::now(); }

template <typename Duration>
  requires std::derived_from<Duration, std::chrono::duration<typename Duration::rep, typename Duration::period>>
double Timer::elapsed() const {
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<Duration>(end - start_).count();
}

template double Timer::elapsed<std::chrono::nanoseconds>() const;
template double Timer::elapsed<std::chrono::microseconds>() const;
template double Timer::elapsed<std::chrono::milliseconds>() const;
template double Timer::elapsed<std::chrono::seconds>() const;
template double Timer::elapsed<std::chrono::minutes>() const;
template double Timer::elapsed<std::chrono::hours>() const;
