#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
public:
  Timer();

  void reset();

  template <typename Duration>
    requires std::derived_from<Duration,
                               std::chrono::duration<typename Duration::rep,
                                                     typename Duration::period>>
  double elapsed() const;

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

#endif
