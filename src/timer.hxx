#ifndef TIMER_HXX
#define TIMER_HXX

#include <time.h>
#include <string>
#include <sstream>

struct Timer
{
  clockid_t clock_id = CLOCK_THREAD_CPUTIME_ID;
  bool running = false;
  double time_start = 0.0;
  double time_elapsed = 0.0;
  int segment_n = 0;

  void start() {
    if (running) return;
    running = true;
    time_start = cur_time();
  }
  void stop() {
    if (!running) return;
    double time_stop = cur_time();
    running = false;
    time_elapsed += time_stop - time_start;
    ++segment_n;
  }
  double elapsed() const {
    return time_elapsed + (running ? cur_time() - time_start : 0);
  }
  std::string elapsed_str() const {
    std::stringstream ss;
    ss << elapsed() << " seconds (" << segment_n << " segments)";
    return ss.str();
  }

private:
  double cur_time() const {
    struct timespec ts;
    clock_gettime(clock_id, &ts);
    return ts.tv_sec + (double) ts.tv_nsec / 1e9;
  }
};

#endif // TIMER_HXX
