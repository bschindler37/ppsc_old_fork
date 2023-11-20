#ifndef _PPSC_TIMER_HPP
#define _PPSC_TIMER_HPP

// -----------------------------------------------------------------------
#include <chrono>
#include <iostream>
// -----------------------------------------------------------------------

#ifdef PPSC_OCA_TIMERS
  #define OCA_TIMER(s) ppsc::Timer tmr(s);
#else
  #define OCA_TIMER(s)
#endif

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

class Timer
{
public:
  Timer(std::string timer_name) : timer_name(timer_name), beg_(clock_::now()) {}
  ~Timer() {
    std::cout << "TIMER: " << elapsed() << " s"
	      << " -- " << timer_name << std::endl;
  }
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_>
      (clock_::now() - beg_).count(); }

private:
  std::string timer_name;
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;

};

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_TIMER_HPP
