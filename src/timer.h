#pragma once

#include <chrono>
#include <iostream>

class Timer
{
public:
  using Clock = std::chrono::high_resolution_clock;

  Timer(const std::string& _name, const bool _silent = false)
      : name(_name), silent(_silent), running(true), start_point(Clock::now()), duration(0.0)
  {}

  ~Timer()
  {
    if(!silent)
    {
      stop();
      std::cout << name << " took " << seconds() << "s.\n";
    }
  }

  void start()
  {
    running = true;
    start_point = Clock::now();
  }

  void stop()
  {
    if(running)
      duration += Clock::now() - start_point;

    running = false;
  }

  double seconds()
  {
    if(running)
      return (duration + Clock::now() - start_point).count();
    else
      return duration.count();
  }

private:
  const std::string name;
  const bool silent;
  bool running;
  typename Clock::time_point start_point;
  std::chrono::duration<double> duration;
};
