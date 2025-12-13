#ifndef microsecond_execution_timer_HPP
#define microsecond_execution_timer_HPP

#include <functional>
#include <chrono> 

// Create a helper function that will return a time of execution in milliseconds (general form copied from stack-exchange)
double microsecond_execution_timer(const std::function<void(void)>& fcn)
{
  //get ready to time our function
  using std::chrono::high_resolution_clock; 
  using std::chrono::duration; 

  auto t0 = high_resolution_clock::now(); 
  fcn(); 
  auto t1 = high_resolution_clock::now(); 

  duration<double, std::micro> elapsed_time = t1 - t0; 

  return elapsed_time.count(); 
} 

#endif 