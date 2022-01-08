#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <vector>

template <typename type>
class timer
{
private:
    std::vector<std::chrono::time_point<std::chrono::_V2::system_clock,std::chrono::duration<long long,std::ratio<1,1000000000>>>> time;
    int counter = 0;
    static int timerCounter;
public:
    timer()
    {
        time.push_back(std::chrono::high_resolution_clock::now());
        timerCounter++;
    }
    ~timer(){timerCounter--;}
    std::chrono::duration<long long,std::ratio<1,1000>>::rep getTime()
    {
        time.push_back(std::chrono::high_resolution_clock::now());
        counter++;
        return std::chrono::duration_cast<type>(time.at(counter) - time.at(0)).count();
    }
    std::chrono::duration<long long,std::ratio<1,1000>>::rep getPointTime()
    {
        time.push_back(std::chrono::high_resolution_clock::now());
        counter++;
        return std::chrono::duration_cast<type>(time.at(counter) - time.at(counter-1)).count();
    }
    static int getTimerCounter()
    {
        return timerCounter;
    }
};

template <typename type>
int timer<type>::timerCounter = 0;

#endif // TIMER_H
