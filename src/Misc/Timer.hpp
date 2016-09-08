/** Provides the Timer class.
 *
 * @file Timer.hpp
 *
 * @author Harrison Steggles
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <chrono>
#include <string>

class Timer {
	using TimePoint = typename std::chrono::time_point<std::chrono::system_clock>;
	using Duration = typename std::chrono::duration<double>;
    public:
		Timer();

		void start();
		void stop();
		void pause();
		void unpause();

		double getTicks();
		std::string formatTime(double time) const;

		bool isStarted();
		bool isPaused();

    private:
		TimePoint m_startTicks;
		TimePoint m_pausedTicks;
		Duration m_pausedDuration;

		bool m_isPaused;
		bool m_isStarted;
};



#endif // TIMER_HPP_
