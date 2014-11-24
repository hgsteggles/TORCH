/*
 * Timer.h
 *
 *  Created on: 19 Oct 2014
 *      Author: harry
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <cstdint>
#include <string>
#include <chrono>

class Timer {
	using TimePoint = typename std::chrono::time_point<std::chrono::system_clock>;
	using Duration = typename std::chrono::duration<double>;
    public:
		//Initializes variables
		Timer();

		//The various clock actions
		void start();
		void stop();
		void pause();
		void unpause();

		//Gets the timer's time
		double getTicks();
		std::string formatTime(double time) const;

		//Checks the status of the timer
		bool isStarted();
		bool isPaused();

    private:
		//The clock time when the timer started
		TimePoint m_startTicks;

		//The ticks stored when the timer was paused
		TimePoint m_pausedTicks;
		Duration m_pausedDuration;

		//The timer status
		bool m_isPaused;
		bool m_isStarted;
};



#endif /* TIMER_H_ */
