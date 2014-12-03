#include "Timer.hpp"

#include <sstream>
#include <chrono>

Timer::Timer()
: m_startTicks(std::chrono::system_clock::now())
, m_pausedTicks(std::chrono::system_clock::now())
, m_pausedDuration(Duration::zero())
, m_isPaused(false)
, m_isStarted(false)
{

}

void Timer::start() {
	m_isStarted = true;
	m_isPaused = false;

	m_startTicks = std::chrono::system_clock::now();
	m_pausedTicks = m_startTicks;
	m_pausedDuration = Duration::zero();
}

void Timer::stop() {
	m_isStarted = false;
	m_isPaused = false;
}

void Timer::pause() {
	//If the timer is running and isn't already paused
	if ( m_isStarted && !m_isPaused ) {
		m_isPaused = true;
		m_pausedTicks = std::chrono::system_clock::now();
	}
}

void Timer::unpause() {
	//If the timer is running and paused
	if ( m_isStarted && m_isPaused ) {
		m_isPaused = false;
		m_pausedDuration += std::chrono::system_clock::now() - m_pausedTicks;
	}
}

double Timer::getTicks() {
	Duration dur = Duration::zero();

	if ( m_isStarted ) {
		if ( m_isPaused )
			dur = m_pausedTicks - m_startTicks - m_pausedDuration;
		else
			dur = std::chrono::system_clock::now() - m_startTicks - m_pausedDuration;
	}

	return dur.count();
}

bool Timer::isStarted() {
	return m_isStarted;
}

bool Timer::isPaused() {
	return m_isPaused && m_isStarted;
}

std::string Timer::formatTime(double time) const {
	int hours = (int)(time/3600.0);
	int minutes = (int)(time/60.0 - 60*hours);
	int seconds = (int)(time - 60*minutes - 3600*hours + 0.5);

	std::stringstream ss;
	ss << hours << "h:" << minutes << "m:" << seconds << "s.";
	return ss.str();
}


