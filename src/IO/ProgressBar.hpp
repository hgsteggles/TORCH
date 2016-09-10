/** Provides the ProgressBar class.
 *
 * @file ProgressBar.hpp
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014 - the first version.
 */

#ifndef PROGRESSBAR_HPP_
#define PROGRESSBAR_HPP_

#include <stdint.h>
#include <chrono>
#include <string>

/**
 * @class ProgressBar
 *
 * @brief UI feature: an indicator of progress is displayed in the terminal.
 *
 * @version 0.8, 24/11/2014
 */
class ProgressBar
{
public:
	ProgressBar(double maxVal, int logPeriodMillis);

	bool timeToUpdate();
	void update(double currentVal);
	void end();

	std::string getFullString();
	std::string getFinalString();
	std::string getBarString();
	std::string getPercentDoneString();
	std::string getTimeLeftString();
	std::string getTimeTakenString();

private:
	typedef std::chrono::steady_clock Clock;
	typedef Clock::time_point time_point;
	typedef Clock::period period;
	typedef std::chrono::duration<float, period> duration;

	std::string durationToString(duration time_left);

	bool isStarting;
	double maxVal;
	double currentPercent;
	int logPeriod;

	int barWidth;

	time_point clockStart;	
	time_point clockLast;
	duration speedEMA;
	duration timeLeft;
	double smoothing;
};

#endif //PROGRESSBAR_HPP_
