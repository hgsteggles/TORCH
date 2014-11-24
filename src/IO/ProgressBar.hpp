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
	ProgressBar(double tmax, int cpoint, const std::string msg, bool debug);
	bool update(double timeCurrent, double& dt_nextCheckpoint, bool output_on);
	void end(bool output_on);
	void reset(double tmax, int cpoint, const std::string msg);
private:
	typedef std::chrono::steady_clock Clock;
	typedef Clock::time_point time_point;
	typedef Clock::period period;
	typedef std::chrono::duration<float, period> duration;
	double timeTotal;
	int checkpoint;
	std::string messageProgress;
	double percentProgress;
	int checkpointProgress;
	time_point clockStart, clockLast;
	duration speedEMA;
	double smoothing;
	bool debugging;
};

#endif //PROGRESSBAR_HPP_
