#include "ProgressBar.hpp"

#include <stdio.h>
#include <iostream>
#include <iomanip>


ProgressBar::ProgressBar(double tmax, int cpoint, const std::string msg, bool debug)
	: timeTotal(tmax)
	, checkpoint(cpoint)
	, messageProgress(msg)
	, percentProgress(0)
	, checkpointProgress(0)
	, clockStart(std::chrono::steady_clock::now())
	, clockLast(clockStart)
	, smoothing(0.5)
	, debugging(debug)
{
	speedEMA = Clock::now()-Clock::now();
}

bool ProgressBar::update(double time_current, double& dt_nextCheckpoint, bool output_on) {
	static bool once = false;
	double percent_current = 100*time_current/timeTotal;
	if (output_on && percent_current - percentProgress >= 0.1) {
		int barWidth = 20;
		int prcnow = (int)(percent_current + 0.5);
		int pos = barWidth * prcnow / 100.0;
		std::cout << messageProgress << ": [";
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] ";

		if (percent_current != 100) {
			duration speed = (Clock::now() - clockLast)/(percent_current-percentProgress);
			if (!once) {
				speedEMA = speed;
				once = true;
			}
			else
				speedEMA = smoothing*speed + (1.0 - smoothing)*speedEMA;
			duration time_left = speedEMA*(100 - percent_current);
			std::chrono::minutes minutes_left = std::chrono::duration_cast<std::chrono::minutes>(time_left);
			std::chrono::seconds seconds_left = std::chrono::duration_cast<std::chrono::seconds>(time_left - minutes_left);
			std::cout << (int)(percent_current*10 + 0.5)/10.0 << " % (";
			std::cout << minutes_left.count() << "m " << seconds_left.count() << "s)                     \r";
		}
		else {
			duration total_time = Clock::now() - clockStart;
			std::chrono::minutes total_minutes = std::chrono::duration_cast<std::chrono::minutes>(total_time);
			std::chrono::seconds total_seconds = std::chrono::duration_cast<std::chrono::seconds>(total_time - total_minutes);
			std::cout << "(time taken = " << total_minutes.count() << "m " << total_seconds.count() << "s)               \r";
		}
		std::cout << std::flush;

		percentProgress = (int)(percent_current*10 + 0.5)/10.0;
		clockLast = Clock::now();

		if (percentProgress == 100 || debugging)
			std::cout << std::endl;
	}

	bool isTimeToUpdate = false;
	if (percent_current - checkpointProgress > checkpoint - 0.000001) {
		checkpointProgress = (int)(percent_current + 0.5);
		isTimeToUpdate = true;
	}

	dt_nextCheckpoint = (checkpointProgress + checkpoint - percent_current)*timeTotal/100.0;

	return isTimeToUpdate;
}
/*
bool ProgressBar::update(double time_current, bool output_on) {
	using namespace std::chrono;
	double percent_current = 100*time_current/timeTotal;
	if (output_on && percent_current - percentProgress >= 1) {
		int barWidth = 20;
		int prcnow = (int)(percent_current + 0.5);
		int pos = barWidth * prcnow / 100.0;
		std::cout << messageProgress << ": [";
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		duration time_taken = Clock::now() - clockStart;
		duration time_left = time_taken * (100/percent_current - 1);
		minutes minutes_left = duration_cast<minutes>(time_left);
		seconds seconds_left = duration_cast<seconds>(time_left - minutes_left);

		std::cout << "] " << prcnow << " % (";
		std::cout << minutes_left.count() << "m " << seconds_left.count() << "s)                     \r";
		std::cout << std::flush;

		percentProgress = prcnow;
		if (percentProgress == 100)
			std::cout << std::endl;
	}
	if (percent_current - checkpointProgress >= checkpoint) {
		checkpointProgress = (int)(percent_current + 0.5);
		return true;
	}
	return false;
}
*/
void ProgressBar::end(bool output_on) {
	if (percentProgress != 100) {
		percentProgress = 100;
		if (output_on) {
			int barWidth = 20;
			int pos = barWidth;
			std::cout << messageProgress << ": [";
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}

			std::cout << "] " << 100 << " % (";
			std::cout << 0 << "m " << 0 << "s)                     \r";
			std::cout << std::endl;
		}
	}
}

void ProgressBar::reset(double tmax, int cpoint, const std::string msg) {
	timeTotal = tmax;
	checkpoint = cpoint;
	messageProgress = msg;
	percentProgress = 0;
	checkpointProgress = 0;
	clockStart = std::chrono::steady_clock::now();
}
