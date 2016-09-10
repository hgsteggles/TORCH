#include "ProgressBar.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>

ProgressBar::ProgressBar(double maxVal, int logPeriodMillis)
	: isStarting(true)
	, maxVal(maxVal)
	, currentPercent(0)
	, logPeriod(logPeriodMillis)
	, barWidth(10)
	, clockStart(Clock::now())
	, clockLast(Clock::now())
	, speedEMA(Clock::now() - Clock::now())
	, timeLeft(Clock::now() - Clock::now())
	, smoothing(0.5)
{

}

bool ProgressBar::timeToUpdate() {
	duration d = Clock::now() - clockLast;

	std::chrono::milliseconds millis = std::chrono::duration_cast<std::chrono::milliseconds>(d);

	return millis.count() > logPeriod;
}

std::string ProgressBar::getFullString() {
	return getBarString() + " " + getPercentDoneString() + " " + getTimeLeftString();
}

std::string ProgressBar::getFinalString() {
	return getBarString() + " " + getPercentDoneString() + " " + getTimeTakenString();
}

std::string ProgressBar::getBarString() {
	int pos = (int)(barWidth * currentPercent / 100.0);

	std::string barString = "[";
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos)
			barString += "=";
		else if (i == pos)
			barString += ">";
		else
			barString += " ";
	}
	barString += "]";

	return barString;
}

std::string ProgressBar::getPercentDoneString() {
	double percentProgress = std::min(100.0, (int)(currentPercent * 10 + 0.5) / 10.0);

	std::stringstream ss;
	ss << std::setw(5) << percentProgress << "%";

	return ss.str();
}

std::string ProgressBar::getTimeLeftString() {
	return durationToString(timeLeft);
}

std::string ProgressBar::getTimeTakenString() {
	return durationToString(Clock::now() - clockStart);
}

std::string ProgressBar::durationToString(duration time_left) {
	std::chrono::hours hours_left = std::chrono::duration_cast<std::chrono::hours>(time_left);
	std::chrono::minutes minutes_left = std::chrono::duration_cast<std::chrono::minutes>(time_left - hours_left);
	std::chrono::seconds seconds_left = std::chrono::duration_cast<std::chrono::seconds>(time_left - hours_left -  minutes_left);

	int h = hours_left.count();
	int m = minutes_left.count();
	int s = seconds_left.count();

	int d = h / 24;
	h = h - 24 * d;

	if (d > 999) {
		d = 999;
		h = 23;
		m = 59;
		s = 59;
	}

	std::stringstream ss;
	ss << std::setw(3) << d << "d:" << std::setw(2) << h << "h:" << std::setw(2) << m << "m:" << std::setw(2) << s << "s";
	
	return ss.str();
}

void ProgressBar::update(double currentVal) {
	double prevPercent = currentPercent;
	currentPercent = std::min(100.0, 100.0 * currentVal / maxVal); 

	if (currentPercent < 100) {
		duration speed = (Clock::now() - clockLast) / (currentPercent - prevPercent);
		if (isStarting) {
			isStarting = false;
			speedEMA = (duration)speed;
		}
		else
			speedEMA = (duration)(smoothing * speed + (1.0 - smoothing) * speedEMA);
		
		timeLeft = speedEMA * (100 - currentPercent);
	}

	clockLast = Clock::now();
}

void ProgressBar::end() {
	if (currentPercent < 100)
		update(1.001 * maxVal);
}
