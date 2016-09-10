#include "Checkpointer.hpp"

#include <algorithm>
#include <iostream>

Checkpointer::Checkpointer(double maxTime, int ncheckpoints)
	: checkpointCount(0)
	, checkpointDelta(100.0 / ncheckpoints)
	, maxTime(maxTime)
{

}

bool Checkpointer::update(double currentTime, double& timeToNextCheckpoint) {
	double currentPercent = std::min(100.0, 100.0 * currentTime / maxTime);

	bool checkpointPassed = false;
	while (currentPercent - checkpointCount * checkpointDelta > checkpointDelta - 1.0e-8) {
		checkpointCount += 1;
		checkpointPassed = true;
	}

	timeToNextCheckpoint = ((checkpointCount + 1) * checkpointDelta - currentPercent) * maxTime / 100.0;

	return checkpointPassed;
}

int Checkpointer::getCount() {
	return checkpointCount;
}
