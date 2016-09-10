
static double dummy_checkpoint;

class Checkpointer {
public:
	Checkpointer(double maxTime, int ncheckpoints);
	bool update(double currentTime, double& timeToNextCheckpoint = dummy_checkpoint);
	int getCount();
private:
	int checkpointCount;
	double checkpointDelta;
	double maxTime;
};
