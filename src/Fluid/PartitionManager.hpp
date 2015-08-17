/*
 * PartitionManager.hpp
 *
 *  Created on: 17 Jun 2015
 *      Author: harry
 */

#ifndef PARTITIONMANAGER_HPP_
#define PARTITIONMANAGER_HPP_

#include "MPI/MPI_Wrapper.hpp"

class PartitionManager {
public:
	~PartitionManager();

	void initialise(int ncells);

	void resetBuffer();
	void addSendItem(double val);
	double getRecvItem();
	void sendData(int destination, SendID tag);
	void recvData(int source, SendID tag);
	void exchangeData(int destination, SendID tag);
	int getBufferCount();
	int getSendCount();
	int getRecvCount();
private:
	unsigned int nelements = 0;
	int m_bufferCount = 0;
	int m_recvCount = 0;
	int m_sendCount = 0;
	double* send_buffer = nullptr;
	double* recv_buffer = nullptr;
};



#endif /* PARTITIONMANAGER_HPP_ */
