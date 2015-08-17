/*
 * PartitionManager.cpp
 *
 *  Created on: 17 Jun 2015
 *      Author: harry
 */

#include "PartitionManager.hpp"

#include <iostream>
#include <stdexcept>


PartitionManager::~PartitionManager() {
	if (send_buffer)
		delete[] send_buffer;
	if (recv_buffer)
		delete[] recv_buffer;
}

void PartitionManager::initialise(int ncells) {
	nelements = ncells;
	send_buffer = new double[nelements];
	recv_buffer = new double[nelements];
}

void PartitionManager::resetBuffer() {
	m_bufferCount = 0;
	m_recvCount = 0;
	m_sendCount = 0;
}

void PartitionManager::addSendItem(double val) {
	if (m_bufferCount < (int)nelements)
		send_buffer[m_sendCount++] = val;
}

double PartitionManager::getRecvItem() {
	if (m_recvCount < m_bufferCount)
		return recv_buffer[m_recvCount++];
	else
		throw std::runtime_error("PartitionManager::getRecvItem(): trying to receive item that doesn't exist.");
}

void PartitionManager::exchangeData(int destination, SendID tag) {
	MPIW::Instance().exchange(send_buffer, recv_buffer, m_sendCount, destination, tag);
	m_bufferCount = m_sendCount;
	m_recvCount = 0;
	m_sendCount = 0;
}

void PartitionManager::sendData(int destination, SendID tag) {
	MPIW::Instance().send(&m_sendCount, 1, destination, tag);
	MPIW::Instance().send(send_buffer, m_sendCount, destination, tag);
	m_sendCount = 0;
	m_recvCount = 0;
	m_bufferCount = 0;
}

void PartitionManager::recvData(int source, SendID tag) {
	m_recvCount = 0;
	MPIW::Instance().receive(&m_bufferCount, 1, source, tag);
	MPIW::Instance().receive(recv_buffer, m_bufferCount, source, tag);
}

int PartitionManager::getBufferCount() {
	return m_bufferCount;
}

int PartitionManager::getSendCount() {
	return m_sendCount;
}

int PartitionManager::getRecvCount() {
	return m_recvCount;
}
