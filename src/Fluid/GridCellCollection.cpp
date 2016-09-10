/*
 * GridCellCollection.cpp
 *
 *  Created on: 17 Jun 2015
 *      Author: harry
 */

#include "GridCellCollection.hpp"
#include "IO/Logger.hpp"

#include <iostream>

void GridCellCollection::start(const std::string& name) {
	std::map<std::string, std::pair<int, int>>::iterator it = guards.find(name);
	if (it == guards.end()) {
		guardsStarted.push_back(name);
		std::pair<int, int> newGuards;
		newGuards.first = cells.size();
		newGuards.second = cells.size();
		guards[name] = newGuards;
	}
}

int GridCellCollection::add() {
	cells.emplace_back();
	for (const std::string& name : guardsStarted)
		guards[name].second += 1;
	return cells.size()-1;
}

void GridCellCollection::stop(const std::string& name) {
	for (int i = guardsStarted.size() - 1; i >= 0; --i) {
		if (guardsStarted[i] == name)
			guardsStarted.erase(guardsStarted.begin() + i);
	}
}

std::vector<GridCell>& GridCellCollection::getCellVector() {
	return cells;
}

Looper GridCellCollection::getIterable(const std::string& name) {
	std::map<std::string, std::pair<int, int>>::iterator it = guards.find(name);
	if (it != guards.end()) {
		std::pair<int, int>& iterGuards = it->second;
		return Looper(cells, iterGuards.first, iterGuards.second);
	}
	else {
		Logger::Instance().print<SeverityType::WARNING>("GridCellCollection::getIteratable: invalid iterable: ", name, '\n');
		return Looper(cells, 0, 0);
	}
}

Looper GridCellCollection::getIterable() {
	return Looper(cells, 0, cells.size());
}

ConstLooper GridCellCollection::getIterable(const std::string& name) const {
	std::map<std::string, std::pair<int, int>>::const_iterator it = guards.find(name);
	if (it != guards.end()) {
		std::pair<int, int>& iterGuards = (std::pair<int, int>&)it->second;
		return ConstLooper(cells, iterGuards.first, iterGuards.second);
	}
	else {
		Logger::Instance().print<SeverityType::WARNING>("GridCellCollection::getIteratable: invalid iterable: ", name, '\n');
		return ConstLooper(cells, 0, 0);
	}
}

ConstLooper GridCellCollection::getIterable() const {
	return ConstLooper(cells, 0, cells.size());
}

ConstLooper::ConstLooper(const std::vector<GridCell>& cells, int startID, int endID)
: m_cells(cells)
, startID(startID)
, endID(endID)
{

}

std::vector<GridCell>::const_iterator ConstLooper::begin() const {
	return std::vector<GridCell>::const_iterator(m_cells.begin()+startID);
}
std::vector<GridCell>::const_iterator ConstLooper::cbegin() const {
	return std::vector<GridCell>::const_iterator(m_cells.cbegin()+startID);
}
std::vector<GridCell>::const_iterator ConstLooper::end() const {
	return std::vector<GridCell>::const_iterator(m_cells.begin()+endID);
}
std::vector<GridCell>::const_iterator ConstLooper::cend() const {
	return std::vector<GridCell>::const_iterator(m_cells.cbegin()+endID);
}

Looper::Looper(std::vector<GridCell>& cells, int startID, int endID)
: m_cells(cells)
, startID(startID)
, endID(endID)
{

}

std::vector<GridCell>::iterator Looper::begin() {
	return std::vector<GridCell>::iterator(m_cells.begin()+startID);
}
std::vector<GridCell>::const_iterator Looper::begin() const {
	return std::vector<GridCell>::const_iterator(m_cells.begin()+startID);
}
std::vector<GridCell>::const_iterator Looper::cbegin() const {
	return std::vector<GridCell>::const_iterator(m_cells.cbegin()+startID);
}
std::vector<GridCell>::iterator Looper::end() {
	return std::vector<GridCell>::iterator(m_cells.begin()+endID);
}
std::vector<GridCell>::const_iterator Looper::end() const {
	return std::vector<GridCell>::const_iterator(m_cells.begin()+endID);
}
std::vector<GridCell>::const_iterator Looper::cend() const {
	return std::vector<GridCell>::const_iterator(m_cells.cbegin()+endID);
}

