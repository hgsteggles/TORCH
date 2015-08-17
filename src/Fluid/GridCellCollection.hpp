/*
 * GridCellCollection.hpp
 *
 *  Created on: 17 Jun 2015
 *      Author: harry
 */

#ifndef GRIDCELLCOLLECTION_HPP_
#define GRIDCELLCOLLECTION_HPP_

#include <map>
#include <string>
#include <vector>

#include "GridCell.hpp"

class Looper;
class ConstLooper;

class GridCellCollection {
public:
	void start(const std::string& name);
	void stop(const std::string& name);
	int add();
	std::vector<GridCell>& getCellVector();

	void addIndexOrder(const std::string& name, const std::vector<int>& order);

	Looper getIterable(const std::string& name);
	Looper getIterable();
	ConstLooper getIterable(const std::string& name) const;
	ConstLooper getIterable() const;

private:
	std::vector<GridCell> cells;
	std::vector<std::string> guardsStarted;
	std::map<std::string, std::pair<int, int>> guards;
};

class ConstLooper {
public:
	ConstLooper(const std::vector<GridCell>& cells, int startID, int endID);

	std::vector<GridCell>::const_iterator begin() const;
	std::vector<GridCell>::const_iterator cbegin() const;
	std::vector<GridCell>::const_iterator end() const;
	std::vector<GridCell>::const_iterator cend() const;
private:
	const std::vector<GridCell>& m_cells;
	int startID;
	int endID;
};

class Looper {
public:
	Looper(std::vector<GridCell>& cells, int startID, int endID);

	std::vector<GridCell>::iterator begin();
	std::vector<GridCell>::const_iterator begin() const;
	std::vector<GridCell>::const_iterator cbegin() const;
	std::vector<GridCell>::iterator end();
	std::vector<GridCell>::const_iterator end() const;
	std::vector<GridCell>::const_iterator cend() const;
private:
	std::vector<GridCell>& m_cells;
	int startID;
	int endID;
};


#endif /* GRIDCELLCOLLECTION_HPP_ */
