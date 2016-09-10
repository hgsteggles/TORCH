/*
 * ParseLua.hpp
 *
 *  Created on: 5 Mar 2015
 *      Author: harry
 */

#ifndef PARSELUA_HPP_
#define PARSELUA_HPP_

#include <string>

#include "Logger.hpp"

#include "selene/include/selene.h"

std::string parseString(sel::Selector selector) {
	std::string str = selector;
	return str;
}

bool exists(sel::Selector selector) {
	std::string str = selector;

	return !(str.compare("") == 0);
}

template<typename T>
void parseLuaVariable(sel::Selector selector, T& var) {
	if (exists(selector))
		var = (T)selector;
	else
		Logger::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default: ", var, '\n');
}

template<>
void parseLuaVariable(sel::Selector selector, std::string& var) {
	if (exists(selector))
		var = parseString(selector);
	else
		Logger::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default: ", var, '\n');
}


#endif /* PARSELUA_HPP_ */
