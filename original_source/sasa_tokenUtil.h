// below code chops strings up into tokens, courtacy of Jens
#ifndef TOKENUTIL_H
#define TOKENUTIL_H

#include <iostream>
#include <vector>
#include <sstream>


class tokenUtil
{
private:
public:
		tokenUtil() ;
		bool tokenize(std::vector<std::string>&, const char*, const char*) ;
} ;
/**
 * This purloined from: http://www.codeguru.com/forum/showthread.php?t=231054
*/
template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
		std::istringstream iss(s) ;
  return !(iss >> f >> t).fail() ;
}
#endif
