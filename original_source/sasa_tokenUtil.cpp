#include "sasa_tokenUtil.h"

using namespace std ;

tokenUtil::tokenUtil()
{}

bool tokenUtil::tokenize(vector<string> &vcr, const char *buf, const char *delimstr)
{
	vcr.clear() ;
  if (!buf || !delimstr) return false ;

  string s(buf) ;
  s += " " ;
  size_t startpos=0, endpos=0 ;

  for(;;)
  {
    startpos = s.find_first_not_of(delimstr,startpos) ;
    endpos = s.find_first_of(delimstr,startpos) ;

		if(endpos <= s.size() && startpos <= s.size())vcr.push_back(s.substr(startpos, endpos-startpos)) ;
    else break ;
    startpos = endpos+1 ;
  }
  return(true) ;
}