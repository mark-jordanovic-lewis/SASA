// Mark Lewis 20-10-11

#include <iostream>
#include "SASATypes.h"

int main()
{
  structSetup build ;
  struct geometry system ;
  int buildErr = build.geometrySetup(&system) ;
  
  
