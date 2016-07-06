#include <stdio>
#include <stdlib>
#include "sasacalcextension.h"
#include "sasa_calculation.h"
#include "sasa_tokenUtil.h"
#include "sasa_transformMatrix.h"
#include "sasa_RandomGenerator.h"

int main(argc, argv[])
{
  FILE * xyz_file, dat_file, atoms_file ;
  xyz_file = fopen('0987553_1704.xyz', 'r') ;
  dat_file = fopen('1millstep-1millop_.dat', 'r') ;
  atoms_file = fopen('1millstep-1millop_.atoms', 'r') ;
