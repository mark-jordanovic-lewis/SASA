/**********************************************************************
  SASACalc Extension

  Copyright (C) 2011 by Mark Lewis
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef SASACALCEXTENSION_H
#define SASACALCEXTENSION_H

#include <avogadro/extension.h>
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

#include <QDialog>
#include <QAction>
#include <QMessageBox>
#include <QtCore/QDebug>
#include <QFileDialog>
#include <QDir>
#include <iostream>
#include <fstream>
#include <sstream>

#include "sasa_geometrySetup.h"
#include "sasa_calculation.h"
#include "ui_sasacalcfiledialog.h"

class SASACalcFileDialog: public QDialog, private Ui::SASACalcFileDialog
{
 public:
  
 private:
  //local variables
  sasa_geometrySetup geomSetup ;
  double result, probeRadius ;
  int nRuns, nTrials ;

};
	  
#endif
