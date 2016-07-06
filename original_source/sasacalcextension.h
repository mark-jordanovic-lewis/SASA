/**********************************************************************
  SASACalc Extension

  Copyright (C) 2011 by Mark Lewis
  
  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.openmolecules.net/>

  Some code is based on Open Babel
  For more information, see <http://openbabel.sourceforge.net/>

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

namespace Avogadro {
	class SASACalcFileDialog: public QDialog, private Ui::SASACalcFileDialog
	{
		// We inherit from QDialog so we need the Q_OBJECT Macro to set up the slot and signal functionality
		Q_OBJECT

	public:
	    //! Constructor
		explicit SASACalcFileDialog( QWidget *parent = 0, Qt::WindowFlags f = 0 ) ;
		~SASACalcFileDialog() ;
	private:
		//locals
		sasa_geometrySetup geomSetup ;
		double result, probeRadius ;
		int nRuns, nTrials ;

	public slots:
		void get_xyzPath() ;
		void populate_fileNames() ;
		void multipleRuns() ;
		void setup() ;
		void calculate() ;
		void quit() ;

	private slots:
		void set_probeDiameter() ;
	};

	  
	class SASACalcExtension : public Extension
	{
		Q_OBJECT
		AVOGADRO_EXTENSION("SASACalc", tr("SASACalc"), tr("SASACalc extension"))

	public:
		//! Constructor
		SASACalcExtension(QObject *parent=0)  ;
		//! Deconstructor
		virtual ~SASACalcExtension() ;

		virtual QList<QAction *> actions() const ; 
		virtual QString menuPath(QAction *action) const ;

		virtual QDockWidget * dockWidget() ;
		virtual QUndoCommand* performAction(QAction *action, GLWidget *widget) ;

		virtual void setMolecule(Molecule *molecule) ;

	private:
		QList<QAction *> m_actions ;
		Molecule *m_molecule ;
		SASACalcFileDialog *m_dialog ;
	};

	
	class SASACalcExtensionFactory : public QObject, public PluginFactory
	{
		Q_OBJECT
		Q_INTERFACES(Avogadro::PluginFactory)
		AVOGADRO_EXTENSION_FACTORY(SASACalcExtension)
	};

} // end namespace Avogadro

#endif
