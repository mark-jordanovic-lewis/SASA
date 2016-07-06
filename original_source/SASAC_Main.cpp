/* *********************************************************************
  SASAC_Text

  Copyright (C) 2011 by Mark Lewis

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ********************************************************************** */
 
#include "sasacalcextension.h"

using namespace std ;
  
int main(void) 
{
  char * running, path = './'  ;
  int trials, runs, probeType ;
  

  while(!strcmp(running, 'quit'))
  {
    cout << "Current file path: " << path << endl ;
    cout << "Open another directory? (y/n)" << endl ;
    if(!strcmp('y', cin())
    {
      cout << "input new dir: " ;
      getline(cin, path) ;
    }
    
    
  }
    
    
  spinBox_trials->setRange(1000, 1000000) ;                       							//manual set a couple of boxes
  spinBox_trials->setSingleStep(1000) ;
  spinBox_runs->setVisible(false) ;
  spinBox_runs->setMinimum(1) ;

  comboBox_probeType->addItem("ArgonAtom") ;
  comboBox_probeType->addItem("CabonDioxide") ;
  comboBox_probeType->addItem("CarbonMonoxide") ;
  comboBox_probeType->addItem("HydrogenMolecule") ;
  comboBox_probeType->addItem("NitrogenMolecule") ;
	  
  connect(pushButton_browse, SIGNAL(clicked()), this, SLOT(get_xyzPath())) ;  			//browse button to xyz_path
  connect(pushButton_browse, SIGNAL(clicked()), this, SLOT(populate_fileNames())) ;	//browse button to populate_filenames
  connect(lineEdit, SIGNAL(editingFinished()), this, SLOT(populate_fileNames())) ; 		//type in browser to populate_filenames

  connect(checkBox_runs, SIGNAL(stateChanged(int)), this, SLOT(multipleRuns())) ; 		//checkbox
		
  connect(comboBox_probeType, SIGNAL(currentIndexChanged(int)), this, SLOT(set_probeDiameter())) ;

  connect(pushButton_setup, SIGNAL(clicked()), this, SLOT(setup())) ;         			//setup button

  connect(pushButton_calculate, SIGNAL(clicked()), this, SLOT(calculate())) ; 			//calculate button

  connect(pushButton_quit, SIGNAL(clicked()), this, SLOT(quit())) ;          				//quit button
  nRuns = 1 ;
}


void SASACalcFileDialog::get_xyzPath() 														//opens file path to xyz (and possibly probe) files
{
  QString path = QFileDialog::getExistingDirectory(this, "Choose Directory", "/home", QFileDialog::ShowDirsOnly) ; //opens dialogue window
  lineEdit->setText(path) ;										//sett dir path to output on browser bar
}


void SASACalcFileDialog::populate_fileNames()													 //populates the comboBoxes
{
  comboBox_xyz->clear() ;                    							//clear entries in combo box when new dir
  comboBox_dat->clear() ;
  comboBox_ff->clear() ;
  QDir searchDIR(lineEdit->displayText()) ;   													//getting and filtering the files in the directory path
  if(searchDIR.exists())                      							//if directory exist
    {
      QStringList xyzFilter, datFilter, ffFilter ;                				//setup filters to find files
      xyzFilter << "*.xyz" ;
      datFilter << "*.dat" ;
      ffFilter << "*.atoms" ;
      QStringList xyzFiles = searchDIR.entryList(xyzFilter) ;    		//make lists of relevant files
      QStringList datFiles = searchDIR.entryList(datFilter) ;
      QStringList ffFiles = searchDIR.entryList(ffFilter) ;

      if(xyzFiles.empty()) SA_TextOutput->append("No relevant data files in " + lineEdit->displayText()) ; //if no report to text window
      else
	{
	  xyzFiles.sort() ;                                       					//sort files by name
	  datFiles.sort() ;
	  ffFiles.sort() ;
	  comboBox_xyz->addItems(xyzFiles) ;                      		//add to comboBoxes
	  comboBox_dat->addItems(datFiles) ;
	  comboBox_ff->addItems(ffFiles) ;
	  SA_TextOutput->append("Drop menus are populated.") ;   	//report to text window
	}
    }
  else                                                                        					//if there is no such directory on the system
    {
      QString invalid = lineEdit->displayText() + " is not a valid directory." ; 
      lineEdit->clear() ;                                                         
      SA_TextOutput->append(invalid) ;
    }
}


void SASACalcFileDialog::multipleRuns()							//reports status of multiple run check box THIS OUGHT TO BE FINISHED >:@
{
  if(checkBox_runs->checkState())
    {
      SA_TextOutput->append("Multiple runs on") ;
      spinBox_runs->setVisible(true) ;
    }
  else
    {
      SA_TextOutput->append("Multiple runs off") ;
      spinBox_runs->setVisible(false) ;
      spinBox_runs->setValue(1) ;
    }
}

void SASACalcFileDialog::set_probeDiameter() 												//sets up probe radii
{
  SA_TextOutput->append("probe change detected") ;
  int pT = comboBox_probeType->currentIndex() ;
  switch(pT)
    {
    case 0	: 	probeRadius = 1.88  ; break ;
    case 1 	: 	probeRadius = 1.128 ; break ;
    case 2	:	probeRadius = 1.163 ; break ;
    case 3	:	probeRadius = 1.2 ; break ;
    case 4	:	probeRadius = 1.82 ; break ;
    }
}
	
void SASACalcFileDialog::setup()
{
  SA_TextOutput->clear() ;
		
  QString path = lineEdit->displayText() + "/" ;
  QString xyzFile = path + comboBox_xyz->currentText() ;
  QString datFile = path + comboBox_dat->currentText() ;
  QString ffFile = path + comboBox_ff->currentText() ;
  cout << "files being accessed" << endl ;
  SA_TextOutput->append("Files being accessed: ") ;
  SA_TextOutput->append(xyzFile) ;
  SA_TextOutput->append(datFile) ;
  SA_TextOutput->append(ffFile) ;
  cout << "files being sent to setup" << endl ;
  geomSetup.Read(xyzFile.toStdString(), datFile.toStdString(), ffFile.toStdString()) ;
}
	

	
void SASACalcFileDialog::calculate()
{
  sasa_calculation casio ;
  QString path = lineEdit->displayText() + "/" ;
  if(checkBox_MD->isChecked())
    {
      stringstream outstream ;
      QString datFile = path + comboBox_dat->currentText() ;
      QString ffFile = path + comboBox_ff->currentText() ;
      int counter = 0, xyz_size = comboBox_xyz->count() ;
      nTrials = spinBox_trials->cleanText().toInt() ;
      nRuns = spinBox_runs->cleanText().toInt() ;
      outstream << path.toStdString() << "/SASAx" << nTrials << "of" << comboBox_probeType->currentText().toStdString() << "_output.rslt" ;
      ofstream MD_Results(outstream.str().c_str()) ;	
      result = 0 ;
			
      for(int calc = 0 ; calc < xyz_size ; calc++) 
	{
	  QString xyzFile = path + comboBox_xyz->itemText(calc) ;
	  geomSetup.Read(xyzFile.toStdString(), datFile.toStdString(), ffFile.toStdString()) ;
	  if(checkBox_runs->checkState()) 
	    {
	      casio.setParams(&nTrials, &probeRadius, &geomSetup.xyzPositions, &geomSetup.abcPositions, &geomSetup.Radii, &geomSetup.ortho_to_nonortho, &geomSetup.nonortho_to_ortho, &geomSetup.basisLengths, &geomSetup.basisAngles) ;
	      for(int i = 0 ; i < nRuns ; i++) result += casio.calculateSA() ;
	      result /= nRuns ;
	    }
	  else
	    {
	      casio.setParams(&nTrials, &probeRadius, &geomSetup.xyzPositions, &geomSetup.abcPositions, &geomSetup.Radii, &geomSetup.ortho_to_nonortho, &geomSetup.nonortho_to_ortho, &geomSetup.basisLengths, &geomSetup.basisAngles) ;
	      result = casio.calculateSA() ;
	    }
				
	  stringstream resultStream ;
	  resultStream << counter++ << " " << result ;
	  MD_Results << resultStream.c_str() << "\n" ;
	}
    }
  else
    {
      SA_TextOutput->append("Calculate was clicked") ;
      casio.setParams(&nTrials, &probeRadius, &geomSetup.xyzPositions, &geomSetup.abcPositions, &geomSetup.Radii, &geomSetup.ortho_to_nonortho, &geomSetup.nonortho_to_ortho, &geomSetup.basisLengths, &geomSetup.basisAngles) ;
      nTrials = spinBox_trials->cleanText().toInt() ;
      nRuns = spinBox_runs->cleanText().toInt() ;
      result = 0 ;
      for(int i = 0 ; i < nRuns ; i++) result += casio.calculateSA() ;
      result /= nRuns ;
      stringstream resultStream, res ;
      res << result ;
      resultStream << "Solvent accessible surface: " << res.str() << " Angstrom squared" ;
      QString resultString(resultStream) ;
      SA_TextOutput->append(resultString) ;
    }
}
  
void SASACalcFileDialog::quit(){SA_TextOutput->append("Quit was clicked") ;}

SASACalcFileDialog::~SASACalcFileDialog(){}
