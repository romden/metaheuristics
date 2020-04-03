#include "Settings.h"


CSettings::CSettings(void)
{
	initialization();
}


CSettings::~CSettings(void)
{
}


void CSettings::initialization(void)
{
	// READ PARAMETER SETTINGS FROM FILE

	std::ifstream fin; // create ifstream object called fin	
	fin.open("param_settings.dat"); // opens the file

	std::string temp; // temp string
	fin >> temp >> CR; // DE crossover probability
	fin >> temp >> F; //  scaling factor 
	fin >> temp >> pm; // mutation probability 
	fin >> temp >> etam; // mutation distribution index
	fin >> temp >> delta; // probability for mating pool
	fin >> temp >> T; // neighborhood size
	fin >> temp >> p; // parameter to compute Lp norm
	fin >> temp >> K; // number of closest to check dominance
	fin >> temp >> mu; // population size
	fin >> temp >> maxEval; // stopping criterion	
	
	fin.close(); // close the file

}


//#include <iostream>
	//std::cout << CR << std::endl;
	//std::cout << F << std::endl;
	//std::cout << pm << std::endl;
	//std::cout << etam << std::endl;
	//std::cout << delta << std::endl;
	//std::cout << T << std::endl;
	//std::cout << mu << std::endl;
	//std::cout << maxEval << std::endl;

	//CR = 1.0;
	//F = 0.5;
	//pm = 1.0/30.0;
	//etam = 20.0;
	//delta = 0.9;
	//T = 20;
	////nr = 2;
	//mu = 300;
	//maxEval = 3E+4;