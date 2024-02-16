#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <TH1F.h>
#include <TApplication.h>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TRandom1.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <iostream.h>
#include <fstream.h>
#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TGraphErrors.h>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TRandom1.h>
#include <sys/stat.h>


void first_function()
{
	int numero1;
	int numero2;
	cout<<"Por favor ingrese dos numeros para sumar"<<endl;
	cout<<"Numero 1"<<endl;
	cin>>numero1;
	cout<<"Numero2"<<endl;
	cin>>numero2;
	cout<<"La suma de "<<numero1<<" y "<<numero2<<" es: "<<numero1+numero2<<endl;
}

double fdex(double name)
{
	double value;
	double euler=2.718281828459045;
	value=(pow(euler,name))+(pow(name,3))+(2*(pow(name,2)))+(-1);
	return value;
}

double fpdex(double num)
{
	double value;
	double euler=2.718281828459045;
	value=(pow(euler,num))+(3*(pow(num,2)))+(4*num);
	return value;
}


void newt_iter()
{
	float ciclo=0;
	float inicio=0;
	cout<<"Ingrese el numero de iteraciones a realizar"<<endl;
	cin>>ciclo;
	cout<<"Ingrese el valor aproximado del valor de la funcion"<<endl;
	cin>>inicio;
	for(int a=0;a<ciclo;a++)
	{
		inicio=inicio-(fdex(inicio))/(fpdex(inicio));
	}
	cout<<"El resultado de "<<ciclo<<" iteraciones para la funcion fue "<<endl;
	std::cout << std::fixed;
    std::cout << std::setprecision(8);
    std::cout << inicio;
}

	
	
	
	
	
	
	
