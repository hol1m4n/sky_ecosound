#include <iostream>
#include <iomanip>

void intensidades()
{
    ifstream data;
	data.open("decibels.txt");
	
	double value[1323000];
	
	
	for(int i=0;i<1323000;i++)
	{
		data>>value[i];
	}
	data.close();
	
	ofstream salida;
	salida.open("intensities.txt");
	
	for(int i=0;i<1323000;i++)
	{
		std::salida<<std::setprecision(1)<<convertir(value[i])<<std::endl;
	}
	
	salida.close();
	
	
	
	
}


double convertir(double db)
{
	
	return pow(10,db/10--3);
	
}