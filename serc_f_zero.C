
//Se incluyen estas librerias para realizar cosas basicas tales como ingreso y salida de datos u operaciones matematicas.
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>


double absoluto(double num)
{
	// Esta funcion entrega el valor absoluto de un numero dado.
	// Mas adelante sera de utilidad para evaluar un loop. 
	
	if(num>0)
	{
		return num; //Si el numero es mayor a cero retorna el mismo numero.
	}
	if(num<0)
	{
		return num*-1;//Si el numero es menor a cero retorna el numero multiplicado por -1.
	}
    if(num==0)
	{
		return 0;//Si el numero es igual a cero retorna cero.
	}
}

double selector(double num1, double num2)
{
	//Esta funcion compara si la diferencia entre dos numeros es mayor a 0.001.
	//De ser afirmativo lo anterior entregara el numero mayor. 
	
	double inlet[3]={0,0,0}; // Se crea este arreglo apra guardar lo numeros que se ingresan.
	inlet[0]=num1; //Valor original del Primer numero que entra.
	inlet[1]=num2; //Valor original del Segundo numero que entra.
	//Convertirlos a valor absoluto para hallar la diferencia.
	num1=absoluto(num1);
	num2=absoluto(num2);

	inlet[2]=num1-num2;
	if(inlet[2]<0)
	{
		inlet[2]=inlet[2]*(-1);
	}
	if(inlet[2]>0.01)
	{
		return inlet[0]; //Si la dierencia es mayor a 0.01 retorna el primer numero de lo contratio entrega 0.
	}
	else
	{
		return 0;
	}
}
	
double fdex(double name) //Esta funcion evalua f(x) a partir de un numero (x) de entrada
{
	double value;
	double euler=2.718281828459045;
	value=(pow(euler,name))+(pow(name,3))+(2*(pow(name,2)))+(-1);
	return value;
}

double fpdex(double num)//Esta funcion evalua f'(x) (la derivada) a partir de un numero (x) de entrada
{
	double value;
	double euler=2.718281828459045;
	value=(pow(euler,num))+(3*(pow(num,2)))+(4*num);
	return value;
}

void main_func() //Esta es la funcion principal del codigo, los datos que el usuario debera propocionar seran los limites de evaluacion (el dominio) y el numero de iteraciones del metodo de Newton. 
{
	int sugg=0; // En esta variable se guarda el numero de iteraciones que tendran en el metodo de Newton.
	int lim[2]={0,0}; //Arreglo que guarda el limite de evaluacion y el intervalo dentro del dominio. 
	cout<<"Ingrese el limite de evaluacion de la recta (valor negativo)"<<endl;
	cin>>lim[0]; //Valor limite de evaluacion
	cout<<"Ingrese el numero de iteraciones a realizar por cada zona o valor estimado de X"<<endl;
	cin>>sugg; //Numero de iteraciones del metodo Newton Rapshon.
	lim[1]=((absoluto(lim[0])*2)/0.001)+1;  //Este es el intervalo de evaluacion de los posibles ceros en la recta espaciados cada 0.001
	double aux=0; //Variable auxiliar, se utilizara mas adelante
	double val=0; //En esta variable se guarda temporalmente el valor de f(x) de una X determinada espaciada 0.001.
	int selec=0; //Esta variable cuenta el numero de valores en que f(x) se aproxima a cero.
	double classy[10]={0,0,0,0,0,0,0,0,0,0}; //En este arreglo se guardan la cantidad de valores que se aproximan lo suficiente a cero.	
	double chosen[10]={0,0,0,0,0,0,0,0,0,0}; //En este arreglo se guardan los valores definitivos que se usaran en la iteracion de Newton.
	
	std::cout << std::fixed; //Estas dos lineas se encargan de fijar el numero de cifras significativas que se muestran en pantalla. 
    std::cout << std::setprecision(5); //Para este caso es 5
	
	cout<<"Zonas o valores del eje X donde f(y) se aproxima a cero"<<endl;
	cout<<"--------------------------"<<endl;
	
	
	//Bucle 1
	//Este bucle se encarga de evaluar la funcion en el intervalo establecido anteriormente cada 0.001 veces en el eje X.
	for(int b=0;b<lim[1];b++) 
	{
		aux=lim[0]+(b*(0.001)); //La variable auxiliar se encarga de establecer X segun el intervalo dado y el espacio entre ellos.
		val=fdex(aux); //Aqui se calcula f(x) con respecto al valor de la variable aux. 
		if((val<=0.001)&&(val>=(-0.001))) // Si el valor de f(x) se encuentra entre 0.001 y -0.001 se toma como criterio de seleccion al momento de iterar con el metodo Newton.
		{
			selec++; // Si se encuentra un valor cercano a cero se tiene en cuenta y se numera.
			classy[selec-1]=aux; //Se guarda el presunto valor de X para el cual f(x) es igual a cero en el arreglo Classy.
			// Se muestran en pantalla los presuntos valores de X para los cuales f(x) es cero.
			
			cout<<"Numero "<<selec<<": "<<classy[selec-1]<<endl; //Cada vez que se calcule un posible candidato se mostrara en pantalla. 
		}	
	}
	
	cout<<"--------------------------"<<endl;
	
	//Una vez identificadas las zonas o los valores donde f(x) se aproxima a cero en X, se escogen valores significativos en cada una de las zonas o bien el numero calculado.
	// El siguiente bucle se encarga de realizar esa tarea repitiendose el numero de veces como valores hayan sido identificados.
	
	double senal[2]={0,0}; //En este arreglo se guardan los valores seleccionados respecto a la tarea anterior.
	int lax=0; //Variable contador, su utilidad se vera mas adelante.
	
	
	//Bucle 2
	for(int b=0;b<selec;b++)
	{
		//Si el valor entre los numeros es mayor a 0.001 entonces se guarda en el arreglo chosen, de lo contrario se promedia y se sigue con el bucle.
		senal[0]=selector(classy[b],classy[b+1]); 
		senal[1]=absoluto(senal[0]);
		if(senal[1]>0) //Este condicional esta dado con base a la funcion selector.
		{
			lax++;
			chosen[lax-1]=senal[0]; //Si la diferencia entre los numeros es mayor a 0.001 el numero se guarda, de lo contrario se suma al siguiente valor y se promedia.
		}
		else
		{
			classy[b+1]=(classy[b]+classy[b+1])/2;
		}
	}
	
	
	
	std::cout << std::fixed; //Estas dos lineas se encargan de fijar el numero de cifras significativas que se muestran en pantalla. 
    std::cout << std::setprecision(30); //Para e
	

	cout<<"Valores calculados donde F(x) es igual a cero"<<endl;
	cout<<"---------------------------------------------"<<endl;
	cout<<"Valor de X	Valor de F(x)	Error absoluto"<<endl;
	
	//Cuando se tienen los posibles valores frontera para alimentar las iteraciones del metodo de Newton se procede a iterar.
	
	//Este condicional indica que si algun valor del arreglo chosen es 0 o f(x) es cero cuando X es cero entonces muestre una fila de ceros.
	if((fdex(0)==0)||(chosen[0]==0)||(chosen[1]==0)||(chosen[2]==0)||(chosen[3]==0)||(chosen[4]==0)||(chosen[5]==0)||(chosen[6]==0)||(chosen[7]==0)||(chosen[8]==0)||(chosen[9]==0))
	{
		cout<<"0	0	0"<<endl;
	}
	
	
	// Si alguno de los valores de la variable chosen es mayor a cero entonces se procede a iterar con el metodo de Newton. 
	if((chosen[0]!=0)||(chosen[1]!=0)||(chosen[2]!=0)||(chosen[3]!=0)||(chosen[4]!=0)||(chosen[5]!=0)||(chosen[6]!=0)||(chosen[7]!=0)||(chosen[8]!=0)||(chosen[9]!=0))
	{
		//Bucle 3
		for(int c=0;c<10;c++) // Se hace el numero de veces que tenga el arreglo chosen.
		{
			//Bucle 3.1
			for(int a=0;a<sugg;a++)
			{
				chosen[c]=chosen[c]-(fdex(chosen[c]))/(fpdex(chosen[c])); //Este bucle realiza iteraciones igual a la variable sugg de cada uno de los valores frontera escogidos en el bucle 2. 
			}
			
			if(absoluto(chosen[c])>0)
			{
				cout<<chosen[c]<<"	"<<fdex(chosen[c])<<"	"<<absoluto(0-fdex(chosen[c]))<<endl; // Mostrar en pantalla los valores con su respectivo error. 
			}
		}
	}
	
	
	int n=100;
	double xvalues[100],yvalues[100];
	double xchosen[3]={-1.723,-0.528,0},ychosen[3]={0,0,0};
	
	for(int i=0;i<n;i++)
	{
		xvalues[i]=-5+(0.1*(i+1));
		yvalues[i]=fdex(-5+(0.1*(i+1)));
	}
	
	c1 = new TCanvas("Grafica","Grafica",1000,1000);
	c1->SetGrid();
	
	mg = new TMultiGraph("mg","f(x)=(e^x)+(x^3)+(2x^2)-1");
	
    gr1 = new TGraph(n,xvalues,yvalues);
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1);
    gr1->GetYaxis()->SetRangeUser(-2,2);
    gr1->GetXaxis()->SetRangeUser(-2,1);

	gr2 = new TGraph(3,xchosen,ychosen);
	
    
	mg->Add(gr1,"ACP");
    mg->Add(gr2,"P");
    
	gr1->SetMarkerSize(1);
	
	gr2->SetMarkerStyle(24);
	gr2->SetMarkerSize(5);
	gr2->SetMarkerColor(4);
	
    mg->Draw("AP");
	mg->GetXaxis()->SetTitle("x");
    mg->GetYaxis()->SetTitle("f(x)");

	gPad->Modified();
	mg->GetXaxis()->SetLimits(-2,1);
	mg->SetMinimum(-2);
	mg->SetMaximum(2);
	
	
}

void grapher() // Esta función simplemente se encarga de plotear los resutados en ROOT
{
	int n=100;
	double xvalues[100],yvalues[100];
	double xchosen[3]={-1.723,-0.528,0},ychosen[3]={0,0,0};
	
	for(int i=0;i<n;i++)
	{
		xvalues[i]=-5+(0.1*(i+1));
		yvalues[i]=fdex(-5+(0.1*(i+1)));
	}
	
	c1 = new TCanvas("Grafica","Grafica",1000,1000);
	c1->SetGrid();
	
	mg = new TMultiGraph("mg","f(x)=(e^x)+(x^3)+(2x^2)-1");
	
    gr1 = new TGraph(n,xvalues,yvalues);
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1);
    gr1->GetYaxis()->SetRangeUser(-2,2);
    gr1->GetXaxis()->SetRangeUser(-2,1);

	gr2 = new TGraph(3,xchosen,ychosen);
	
    
	mg->Add(gr1,"ACP");
    mg->Add(gr2,"P");
    
	gr1->SetMarkerSize(1);
	
	gr2->SetMarkerStyle(24);
	gr2->SetMarkerSize(5);
	gr2->SetMarkerColor(4);
	
    mg->Draw("AP");
	mg->GetXaxis()->SetTitle("x");
    mg->GetYaxis()->SetTitle("f(x)");

	gPad->Modified();
	mg->GetXaxis()->SetLimits(-2,1);
	mg->SetMinimum(-2);
	mg->SetMaximum(2);

}





