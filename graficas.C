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

void grapher()
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

	/*
c1 = new TCanvas("Grafica","Grafica",1000,1000);
	c1->Divide(1,1);
	c1->cd(1);
    main_graph = new TGraph();
	for(int a=-1;a<2;a++) //Llenar
	{
	    main_graph->SetPoint(a+10,a,fdex(a));
		cout<<fdex(a)<<endl;
	}

	main_graph->SetLineColor(kBlue-6);
	main_graph->SetLineWidth(1);
	main_graph->Draw("AC*");
	c1->Update();
	*/
}

double fdex(double name) //Esta funcion evalua f(x) a partir de un numero (x) de entrada
{
	double value;
	double euler=2.718281828459045;
	value=(pow(euler,name))+(pow(name,3))+(2*(pow(name,2)))+(-1);
	return value;
}





