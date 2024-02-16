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

void elephant()
{
    double mag_v[513][88];
	
    ifstream rawdata;
	rawdata.open("animal_spec.txt");
    
    int t=0;
	
	while(t!=88)
	{
	
	for(int i=0;i<513;i++)
	{
		rawdata>>mag_v[i][t];
		//cout<<mag_v[i][t]<<endl;
	}
    t++;
	}
	rawdata.close();
	
	
	double rip[135432];
	double time[45144];
	double frequency[45144];
	double amplitude[45144];
	
	
	ifstream data;
	data.open("zoo_rawdata.txt");
	
	for(int i=0;i<135432;i++)
	{
        data>>rip[i];
	}
	data.close();
	
	
	for(int i=0;i<45144;i++)
	{
        frequency[i]=rip[i];
		time[i]=rip[i+45144];
		amplitude[i]=rip[i+90288];
	}
	

	plot = new TCanvas("Plot","Animal Spectrogram",600,400);
    tridihisto  = new TH2F("LEGO","Elephant Spectrogram 3D Representation",87,-1,87,512,-1,512);
   
   for(int i=0;i<45144;i++)
   {
	   tridihisto->Fill(time[i],frequency[i],amplitude[i]);
   }
    
  
	//gStyle->SetPalette(kBird);
    tridihisto->Draw("CONT4Z");
	//tridihisto ->Draw("LEGO2Z");
	//tridihisto ->Draw("SURF2");
	//tridihisto ->Draw("COLZ");
	//tridihisto ->Draw("scat=0.5");
	//tridihisto ->Draw("E");  //FOR MANIPULATION
	
	tridihisto->GetXaxis()->SetTitle("Time(s)-bs 1.95");
    tridihisto->GetYaxis()->SetTitle("Frequency(Hz)-bs 0.51");
	tridihisto->GetZaxis()->SetTitle("Amplitude(dB)");
	
	
	
}






