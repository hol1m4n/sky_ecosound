#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sndfile.h>
#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <TH1F.h>
#include "TGraphErrors.h"
#include <TApplication.h>
#include <sstream>
#include <TFile.h>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include "TLegend.h"
#include "TProfile.h"
#include <TRandom1.h>



using namespace std;

Double_t convert_dB_to_Intensity(Double_t dB, Int_t b)
{
   return pow(10,dB/10-b);// dB=10log(la/lb) with la=sound intensity and lb=reference power = pow(10,-b) (Note: normally lb=threshold of human hearing 10⁻¹²W/m²)
}
Double_t convert_Intensity_to_dB(Double_t I, Int_t b)
{
   return 10*log(I/pow(10,-b));// dB=10log(la/lb) with la=sound intensity and lb=reference power = pow(10,-b) (Note: normally lb=threshold of human hearing 10⁻¹²W/m²)
}
Int_t
main (Int_t argc, char * argv [])
{//1 
    
    char         *infilename;
    stringstream ss;
    string s_infilename = "CanyonWren.wav";
    infilename = &s_infilename[0u];

  
    SNDFILE    *infile = NULL;
    SF_INFO    sfinfo;
    memset (&sfinfo, 0, sizeof (sfinfo)) ;


    //begin: reviewing that run parameters are OK
    if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
      {    
        printf ("Not able to open input file %s.\n", infilename) ;
        puts (sf_strerror (NULL)) ;
        return 1 ;
      }
    //end: reviewing that run parameters are OK
  

    ofstream output;
    output.open("WAV_in_text.txt");   

    ofstream output1;
    output1.open("FT_WAV.txt");   


    Int_t N=1024;
    float buf[N];
    Int_t readcount;         
    Int_t ref_power=-3;


    Double_t *in = new Double_t[2*(N/2+1)];
      
    Double_t *re = new Double_t[N/2+1];
    Double_t *im = new Double_t[N/2+1];

             
    Double_t *mag = new Double_t[N/2+1];
    Double_t *pha = new Double_t[N/2+1];

    Double_t mean_mag[N/2+1];

    for (int i=0; i<N/2+1; i++)
      {
         mean_mag[i]=0;
      }

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(111111111);
    gStyle->SetPalette(60);
    gStyle->SetOptFit(1);

    TApplication* rootapp = new TApplication("example",&argc, argv);
    TFile * outf = new TFile("Plots_Hollman.root", "RECREATE");
        
    int count=0;
 
    while ((readcount = sf_readf_float(infile, buf, N) )>0 ) 
      {
         count++;
         //Make our own TVirtualFFT object (using option "K")
         //Third parameter (option) consists of 3 parts:
         //-transform type:
         // real input/complex output in our case
         //-transform flag: 
         // the amount of time spent in planning
         // the transform (see TVirtualFFT class description)
         //-to create a new TVirtualFFT object (option "K") or use the global (default)
            
         TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &N, "R2C ES K");
         if (!fft_own) return 1;
         for (Int_t i=0; i<N; i++)
         {
            in[i] =  convert_dB_to_Intensity(buf[i],ref_power);
         }
         fft_own->SetPoints(in);
         fft_own->Transform();
         fft_own->GetPointsComplex(re,im);
         delete  fft_own;
         for (Int_t i=0; i<N/2+1;i++)
         {
           mag[i]=sqrt(pow(re[i],2)+pow(im[i],2));
           pha[i]=atan(im[i]/re[i]);
           mean_mag[i]+= mag[i];
         }
     
       for (int i=0; i<N; i++)
          {
            output<<buf[i]<<" "<<convert_dB_to_Intensity(buf[i], ref_power)<<endl; 
          }

       for (int i=0; i<N/2+1; i++)
          {
            output1<<mag[i]<<" "<<pha[i]<<endl; 
          }
      }


      for (int i=0; i<N/2+1; i++)
        {
          mean_mag[i]=mean_mag[i]/count; 
        }
    
    ifstream prueba;
    prueba.open("TF_Raven2.txt");
    Double_t tmp=0;


    TH1D * h_Magnitude;
    TH1D * h_Raven;
    Int_t nbins= N/2;
    h_Magnitude = new TH1D("Spectrum's Magnitude","Spectrum's Magnitude", nbins, 0, (nbins-1));
    h_Raven = new TH1D("Spectrum's Magnitude Raven","Spectrum's Magnitude Raven", nbins, 0, (nbins-1));
    for (Int_t i=0; i<nbins+1; i++)
      {
         prueba>>tmp;
         h_Raven->SetBinContent(i,tmp -86);    
         h_Magnitude->SetBinContent(i,convert_Intensity_to_dB(mean_mag[i],ref_power));      
      }

    TCanvas *c = new TCanvas("Canvas","Spectrum's Magnitudes Comparison",0,0,600,600);
    gStyle->SetOptStat(0);
    //cy->Divide(1,3);
    //cy->cd(1);
    h_Magnitude->SetLineColor(kBlue);
    h_Raven->SetLineColor(kRed);
    h_Magnitude->Draw();    
    h_Raven->Draw("same");
    c->Update();
    //gSystem->mkdir("Plots");
    //gSystem->mkdir(("Plots/"+month[m]).c_str());
    c->SaveAs("Raven.png");
    c->Clear();
    
    output.close();
    output1.close();
    outf->Write();
    return 0 ;
   
}//1
