#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sndfile.h>
#include "TH1D.h"
#include "TH2D.h"
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
        Int_t ref_power=-3;
        Int_t Na=1024;//Number of points used for FT of audio data. 
        Int_t Na1=1024;//Normally is equal to Na, the only difference is when it is the last buffer of one point and you have to skip some values in order to maintain the number of points. 
        Int_t Ni=4096;//Number of points used for FT of indexes data. Note:To choose as high as possible this value.  
        Double_t scale_y=1; //must be sqrt(N);
        Double_t scale_x=1; //must be ((Double_t)(96000))/((Double_t)(N));

         
        
        TApplication* rootapp = new TApplication("example",&argc, argv);
        TFile * outf = new TFile("Plots.root", "RECREATE");
        char         *infilename;


        float buf[Na];  

        Int_t readcount=0;

        Int_t numberOfHours=71;//must be the same as the number of files, each file must be per hour. Doesn't matter if it is the total hour in order to run the code, however the periodicity is better analyzed if there are no holes.

        Int_t nbina= Na/2;        

        Int_t numberOfIndexes=2;
        TH1D * h_FT_Audio[Ni+1];//one histogram for each point. Note: the number of points is nbini+1.
        TH1D * h_Indexes[numberOfIndexes];
        TH1D * h_FT_Indexes[numberOfIndexes];


        string index[numberOfIndexes];
        index[0]="RMS";
        index[1]="Entropy";

        string file[numberOfHours];
        file[0]="070319_09.wav";
        file[1]="070319_10.wav";
        file[2]="070319_11.wav";
        file[3]="070319_12.wav";
        file[4]="070319_13.wav";
        file[5]="070319_14.wav";
        file[6]="070319_15.wav";
        file[7]="070319_16.wav";
        file[8]="070319_17.wav";
        file[9]="070319_18.wav";
        file[10]="070319_19.wav";
        file[11]="070319_20.wav";
        file[12]="070319_21.wav";
        file[13]="070319_22.wav";
        file[14]="070319_23.wav";
        file[15]="080319_00.wav";
        file[16]="080319_01.wav";
        file[17]="080319_02.wav";
        file[18]="080319_03.wav";
        file[19]="080319_04.wav";
        file[20]="080319_05.wav";
        file[21]="080319_06.wav";
        file[22]="080319_07.wav";
        file[23]="080319_08.wav";
        file[24]="080319_09.wav";
        file[25]="080319_10.wav";
        file[26]="080319_11.wav";
        file[27]="080319_12.wav";
        file[28]="080319_13.wav";
        file[29]="080319_14.wav";
        file[30]="080319_15.wav";
        file[31]="080319_16.wav";
        file[32]="080319_17.wav";
        file[33]="080319_18.wav";
        file[34]="080319_19.wav";
        file[35]="080319_20.wav";
        file[36]="080319_21.wav";
        file[37]="080319_22.wav";
        file[38]="080319_23.wav";
        file[39]="090319_00.wav";
        file[40]="090319_01.wav";
        file[41]="090319_02.wav";
        file[42]="090319_03.wav";
        file[43]="090319_04.wav";
        file[44]="090319_05.wav";
        file[45]="090319_06.wav";
        file[46]="090319_07.wav";
        file[47]="090319_08.wav";
        file[48]="090319_09.wav";
        file[49]="090319_10.wav";
        file[50]="090319_11.wav";
        file[51]="090319_12.wav";
        file[52]="090319_13.wav";
        file[53]="090319_14.wav";
        file[54]="090319_15.wav";
        file[55]="090319_16.wav";
        file[56]="090319_17.wav";
        file[57]="090319_18.wav";
        file[58]="090319_19.wav";
        file[59]="090319_20.wav";
        file[60]="090319_21.wav";
        file[61]="090319_22.wav";
        file[62]="090319_23.wav";
        file[63]="100319_00.wav";
        file[64]="100319_01.wav";
        file[65]="100319_02.wav";
        file[66]="100319_03.wav";
        file[67]="100319_04.wav";
        file[68]="100319_05.wav";
        file[69]="100319_06.wav";
        file[70]="100319_07.wav";

        


        
        for (Int_t p=0; p<Ni; p++)//the number of points is defined with Ni
          {
             stringstream ss;
             ss<<p;
             string str=ss.str();
             h_FT_Audio[p] = new TH1D(("FT_Audio for Point = "+ str).c_str(), ("FT_Audio for Point = "+ str).c_str(), nbina/2, 0, (nbina/2-1));
             /*h_FT_Signal[b]->SetOption("");//Surf3 Box1 Text");
             h_FT_Signal[b]->SetFillColor(42);
             h_FT_Signal[b]->SetStats(kFALSE);*/
          }

        for (Int_t i=0; i<numberOfIndexes; i++)
          {
               h_Indexes[i] = new TH1D(("Index ("+ index[i]+")").c_str(), ("Index ("+ index[i]+")").c_str(), Ni, 0, Ni-1);
               /*h_FT_Indexes[h]->SetOption("");//Surf3 Box1 Text");
               h_FT_Indexes[h]->SetFillColor(42);
               h_FT_Indexes[h]->SetStats(kFALSE);*/
          }

        
        

        

        Double_t *full_mag = new Double_t[nbina+1];
        Double_t *in = new Double_t[2*(nbina+1)];
      
        Double_t *re = new Double_t[nbina+1];
        Double_t *im = new Double_t[nbina+1];

             
        Double_t *mag = new Double_t[nbina+1];
        //Double_t *pha = new Double_t[nbina+1];


        for (Int_t b=0; b<nbina+1;b++)
         {
           full_mag[b]=0;
         }


       Int_t p=0; //number of the point 
      
    
       Double_t framesInOnePoint= ((Double_t)numberOfHours*60*60*96000)/((Double_t)Ni);
       cout<<"Frames in one point="<<framesInOnePoint<<endl;
       Int_t frames_in_this_point=0;
            

       for (Int_t h=0; h<numberOfHours; h++) //the number of files is equal to the number of hours by free decision.
        {//2 
            SNDFILE    *infile = NULL;///
            SF_INFO    sfinfo;///
            memset (&sfinfo, 0, sizeof (sfinfo)) ;///
            
        
            
            stringstream ss;
            ss<<file[h];
            string str=ss.str();
            infilename = &str[0u];
            //begin: reviewing that run parameters are OK
            if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
            {    
               printf ("Not able to open input file %s.\n", infilename) ;
               puts (sf_strerror (NULL)) ;
               return 1 ;
            }
            //end: reviewing that run parameters are OK
            cout<<"FILE: "<<str<<" PUNTO: "<<p<<" FRAME: "<<frames_in_this_point<<endl;
            cout<<"file size="<<sfinfo.frames<<endl;

            if (sfinfo.frames != (60*60*96000)) 
            {
              cout<<"There is not match with the number of frames in one hour. Please review it!"<<endl;
              //return 1;
            }

            //current_length=sfinfo.frames;
            int tmp=0; 
            while  ((readcount = sf_readf_float(infile, buf, Na1) )>0 ) 
            {//3
              tmp+=readcount;
              
              //if (p==Ni) cout<<tmp<<" "<<readcount<<" "<<frames_in_this_point<<endl;
              if (readcount == Na)  //In order to only include streams with the same lenght,  because maybe the last one has a different one.
                {//4    
                  ///begin:transform
                  /*for (Int_t b=0; b<Na; b++)
                  {
                     for (Int_t ch=0; ch<2; ch++)
                     {
                       if (ch==1) in[b] =  convert_dB_to_Intensity(buf[b * 2 + ch],ref_power);
                       cout<<"in: "<<in[b]<<endl;
                       cout<<"buf: "<<buf[b * 2 ]<<endl;
                       cout<<"buf2: "<<buf[b * 2 + 1 ]<<endl;
                     }
                  }*///for stereo files
                  frames_in_this_point=frames_in_this_point + readcount;
                  for (Int_t b=0; b<Na; b++)
                  {
                    in[b] =  convert_dB_to_Intensity(buf[b],ref_power);
                  }
                  //Make our own TVirtualFFT object (using option "K")
                  //Third parameter (option) consists of 3 parts:
                  //-transform type:
                  // real input/complex output in our case
                  //-transform flag: 
                  // the amount of time spent in planning
                  // the transform (see TVirtualFFT class description)
                  //-to create a new TVirtualFFT object (option "K") or use the global (default)
            
                  TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &Na, "R2C ES K");
                  if (!fft_own) return 1;
                  fft_own->SetPoints(in);
                  fft_own->Transform();
                  fft_own->GetPointsComplex(re, im);
                  delete  fft_own;
                  for (Int_t b=0; b<nbina+1;b++)
                  {
                    mag[b]=sqrt(pow(re[b],2)+pow(im[b],2));
                    full_mag[b]+=mag[b];
                    //cout<<"RE : "<<re[b]<<" IM: "<<im[b]<<endl;
                    //pha=atan(im[b]/re[b]);
                  }
                  if (p>Ni-1) cout<<"saliendo de TF con p="<<p<<" y tmp="<<tmp<<endl;
                }//4

              if ((framesInOnePoint - frames_in_this_point)>=0 && (framesInOnePoint - frames_in_this_point)<Na-1) //&& p<Ni) //frames in one point  
              {
               Na1=framesInOnePoint - frames_in_this_point;
               cout<<"Na1="<<Na1<<endl;
               //if (p=Ni/2) cout<<"Frames: "<<frames_in_this_point<<" "<<framesInOnePoint<<endl;
               Double_t entropy=0;
               Double_t area_freq=0;
               for (Int_t b=0; b<nbina/2;b++)//nbina+1
                {
                  h_FT_Audio[p]->SetBinContent(b,convert_Intensity_to_dB(full_mag[b],ref_power));
                  //cout<<"MAG: "<<full_mag[b]<<endl;
                  area_freq+=full_mag[b];
                }
                h_Indexes[0]->SetBinContent(p,h_FT_Audio[p]->GetRMS()); ////acá se debe llenar h_indexes para los otros índices.       
                for (Int_t b=0; b<nbina/2;b++)//nbina+1
                {
                  entropy+=(-full_mag[b]/area_freq*log(full_mag[b]/area_freq)/log(2));  
                }
                                
                h_Indexes[1]->SetBinContent(p,entropy);
                cout<<"Punto: "<<p<<" Total frames: "<<tmp<<endl;
                //cout<<"RMS: "<<h_FT_Audio[p]->GetRMS()<<endl;                   
                p++;
                if (p>Ni-2) cout<<"PUNTO: "<<p<<" FRAMES: "<<frames_in_this_point<<" Total frames="<<tmp<<endl; 
                for (Int_t b=0; b<nbina+1;b++)
                 {
                   full_mag[b]=0;
                 } 
                frames_in_this_point=0;
              }
              else 
              {
                Na1= Na;
              }

              //if (p>Ni-1) cout<<"volviendo al while con p="<<p<<" y tmp="<<tmp<<endl;

            }//3
               
            cout<<"saliendo del while"<<endl;
            sf_close (infile);
            cout<<"cerrando archivo"<<endl;

        }//2

        cout<<"SALIO"<<endl;
                             
               
       delete [] full_mag;
       delete [] in;
       delete [] re;
       delete [] im;
       delete [] mag;
       //delete [] pha;
        
        Int_t nbini= Ni/2;        

        for (Int_t i=0; i<numberOfIndexes; i++)
          {
             h_FT_Indexes[i] = new TH1D(("FT_Index ("+ index[i]+")").c_str(), ("FT_Index ("+ index[i]+")").c_str(), nbini, 0, (nbini-1));
          }


        Double_t *full_mag_i = new Double_t[nbini+1];
        Double_t *in_i = new Double_t[2*(nbini+1)];
      
        Double_t *re_i = new Double_t[nbini+1];
        Double_t *im_i = new Double_t[nbini+1];

             
        Double_t *mag_i = new Double_t[nbini+1];
        //Double_t *pha_i = new Double_t[nbini+1];

        
        for (Int_t i=0; i<numberOfIndexes; i++)
        {
          cout<<i<<endl;
          for (Int_t bi=0; bi<p; bi++)
          {
            in_i[bi] =  h_Indexes[i]->GetBinContent(bi);
          }
          cout<<"cargo info para TF"<<endl;
          TVirtualFFT *fft_own_i = TVirtualFFT::FFT(1, &Ni, "R2C ES K");
          if (!fft_own_i) return 1;
          fft_own_i->SetPoints(in_i);
          fft_own_i->Transform();
          fft_own_i->GetPointsComplex(re_i, im_i);
          delete  fft_own_i;
          for (Int_t bi=0; bi<nbini+1;bi++)
           {
             mag_i[bi]=sqrt(pow(re_i[bi],2)+pow(im_i[bi],2));
             h_FT_Indexes[i]->SetBinContent(bi,mag_i[bi]);               
           }
          cout<<"paso"<<endl;
        }

          delete [] in_i;
          delete [] re_i;
          delete [] im_i;
          delete [] mag_i;
         


/*
         TRandom1* myrand = new TRandom1();

         TH1F* myhist = new TH1F("stats","",100,0,10);

         for(Int_t i=0;i<10000;++i) {

             myhist->Fill(myrand->Gaus(5,1));

          }
          //f1->SetParameter(0,1);
          f1->SetParameter(1,4);
          f1->SetParameter(2,2);
          
          myhist->Fit(f1,"R");
          
          myhist->Draw(); 
          f1->Draw("same"); */ //esto no se usa es solo para una prueba

        gROOT->SetStyle("Plain");
        //gStyle->SetOptStat(000000000);//(111111111);
        gStyle->SetPalette(60);
        gStyle->SetOptFit(1);


        outf->Write();

         // rootapp->Run();
        return 0 ;
}//1 /* main */

