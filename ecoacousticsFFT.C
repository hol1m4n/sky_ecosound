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
        Int_t ref_power=-3;
        Int_t initial_Month=2;
        Int_t final_Month=2; //change these values according to available data 
        Int_t N=1024;//It has been choosed in the following way:  delta_t * delta_f > 1, then we choose the minimum frequency resolution as 50Hz. This implies that delta_t>1/50 and then N>1/50*44100 where 44100 is the sample rate. This gives N= 882 aprox and we decided N=1024
        Double_t scale_y=1; //must be sqrt(N);
        Double_t scale_x=1; //must be ((Double_t)(44100))/((Double_t)(N));
        Int_t N_buffers_for_test=100;// number of buffers (per file) used to testing
        Int_t N_files_for_test=20; //number of files used to testing

        bool test_buffers=true;//true; //define it as false if you don't want to run in test mode for a maximum of buffers per file
        bool test_files=true;//true; //define it as false if you don't want to run in test mode for a maximum of files
        
              
        gROOT->SetStyle("Plain");
        gStyle->SetOptStat(111111111);
        gStyle->SetPalette(60);
        gStyle->SetOptFit(1);

        TApplication* rootapp = new TApplication("example",&argc, argv);
        TFile * outf = new TFile("Plots.root", "RECREATE");
        char         *infilename;

        ofstream test_output_Mag_h;
        test_output_Mag_h.open("MagnitudePerHour.txt");

        ofstream test_output_Mag_m;
        test_output_Mag_m.open("MagnitudePerMonth.txt");

        ofstream test_output;
        test_output.open("test.txt");
 
        ofstream test_output1;
        test_output1.open("test1.txt");

        ofstream output;
        output.open("EntropyAndRMS.txt");

        float buf[N];
        Int_t readcount;
   
        string filename[11][2];
        filename[0][0]= "74";
        filename[0][1]="100";
        filename[1][0]="101";
        filename[1][1]="126";
        filename[2][0]="127";
        filename[2][1]="151";         
        filename[3][0]="154";
        filename[3][1]="179";
        filename[4][0]="180";
        filename[4][1]="203";
        filename[5][0]="204";
        filename[5][1]="228";
        filename[6][0]="229";
        filename[6][1]="247";
        filename[7][0]="248";
        filename[7][1]="268";
        filename[8][0]="270";
        filename[8][1]="300";
        filename[9][0]="301";
        filename[9][1]="319";
        filename[10][0]="320";
        filename[10][1]="339";

        Int_t start_time[18][2];
        start_time[0][0]=74;
        start_time[0][1]=6; 
        start_time[1][0]=101;
        start_time[1][1]=6; 
        start_time[2][0]=127;
        start_time[2][1]=18;
        start_time[3][0]=154;
        start_time[3][1]=7;
        start_time[4][0]=176;
        start_time[4][1]=6;
        start_time[5][0]=180;
        start_time[5][1]=18;
        start_time[6][0]=189;
        start_time[6][1]=6; 
        start_time[7][0]=204;
        start_time[7][1]=19; 
        start_time[8][0]=211;
        start_time[8][1]=6; 
        start_time[9][0]=219;
        start_time[9][1]=16;
        start_time[10][0]=229;
        start_time[10][1]=18; 
        start_time[11][0]=238;
        start_time[11][1]=6;
        start_time[12][0]=248;
        start_time[12][1]=18; 
        start_time[13][0]=259;
        start_time[13][1]=6;
        start_time[14][0]=270;
        start_time[14][1]=18; 
        start_time[15][0]=280;
        start_time[15][1]=6; 
        start_time[16][0]=301;
        start_time[16][1]=6;
        start_time[17][0]=320;
        start_time[17][1]=6; 

        //allocating variables to store results
        Double_t Magnitude_Spectrum_Hour[11][24][N/2+1];
        Double_t Magnitude_Spectrum_Month[11][N/2+1];
        Double_t Magnitude_Spectrum_Year[N/2+1];

        Double_t Standard_Deviation_Of_Magnitude_Spectrum_Hour[11][24][N/2+1];
        Double_t Positive_Mean_Of_Magnitude_Spectrum_Hour[11][24][N/2+1];//mean of the signals that are above the Magnitude_Spectrum
        Double_t Negative_Mean_Of_Magnitude_Spectrum_Hour[11][24][N/2+1];//mean of the signals that are below the Magnitude_Spectrum

        Double_t Standard_Deviation_Of_Magnitude_Spectrum_Month[11][N/2+1];
        Double_t Positive_Mean_Of_Magnitude_Spectrum_Month[11][N/2+1];
        Double_t Negative_Mean_Of_Magnitude_Spectrum_Month[11][N/2+1];

        Double_t Standard_Deviation_Of_Magnitude_Spectrum_Year[N/2+1];
        Double_t Positive_Mean_Of_Magnitude_Spectrum_Year[N/2+1];
        Double_t Negative_Mean_Of_Magnitude_Spectrum_Year[N/2+1];


        Double_t tmp; //dummy variable


        //Initial values in arrays
        for (Int_t fr=0; fr<N/2+1; fr++)
        {
           for (Int_t m=0; m<11; m++)
           {
              for (Int_t h=0; h<24; h++)
              {
                Magnitude_Spectrum_Hour[m][h][fr]=0;
                Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][fr]=0;
                Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr]=0;
                Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr]=0;
              }
              Magnitude_Spectrum_Month[m][fr]=0;
              Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr]=0;
              Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr]=0;
              Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr]=0;
           }  
           Magnitude_Spectrum_Year[fr]=0;
           Standard_Deviation_Of_Magnitude_Spectrum_Year[fr]=0;
           Positive_Mean_Of_Magnitude_Spectrum_Year[fr]=0;
           Negative_Mean_Of_Magnitude_Spectrum_Year[fr]=0;
        } 



        Int_t initial_file=0;
        Int_t final_file=0;

        Int_t next_change_index=0;

        Int_t previous_length=0;
        Int_t current_length=0;
        bool previous_length_is_normal=true;
        bool current_length_is_normal=true;
  
        Int_t cycle=0;
        Int_t h=0;

        Int_t number_of_samples_h[11][24];
        Int_t number_of_samples_m[11];
        Int_t number_of_samples_y=0;
        Int_t number_of_positive_samples_h[11][24][N/2+1];
        Int_t number_of_negative_samples_h[11][24][N/2+1];
        Int_t number_of_positive_samples_m[11][N/2+1];
        Int_t number_of_negative_samples_m[11][N/2+1];
        Int_t number_of_positive_samples_y[N/2+1];
        Int_t number_of_negative_samples_y[N/2+1];
        
        for (Int_t fr=0; fr<N/2+1;fr++)
        {
          number_of_positive_samples_y[fr]=0;
          number_of_negative_samples_y[fr]=0;
          for (Int_t m=0; m <11; m++)    
          {
            number_of_positive_samples_m[m][fr]=0;
            number_of_negative_samples_m[m][fr]=0;
            for (Int_t h=0; h<24; h++)
            {
              number_of_positive_samples_h[m][h][fr]=0;
              number_of_negative_samples_h[m][h][fr]=0;
            }
          }
        }
        for (Int_t m=0; m <11; m++)    
        {
          number_of_samples_m[m]=0;
          for (Int_t h=0; h<24; h++)
          {
             number_of_samples_h[m][h]=0;
          }
        }
       
        Int_t frames_in_this_hour=0; 
        for (Int_t proc=1; proc<=2; proc++)//in order to first calculate the mean value, it is processed the information and then it is done again but to calculate the positive and negative error.
        {//proc  
        for (Int_t m=initial_Month; m<=final_Month; m++) 
        {//2 
          initial_file = atoi((filename[m-2][0]).c_str());
          final_file = atoi((filename[m-2][1]).c_str());
          Int_t files_analyzed=0;
          bool condition_test_over_files;
          if (!test_files) condition_test_over_files=true;
            else condition_test_over_files= (files_analyzed<N_files_for_test);   
            
          for (Int_t f=initial_file;((f<=final_file) && condition_test_over_files);f++) 
          {//3
            SNDFILE    *infile = NULL;///
            SF_INFO    sfinfo;///
            memset (&sfinfo, 0, sizeof (sfinfo)) ;///
            Double_t *in = new Double_t[2*(N/2+1)];
      
            Double_t *re_full_tmp_1 = new Double_t[N/2+1];
            Double_t *im_full_tmp_1 = new Double_t[N/2+1];

             
            Double_t *mag_full_tmp = new Double_t[N/2+1];
            Double_t *pha_full_tmp = new Double_t[N/2+1];

            if (f==113 || f==258) f++; //the files 113.wav y 258.wav don't exist
            stringstream ss;
            ss<<f;
            string str=ss.str();
            string s_infilename = str+".wav";
            infilename = &s_infilename[0u];
            //begin: reviewing that run parameters are OK
            if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
            {    
               printf ("Not able to open input file %s.\n", infilename) ;
               puts (sf_strerror (NULL)) ;
               return 1 ;
            }
            //end: reviewing that run parameters are OK
            cout<<s_infilename<<endl;
            if (test_files) cout<<s_infilename<<endl;
            if (test_files) test_output1<<s_infilename<<endl;
            files_analyzed++;
            if (test_files) cout<< start_time[next_change_index][0]<<" "<<f<<endl;
            if ( start_time[next_change_index][0]==f  )
            {
               h =  start_time[next_change_index][1] - 1;
               if (test_files) cout<<"hora="<<h<<endl;
               frames_in_this_hour=0;
               next_change_index++;
               if (test_files) test_output1<<"entro"<<start_time[next_change_index][0]<<endl;
            }
            previous_length=current_length;
            previous_length_is_normal=true;
            if (previous_length<134152200)
            {
              previous_length_is_normal=false;
            }

            current_length=sfinfo.frames;
            if (test_files) cout<<"sfinfo.frames="<<sfinfo.frames<<endl;
            current_length_is_normal=true;
            if (current_length<134152200)
            {
              current_length_is_normal=false;
              if (test_files) cout<<"the file "<<f<<" has a length of "<<current_length<<endl;
            }

    
            cycle=0;
            Int_t buffers_analyzed=0;    
            bool condition_test_over_buffers;
            if (!test_buffers) condition_test_over_buffers=true;
            else condition_test_over_buffers= (buffers_analyzed<N_buffers_for_test);   
            while ( ( (readcount = sf_readf_float(infile, buf, N) )>0 ) && (condition_test_over_buffers) )
            {//4
              cycle = cycle + readcount;    
              frames_in_this_hour=frames_in_this_hour + readcount;
              if (frames_in_this_hour>=158760000)//frames in one hour (60*60*44100=158760000)
              {
                 if (h<24) 
                 { 
                    h++;
                 }
                 if (h==24)
                 {
                    h=0;
                 }
                 frames_in_this_hour=0;
                 if (test_files) cout<<"hora="<<h<<endl;
              }
            
              if ( (previous_length_is_normal || cycle>13230000)  && (current_length_is_normal || cycle<current_length-13230000) ) // the first 5 minutes (300*44100=13230000) are discarded if the file is the first one or if it follows one with length less than usual (50*42*44100 = 134152200). In the same way the last 5 minutes are discarded if the file doesn't have the usual length.
              {//5 
                // cout<<"entro y cargó"<<endl;
                 
                 
                 buffers_analyzed=buffers_analyzed + readcount; 

                 
                 //cout<<"frames_in_this_hour="<<frames_in_this_hour<<endl;
                 
        
                 if (test_buffers) cout<<"entro1"<<endl;
                 //cout<<"entro1"<<endl;
                 if (readcount == N)  //In order to only include streams with the same lenght,  because maybe the last one has a different one.
                 {//6    
                   ///begin:transform
                   //cout<<"entro2"<<endl;
                   for (Int_t i=0; i<N; i++)
                   {
                     in[i] =  convert_dB_to_Intensity(buf[i],ref_power);
                   }
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
                   fft_own->SetPoints(in);
                   fft_own->Transform();
                   fft_own->GetPointsComplex(re_full_tmp_1,im_full_tmp_1);
                   delete  fft_own;
                   for (Int_t i=0; i<N/2+1;i++)
                   {
                     mag_full_tmp[i]=sqrt(pow(re_full_tmp_1[i],2)+pow(im_full_tmp_1[i],2));
                     pha_full_tmp[i]=atan(im_full_tmp_1[i]/re_full_tmp_1[i]);
                   }
                   if (test_buffers) cout<<"entro2"<<endl;
                   ///end: transform
                   ///begin: store mean a error values per hour and month
                   
                   if (test_files) test_output1<<h<<endl;
                   if (proc==1) // if the information is processed for the first time, then the mean value is calculated
                   {
                   number_of_samples_h[m-2][h] = number_of_samples_h[m-2][h] + 1;
                   for (Int_t i=0; i<N/2+1; i++)
                   {
                     Magnitude_Spectrum_Hour[m-2][h][i] = Magnitude_Spectrum_Hour[m-2][h][i]+mag_full_tmp[i];
                    
                     if (test_buffers) cout<<"entro3"<<endl;

                   }
                   }
                   if (proc==2) // if the information is processed for second time, then the positive and negative errors are calculated.
                   {
                   for (Int_t i=0; i<N/2+1; i++)
                   {
                     tmp = mag_full_tmp[i] - Magnitude_Spectrum_Hour[m-2][h][i]/((Double_t)(number_of_samples_h[m-2][h]));
                     //cout<<"tmp="<<tmp<<endl;
                     Standard_Deviation_Of_Magnitude_Spectrum_Hour[m-2][h][i] = Standard_Deviation_Of_Magnitude_Spectrum_Hour[m-2][h][i] + pow(tmp,2);//Note:In fact, this is not the standard deviation but the squared standard deviation. This is fixed below.
                     if (tmp>0)
                     {
                       number_of_positive_samples_h[m-2][h][i] = number_of_positive_samples_h[m-2][h][i] + 1;
                       Positive_Mean_Of_Magnitude_Spectrum_Hour[m-2][h][i] = Positive_Mean_Of_Magnitude_Spectrum_Hour[m-2][h][i] + mag_full_tmp[i];  
                     }
                     else
                     {
                       number_of_negative_samples_h[m-2][h][i] = number_of_negative_samples_h[m-2][h][i] + 1;
                       Negative_Mean_Of_Magnitude_Spectrum_Hour[m-2][h][i] = Negative_Mean_Of_Magnitude_Spectrum_Hour[m-2][h][i] + mag_full_tmp[i];  
                     }
                   }
                   }
                   ///end: store mean a error values per hour and month
                 }//6
              }//5
              //cout<<"entro3"<<endl;
              if (!test_buffers) condition_test_over_buffers=true;
              else condition_test_over_buffers= (buffers_analyzed<N_buffers_for_test);   
            
            }//4
            if (!test_files) condition_test_over_files=true;
            else condition_test_over_files= (files_analyzed<N_files_for_test);   
            sf_close (infile);
            delete [] in;
      
            delete [] re_full_tmp_1;
            delete [] im_full_tmp_1;

             
            delete [] mag_full_tmp;
            delete [] pha_full_tmp;
            
            

          }//3
        }//2
        
        }//proc


        //Calculating data for months and taking the square root of Errors for hours
        for (Int_t fr=0; fr<N/2+1; fr++)
        {
           for (Int_t m=initial_Month-2; m<=final_Month-2; m++)
           {
              for (h=0; h<24; h++)
              {
                Magnitude_Spectrum_Month[m][fr]=Magnitude_Spectrum_Month[m][fr] +  Magnitude_Spectrum_Hour[m][h][fr];
                //cout<<"M="<<Magnitude_Spectrum_Hour[m][h][fr]<<endl; 
                Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr]=Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr] +  Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][fr];//Note: As earlier was mentioned it is the squared standard deviation, it is fixed below.
                Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr]=Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr] +  Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr];
                //cout<<"PM="<<Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr]<<endl;
                Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr]=Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr] +  Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr];
                //cout<<"NM="<<Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr]<<endl;

                //Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][fr]=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][fr]);//Now it is the standard deviation!!!
                //cout<<"SD="<<Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][fr]<<endl; 
                if (fr==0) number_of_samples_m[m]=number_of_samples_m[m]+number_of_samples_h[m][h]; //the total number of samples doesn't depend on fr
                number_of_positive_samples_m[m][fr]=number_of_positive_samples_m[m][fr]+number_of_positive_samples_h[m][h][fr];
                number_of_negative_samples_m[m][fr]=number_of_negative_samples_m[m][fr]+number_of_negative_samples_h[m][h][fr];         
              }
              

              Magnitude_Spectrum_Year[fr]=Magnitude_Spectrum_Year[fr] +  Magnitude_Spectrum_Month[m][fr];
              Standard_Deviation_Of_Magnitude_Spectrum_Year[fr]=Standard_Deviation_Of_Magnitude_Spectrum_Year[fr] +  Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr];//Note: As earlier was mentioned it is the squared standard deviation, it is fixed below.
              Positive_Mean_Of_Magnitude_Spectrum_Year[fr]=Positive_Mean_Of_Magnitude_Spectrum_Year[fr] +  Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr];
              Negative_Mean_Of_Magnitude_Spectrum_Year[fr]=Negative_Mean_Of_Magnitude_Spectrum_Year[fr] +  Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr];
              //cout<<"NM="<<Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][fr]<<endl;





              //Magnitude_Spectrum_Month[m][fr]=Magnitude_Spectrum_Month[m][fr]/((Double_t)24); //24 hours
              //Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr]=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr]);//Now it is the standard deviation!!!
              //cout<<"SD="<<Standard_Deviation_Of_Magnitude_Spectrum_Month[m][fr]<<endl;   
              //Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr]=Positive_Mean_Of_Magnitude_Spectrum_Month[m][fr]/((Double_t)24); //24 hours
              //Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr]=Negative_Mean_Of_Magnitude_Spectrum_Month[m][fr]/((Double_t)24); //24 hours            

              if (fr==0) number_of_samples_y=number_of_samples_y+number_of_samples_m[m];//the total number of samples doesn't depend on fr
              number_of_positive_samples_y[fr]=number_of_positive_samples_y[fr]+number_of_positive_samples_m[m][fr];
              number_of_negative_samples_y[fr]=number_of_negative_samples_y[fr]+number_of_negative_samples_m[m][fr];
         
           }
           
           
           //Magnitude_Spectrum_Year[fr]=Magnitude_Spectrum_Year[fr]/((Double_t)(24*(final_Month-initial_Month+1))); //24 hours and number of months analyzed
           //Standard_Deviation_Of_Magnitude_Spectrum_Year[fr]=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Year[fr]);//Now it is the standard deviation!!!
           //cout<<"SD="<<Standard_Deviation_Of_Magnitude_Spectrum_Year[fr]<<endl;   
           //Positive_Mean_Of_Magnitude_Spectrum_Year[fr]=Positive_Mean_Of_Magnitude_Spectrum_Year[fr]/((Double_t)(24*(final_Month-initial_Month+1))); //24 hours and number of months analyzed
           //Negative_Mean_Of_Magnitude_Spectrum_Year[fr]=Negative_Mean_Of_Magnitude_Spectrum_Year[fr]/((Double_t)(24*(final_Month-initial_Month+1))); //24 hours and number of months analyzed        
         
        }
        if (test_buffers) 
        {
           cout<<"entro4"<<endl;
           for (Int_t i=0;i<N/2+1;i++)
           {   
             test_output<<Magnitude_Spectrum_Hour[initial_Month-2][6][i]<<endl;
           }
        }


////////////////Nota:  si se desea hacer un fit descomentar las que están comentadas como ////
        cout<<"entrando a histos"<<endl;
        TH1D * h_Magnitude_h[11][24];
        TH1D * h_Magnitude_m[11];
        TH1D * h_Magnitude_y;
        TH1D * h_SD_Magnitude_h[11][24];
        TH1D * h_SD_Magnitude_m[11];
        TH1D * h_SD_Magnitude_y;
        TH1D * h_Positive_Mean_Magnitude_h[11][24];
        TH1D * h_Positive_Mean_Magnitude_m[11];
        TH1D * h_Positive_Mean_Magnitude_y;
        TH1D * h_Negative_Mean_Magnitude_h[11][24];
        TH1D * h_Negative_Mean_Magnitude_m[11];
        TH1D * h_Negative_Mean_Magnitude_y;
        TH1D * h_Discriminant_h[11][24];
        TH1D * h_Discriminant_m[11];
        TH1D * h_Discriminant_y;
 
       


        Double_t Entropy_h[11][24];
        Double_t Entropy_m[11];

        Double_t RMS_h[11][24];
        Double_t RMS_m[11];

        for (Int_t m=0; m<11; m++)
        {
          for (Int_t h=0; h<24; h++)
          {
            Entropy_h[m][h]=0;
            RMS_h[m][h]=0;
          }
          Entropy_m[m]=0;
          RMS_m[m]=0;
        }


        string month[11];
        month[0]="Feb";
        month[1]="Mar";
        month[2]="Apr";
        month[3]="May";
        month[4]="Jun";
        month[5]="Jul";
        month[6]="Aug";
        month[7]="Sep";
        month[8]="Oct";
        month[9]="Nov";
        month[10]="Dec";

        Int_t hour;

        Double_t X[N/2+1];
        Double_t errorX[N/2+1];
        Int_t nbins= N/2;
        ////TF1* f1 = new TF1("f1","gaus(0)+expo(3)+gaus(5)", 0, 400);//(N/2-1)*scale_x);//gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0. expo(3) is a substitute for exp([3]+[4]*x)
        //TF1* f1 = new TF1("f1",fitfunc,0, (N/2-1)*scale_x,3); In case a specific fit function is needed
        ////TF1* f2 = new TF1("f2","gaus(0)+expo(3)+gaus(5)",2000*((Double_t)(N))/((Double_t)(44100))*scale_x,8000*((Double_t)(N))/((Double_t)(44100))*scale_x);// by paper, the anthrophony is between 0.2 and 2 kHz, and biophony between 2 and 8 kHz
        ////f1->SetLineColor(2);
        ////f2->SetLineColor(3);
        ////Double_t parameters[3];


        for (Int_t i=0; i<N/2+1; i++)
          {
            X[i]=i*scale_x; //NOTE: for "real" frequencies you have to scale the x-axes range with the range of your function (in this case 1 sec is 44100 samples and thus N samples are N/44100) 
            errorX[i]=scale_x; //using uncertainty relation
          }
       
        if (test_buffers) cout<<"entro6"<<endl;

        for (Int_t m=initial_Month-2; m<=final_Month-2; m++)
        {
          for (h=0; h<24; h++)
          {
              
             hour= h+1;
             stringstream ss;
             ss<<hour;
             string str=ss.str();
             h_Magnitude_h[m][h] = new TH1D(("Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), ("Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), nbins, 0, (nbins-1)*scale_x);
             h_SD_Magnitude_h[m][h] = new TH1D(("Standard Deviation over Mean Value for Month = "+month[m]+", Hour = "+str).c_str(), ("Standard Deviation over Mean Value for Month = "+month[m]+", Hour = "+str).c_str(), nbins, 0, (nbins-1)*scale_x);
             h_Positive_Mean_Magnitude_h[m][h] = new TH1D(("Positive Mean of Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), ("Positive Mean of Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), nbins, 0, (nbins-1)*scale_x);
             h_Negative_Mean_Magnitude_h[m][h] = new TH1D(("Negative Mean of Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), ("Negative Mean of Spectrum's Magnitude Month = "+month[m]+", Hour = "+str).c_str(), nbins, 0, (nbins-1)*scale_x);
             h_Discriminant_h[m][h] = new TH1D(("Spectrum's Discriminant Month = "+month[m]+", Hour = "+str).c_str(), ("Spectrum's Discriminant Month = "+month[m]+", Hour = "+str).c_str(), nbins, 0, (nbins-1)*scale_x);
            
             Double_t normalization_factor=0;

             ofstream datah;
             datah.open(("Data/"+month[m]+"/Magnitudes_"+month[m]+"_"+str+".txt").c_str());
             datah<<"Media Media_Superior Media_Inferior Desviación_Estandar Discriminante"<<endl; 
         
             for (Int_t i=0; i<nbins+1; i++)
             {
                //h_Magnitude_h[m][h]->SetBinContent(i,convert_Intensity_to_dB(Magnitude_Spectrum_Hour[m][h][i]/((Double_t)scale_y),ref_power));//y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
                Double_t magnitude_h=0;
                Double_t sd_h=0;
                Double_t positive_magnitude_h=0;
                Double_t negative_magnitude_h=0;
                Double_t discriminant_h=0;
                if (number_of_samples_h[m][h]!=0) 
                {
                  magnitude_h=Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h]));
                  h_Magnitude_h[m][h]->SetBinContent(i,convert_Intensity_to_dB(magnitude_h,ref_power));
                }
                if (Magnitude_Spectrum_Hour[m][h][i]!=0)
                { 
                  sd_h=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][i])/((Double_t)(Magnitude_Spectrum_Hour[m][h][i]));
                  h_SD_Magnitude_h[m][h]->SetBinContent(i,sd_h);
                }


                if (number_of_positive_samples_h[m][h][i]!=0) 
                {
                  positive_magnitude_h=Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_positive_samples_h[m][h][i]));
                  h_Positive_Mean_Magnitude_h[m][h]->SetBinContent(i,convert_Intensity_to_dB(positive_magnitude_h,ref_power));
                  cout<<"MAYOR"<<endl;
                }


                if (number_of_negative_samples_h[m][h][i]!=0)
                {
                  negative_magnitude_h=Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_negative_samples_h[m][h][i]));
                  h_Negative_Mean_Magnitude_h[m][h]->SetBinContent(i,convert_Intensity_to_dB(negative_magnitude_h,ref_power));
                }

                if (number_of_negative_samples_h[m][h][i]*number_of_positive_samples_h[m][h][i]*number_of_samples_h[m][h]!=0)
                {
                  discriminant_h=(positive_magnitude_h-negative_magnitude_h)/( (Double_t) magnitude_h);
                  h_Discriminant_h[m][h]->SetBinContent(i,discriminant_h);
                }


                datah<<magnitude_h<<" "<<positive_magnitude_h<<" "<<negative_magnitude_h<<" "<<sd_h<<" "<<discriminant_h<<endl;


                if (test_files) cout<<"Plot_M:"<<m<<" "<<h<<" "<<Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h]))<<endl;

                if (test_files) cout<<"Plot_SD:"<<m<<" "<<h<<" "<<Standard_Deviation_Of_Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(Magnitude_Spectrum_Hour[m][h][i]*number_of_samples_h[m][h]))<<endl;

                if (test_files) cout<<"Plot_PM:"<<m<<" "<<h<<" "<<(Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][i] - Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])) * number_of_positive_samples_h[m][h][i])/((Double_t)(number_of_samples_h[m][h]))+Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h]))<<endl;

                if (test_files) cout<<"Plot_NM:"<<m<<" "<<h<<" "<<(Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][i] - Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])) * number_of_negative_samples_h[m][h][i])/((Double_t)(number_of_samples_h[m][h]))+Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h]))<<endl;

                if (test_files) cout<<"Plot_Discriminant:"<<m<<" "<<h<<" "<<((Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][i] - Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])) * number_of_positive_samples_h[m][h][i])/((Double_t)(number_of_samples_h[m][h]))+Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])) - (Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][i] - Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])) * number_of_negative_samples_h[m][h][i])/((Double_t)(number_of_samples_h[m][h]))+Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h])))/Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(number_of_samples_h[m][h]))<<endl;

                //cout<<"Positive - Negative: "<<Positive_Mean_Of_Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(scale_y))<<"-"<<Negative_Mean_Of_Magnitude_Spectrum_Hour[m][h][i]/((Double_t)(scale_y));
                 

                
                Entropy_h[m][h]=Entropy_h[m][h] -  Magnitude_Spectrum_Hour[m][h][i]/((Double_t)scale_y) *  log(Magnitude_Spectrum_Hour[m][h][i]/((Double_t)scale_y));// / nbins; quité el dividido porque no veo de donde sale, pero es mejor revisarlo bien, creo que en la ecuación de entropía la magnitud debe estar normalizada
                normalization_factor = normalization_factor + Magnitude_Spectrum_Hour[m][h][i]/((Double_t)scale_y);
             }
             
             datah.close();
             Entropy_h[m][h]=Entropy_h[m][h] / normalization_factor;

             ////TCanvas c1("c1","Magnitude",200,10,600,400); 
             ////f1->SetParameter(1,50);
             ////f1->SetParameter(2,150);
             ////f1->SetParameter(3,-5.5);
             ////f1->SetParameter(4,-0.1);
             ////f1->SetParameter(6,100);
             ////f1->SetParameter(7,5);
             
             ////parameters[0]=f1->GetParameter(0);
             ////parameters[1]=f1->GetParameter(1);
             ////parameters[2]=f1->GetParameter(2);
             ////if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             ////test_output_Mag_h<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             ////h_Magnitude_h[m][h]->Fit(f1,"R");
             
             TCanvas *ch = new TCanvas("ch","Spectrum's Magnitudes Per Hour",0,0,600,600);
             
             
             gStyle->SetOptStat(0);
             ch->Divide(1,3);
             ch->cd(1);
             h_Magnitude_h[m][h]->Draw();
             h_Positive_Mean_Magnitude_h[m][h]->SetLineColor(kBlue);
             h_Negative_Mean_Magnitude_h[m][h]->SetLineColor(kRed);
             h_Positive_Mean_Magnitude_h[m][h]->Draw("same");
             h_Negative_Mean_Magnitude_h[m][h]->Draw("same");
             ch->cd(2);
             h_SD_Magnitude_h[m][h]->Draw();
             ch->cd(3);
             h_Discriminant_h[m][h]->SetLineColor(kGreen);
             h_Discriminant_h[m][h]->Draw();
             ch->Update();
             //gSystem->mkdir("Plots");
             //gSystem->mkdir(("Plots/"+month[m]).c_str());
             ch->SaveAs(("Plots/"+month[m]+"/Magnitudes_"+month[m]+"_"+str+".png").c_str());
             ch->Clear();
             
             RMS_h[m][h] = h_Magnitude_h[m][h]->GetRMS();
             output<<Entropy_h[m][h]<<" "<<RMS_h[m][h]<<endl;
             ////f1->Draw("same");
             if (test_buffers) cout<<"entro7"<<endl;

             /*//same fit but filtering anthropogenic sources <2000Hz
             TCanvas c2("c2","Magnitude Biophony",200,10,600,400);
             f2->SetParameter(1,50);
             f2->SetParameter(2,150);
             f2->SetParameter(3,-5.5);
             f2->SetParameter(4,-0.1);
             f2->SetParameter(6,100);
             f2->SetParameter(7,5);
             h_Magnitude_h[m][h]->Fit(f2,"R");             
             parameters[0]=f2->GetParameter(0);
             parameters[1]=f2->GetParameter(1);
             parameters[2]=f2->GetParameter(2);
             if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             test_output_Mag_h<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             h_Magnitude_h[m][h]->Draw();
             f2->Draw("same");
             if (test_buffers) cout<<"entro8"<<endl;

             /*f1->SetParameter(0,1);
             f1->SetParameter(1,1);
             f1->SetParameter(2,1);
             grafica = new TGraphErrors(nbins, X, Phase_Spectrum_Hour[m][h], errorX, Error_In_Phase_Spectrum_Hour[m][h]);
             grafica->Fit(f1,"R");
             parameters[0]=f1->GetParameter(0);
             parameters[1]=f1->GetParameter(1);
             parameters[2]=f1->GetParameter(2);
             if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             test_output_Pha_h<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             //h_Phase_h[m][h]->GetXaxis()->SetRange(2000,(N/2-1)*scale_x);
             h_Phase_h[m][h]->Draw();
             f1->Draw("same");

             //same fit but filtering anthropogenic sources <2000Hz
             f2->SetParameter(0,1);
             f2->SetParameter(1,1);
             f2->SetParameter(2,1);
             grafica->Fit(f2,"R");
             parameters[0]=f2->GetParameter(0);
             parameters[1]=f2->GetParameter(1);
             parameters[2]=f2->GetParameter(2);
             if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             test_output_Pha_h<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
             f2->Draw("same");*/

          } 
          h_Magnitude_m[m] = new TH1D(("Spectrum's Magnitude Month = "+month[m]).c_str(), ("Spectrum's Magnitude Month = "+month[m]).c_str(), nbins, 0, (nbins-1)*scale_x);
          h_SD_Magnitude_m[m] = new TH1D(("Standard Deviation_over Mean Value for Month = "+month[m]).c_str(), ("Standard Deviation_over Mean Value for Month = "+month[m]).c_str(), nbins, 0, (nbins-1)*scale_x);
          h_Positive_Mean_Magnitude_m[m] = new TH1D(("Positive_Standard Deviation_Spectrum's Magnitude Month = "+month[m]).c_str(), ("Positive_Standard_Deviation_Spectrum's Magnitude Month = "+month[m]).c_str(), nbins, 0, (nbins-1)*scale_x);
          h_Negative_Mean_Magnitude_m[m] = new TH1D(("Negative_Standard Deviation_Spectrum's Magnitude Month = "+month[m]).c_str(), ("Negative_Standard_Deviation_Spectrum's Magnitude Month = "+month[m]).c_str(), nbins, 0, (nbins-1)*scale_x);
          h_Discriminant_m[m] = new TH1D(("Spectrum's Discriminant Month = "+month[m]).c_str(), ("Spectrum's Discriminant Month = "+month[m]).c_str(), nbins, 0, (nbins-1)*scale_x);

          Double_t normalization_factor=0;
          ofstream datam;
          datam.open(("Data/"+month[m]+"/Magnitudes_"+month[m]+".txt").c_str());
          datam<<"Media Media_Superior Media_Inferior Desviación_Estandar Discriminante"<<endl;  
         
          for (Int_t i=0; i<nbins + 1; i++)
             {


         
                if (test_files) cout<<"Plot:"<<m<<" "<<Magnitude_Spectrum_Month[m][i]<<endl;
                if (test_files) cout<<"Plot_SD:"<<m<<" "<<Standard_Deviation_Of_Magnitude_Spectrum_Month[m][i]/((Double_t)(Magnitude_Spectrum_Month[m][i]))<<endl;
                if (test_files) cout<<"Plot_EP:"<<m<<" "<<(Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m]))<<endl;
                if (test_files) cout<<"Plot_EN:"<<m<<" "<<(Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m]))<<endl;
                if (test_files) cout<<"Plot_Discriminant:"<<m<<" "<<((Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) - (Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])))/Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m]))<<endl;


                Double_t magnitude_m=0;
                Double_t sd_m=0;
                Double_t positive_magnitude_m=0;
                Double_t negative_magnitude_m=0;
                Double_t discriminant_m=0;
                if (number_of_samples_m[m]!=0) 
                {
                  magnitude_m=Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m]));
                  h_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB(magnitude_m,ref_power));
                }
                if (Magnitude_Spectrum_Month[m][i]!=0)
                {
                  sd_m=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Month[m][i])/((Double_t)(Magnitude_Spectrum_Month[m][i]));
                  h_SD_Magnitude_m[m]->SetBinContent(i,sd_m);
                  cout<<"sd"<<sd_m<<endl;
                }


                if (number_of_positive_samples_m[m][i]!=0) 
                {
                  positive_magnitude_m=Positive_Mean_Of_Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_positive_samples_m[m][i]));
                  h_Positive_Mean_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB(positive_magnitude_m,ref_power));
                  cout<<"MAYOR_M"<<endl;
                }


                if (number_of_negative_samples_m[m][i]!=0)
                {
                  negative_magnitude_m=Negative_Mean_Of_Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_negative_samples_m[m][i]));
                  h_Negative_Mean_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB(negative_magnitude_m,ref_power));
                }

                if (number_of_negative_samples_m[m][i]*number_of_positive_samples_m[m][i]*number_of_samples_m[m]!=0)
                {
                  discriminant_m=(positive_magnitude_m-negative_magnitude_m)/( (Double_t) magnitude_m);
                  h_Discriminant_m[m]->SetBinContent(i,discriminant_m);
                }


                datam<<magnitude_m<<" "<<positive_magnitude_m<<" "<<negative_magnitude_m<<" "<<sd_m<<" "<<discriminant_m<<endl;



                /*if (number_of_samples_m[m]!=0)
                {
                   h_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB(Magnitude_Spectrum_Month[m][i],ref_power));         
                   h_SD_Magnitude_m[m]->SetBinContent(i,Standard_Deviation_Of_Magnitude_Spectrum_Month[m][i]/((Double_t)(Magnitude_Spectrum_Month[m][i])));
                }
                if (number_of_positive_samples_m[m][i]!=0)  h_Positive_Mean_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB((Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])),ref_power));
                if (number_of_negative_samples_m[m][i]!=0)  h_Negative_Mean_Magnitude_m[m]->SetBinContent(i,convert_Intensity_to_dB((Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])),ref_power));
                if (number_of_negative_samples_m[m][i]*number_of_positive_samples_m[m][i]*number_of_samples_m[m]!=0) h_Discriminant_m[m]->SetBinContent(i,((Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) - (Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])))/Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])));*/

                Entropy_m[m]=Entropy_m[m] -  Magnitude_Spectrum_Month[m][i]/((Double_t)scale_y) *  log(Magnitude_Spectrum_Month[m][i]/((Double_t)scale_y)) / nbins;
                normalization_factor = normalization_factor + Magnitude_Spectrum_Month[m][i]/((Double_t)scale_y);

                //datam<<convert_Intensity_to_dB((Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])),ref_power)<<" "<<convert_Intensity_to_dB(Magnitude_Spectrum_Month[m][i],ref_power)<<" "<<convert_Intensity_to_dB((Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])),ref_power)<<" "<<Standard_Deviation_Of_Magnitude_Spectrum_Month[m][i]/((Double_t)(Magnitude_Spectrum_Month[m][i]))<<" "<<((Positive_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_positive_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) - (Negative_Mean_Of_Magnitude_Spectrum_Month[m][i] - Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])) * number_of_negative_samples_m[m][i])/((Double_t)(number_of_samples_m[m]))+Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m])))/Magnitude_Spectrum_Month[m][i]/((Double_t)(number_of_samples_m[m]))<<endl;
             }
          datam.close();
          Entropy_m[m]=Entropy_m[m] / normalization_factor;
          if (test_buffers) 
          {
            cout<<"entro4"<<endl;
            test_output<<"entro28"<<endl;
            for (Int_t i=0;i<N/2+1;i++)
            {    
              test_output<<Magnitude_Spectrum_Hour[initial_Month-2][6][i]<<endl;
            }
          }
          ////TCanvas c3("c3","Magnitude Month",200,10,600,400);
         
          ////f1->SetParameter(1,50);
          ////f1->SetParameter(2,150);
          ////f1->SetParameter(3,-5.5);
          ////f1->SetParameter(4,-0.1);
          ////f1->SetParameter(6,100);
          ////f1->SetParameter(7,5);

          ////if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          ////test_output_Mag_m<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          
          ////h_Magnitude_m[m]->Fit(f1,"R");

          TCanvas *cm = new TCanvas("ch","Spectrum's Magnitudes Per Month",0,0,600,600);
          gStyle->SetOptStat(0);
          cm->Divide(1,3);
          cm->cd(1);
          h_Magnitude_m[m]->Draw();
          h_Positive_Mean_Magnitude_m[m]->SetLineColor(kBlue);
          h_Negative_Mean_Magnitude_m[m]->SetLineColor(kRed);
          h_Positive_Mean_Magnitude_m[m]->Draw("same");
          h_Negative_Mean_Magnitude_m[m]->Draw("same");
          cm->cd(2);
          h_SD_Magnitude_m[m]->Draw();
          cm->cd(3);          
          h_Discriminant_m[m]->SetLineColor(kGreen);
          h_Discriminant_m[m]->Draw("same");
          cm->Update();
          //gSystem->mkdir("Plots");
          //gSystem->mkdir(("Plots/"+month[m]).c_str());
          cm->SaveAs(("Plots/"+month[m]+"/Magnitudes_"+month[m]+".png").c_str());
          cm->Clear();
 
          RMS_m[m] = h_Magnitude_m[m]->GetRMS();
          output<<Entropy_m[m]<<" "<<RMS_m[m]<<endl;
          ////f1->Draw("same");
          ////c3.SaveAs("prueba.png");
          
          /*TCanvas c4("c4","Magnitude Month - Biophony",200,10,600,400);
          f2->SetParameter(1,50);
          f2->SetParameter(2,150);
          f2->SetParameter(3,-5.5);
          f2->SetParameter(4,-0.1);
          f2->SetParameter(6,100);
          f2->SetParameter(7,5);
          h_Magnitude_m[m]->Fit(f2,"R");
          h_Magnitude_m[m]->Draw();
          f2->Draw("same");
          c4.SaveAs("prueba1.png");
         

          if (test_buffers) cout<<"entro9"<<endl;

          //same fit but filtering anthropogenic sources <2000Hz  //comentado para pruebas
          /*f2->SetParameter(0,1);
          f2->SetParameter(1,1);
          f2->SetParameter(2,1);
          grafica->Fit("f2","R");
          parameters[0]=f2->GetParameter(0);
          parameters[1]=f2->GetParameter(1);
          parameters[2]=f2->GetParameter(2);
          if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          test_output_Mag_m<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          f2->Draw("same");


          /*f1->SetParameter(0,1);
          f1->SetParameter(1,1);
          f1->SetParameter(2,1);
          grafica = new TGraphErrors(nbins, X, Phase_Spectrum_Month[m], errorX, Error_In_Phase_Spectrum_Month[m]);
          grafica->Fit(f1,"R");
          parameters[0]=f1->GetParameter(0);
          parameters[1]=f1->GetParameter(1);
          parameters[2]=f1->GetParameter(2);
          if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          test_output_Pha_m<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          //h_Phase_m[m]->GetXaxis()->SetRange(2000,(N/2-1)*scale_x);
          h_Phase_m[m]->Draw();
          f1->Draw("same");

          //same fit but filtering anthropogenic sources <2000Hz
          f2->SetParameter(0,1);
          f2->SetParameter(1,1);
          f2->SetParameter(2,1);
          grafica->Fit(f2,"R");
          parameters[0]=f2->GetParameter(0);
          parameters[1]=f2->GetParameter(1);
          parameters[2]=f2->GetParameter(2);
          if (test_files) cout<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          test_output_Pha_m<<m<<" "<<h<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<endl;
          f2->Draw("same");*/

           
        }
        h_Magnitude_y = new TH1D("Spectrum's Magnitude Year", "Spectrum's Magnitude Year", nbins, 0, (nbins-1)*scale_x);
        h_SD_Magnitude_y = new TH1D("Standard Deviation over Mean Value for the whole year", "Standard Deviation over Mean Value for the whole year", nbins, 0, (nbins-1)*scale_x);
        h_Positive_Mean_Magnitude_y = new TH1D("Positive_Standard Deviation_Spectrum's Magnitude Year", "Positive_Standard_Deviation_Spectrum's Magnitude Year", nbins, 0, (nbins-1)*scale_x);
        h_Negative_Mean_Magnitude_y = new TH1D("Negative_Standard Deviation_Spectrum's Magnitude Year", "Negative_Standard_Deviation_Spectrum's Magnitude Year", nbins, 0, (nbins-1)*scale_x);
        h_Discriminant_y = new TH1D("Spectrum's Discriminant Year", "Spectrum's Discriminant Year", nbins, 0, (nbins-1)*scale_x);
        
        Double_t discriminant[N/2+1];
        ofstream datay;
        datay.open("Data/Magnitudes.txt");
        datay<<"Media Media_Superior Media_Inferior Desviación_Estandar Discriminante"<<endl; 
        
          
        for (Int_t i=0; i<nbins + 1; i++)
        {

          //      if (number_of_samples_y!=0)  h_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB(Magnitude_Spectrum_Year[i],ref_power));
                if (test_files) cout<<"Plot:"<<" "<<Magnitude_Spectrum_Year[i]<<endl;
                if (test_files) cout<<"Plot_E:"<<" "<<Standard_Deviation_Of_Magnitude_Spectrum_Year[i]/((Double_t)(Magnitude_Spectrum_Year[i]))<<endl;
                if (test_files) cout<<"Plot_EP:"<<" "<<(Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y))<<endl;
                if (test_files) cout<<"Plot_EN:"<<" "<<(Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y))<<endl;
                if (test_files) cout<<"Plot_Discriminant: "<<((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) - (Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)))/Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y))<<endl;
     


                Double_t magnitude_y=0;
                Double_t sd_y=0;
                Double_t positive_magnitude_y=0;
                Double_t negative_magnitude_y=0;
                Double_t discriminant_y=0;
                if (number_of_samples_y!=0) 
                {
                  magnitude_y=Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y));
                  h_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB(magnitude_y,ref_power));
                }
                if (Magnitude_Spectrum_Year[i]!=0)
                {
                  sd_y=sqrt(Standard_Deviation_Of_Magnitude_Spectrum_Year[i])/((Double_t)(Magnitude_Spectrum_Year[i]));
                  h_SD_Magnitude_y->SetBinContent(i,sd_y);
                }


                if (number_of_positive_samples_y[i]!=0) 
                {
                  positive_magnitude_y=Positive_Mean_Of_Magnitude_Spectrum_Year[i]/((Double_t)(number_of_positive_samples_y[i]));
                  h_Positive_Mean_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB(positive_magnitude_y,ref_power));
                  cout<<"MAyor_Y"<<endl;
                }


                if (number_of_negative_samples_y[i]!=0)
                {
                  negative_magnitude_y=Negative_Mean_Of_Magnitude_Spectrum_Year[i]/((Double_t)(number_of_negative_samples_y[i]));
                  h_Negative_Mean_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB(negative_magnitude_y,ref_power));
                }

                if (number_of_negative_samples_y[i]*number_of_positive_samples_y[i]*number_of_samples_y!=0)
                {
                  discriminant_y=(positive_magnitude_y-negative_magnitude_y)/( (Double_t) magnitude_y);
                  h_Discriminant_y->SetBinContent(i,discriminant_y);
                  discriminant[i] = discriminant_y;
                }


                datay<<magnitude_y<<" "<<positive_magnitude_y<<" "<<negative_magnitude_y<<" "<<sd_y<<" "<<discriminant_y<<endl;

                
            /*    if (number_of_samples_y!=0) 
                {
                  h_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB(Magnitude_Spectrum_Year[i],ref_power));
                  h_SD_Magnitude_y->SetBinContent(i,Standard_Deviation_Of_Magnitude_Spectrum_Year[i]/((Double_t)(Magnitude_Spectrum_Year[i])));
                }

                if (number_of_positive_samples_y[i]!=0)  h_Positive_Mean_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)),ref_power));

                if (number_of_negative_samples_y[i]!=0)  h_Negative_Mean_Magnitude_y->SetBinContent(i,convert_Intensity_to_dB((Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)),ref_power));

                if (number_of_negative_samples_y[i]*number_of_positive_samples_y[i]*number_of_samples_y!=0)  h_Discriminant_y->SetBinContent(i,((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) - (Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)))/Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)));

                discriminant[i] = ((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) - (Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)))/Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y));
                if (test_files) cout<<"Discriminant per year: "<<discriminant[i]<<endl;
                datay<<convert_Intensity_to_dB((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)),ref_power)<<" "<<convert_Intensity_to_dB(Magnitude_Spectrum_Year[i],ref_power)<<" "<<convert_Intensity_to_dB((Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)),ref_power)<<" "<<Standard_Deviation_Of_Magnitude_Spectrum_Year[i]/((Double_t)(Magnitude_Spectrum_Year[i]))<<" "<<((Positive_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_positive_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) - (Negative_Mean_Of_Magnitude_Spectrum_Year[i] - Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)) * number_of_negative_samples_y[i])/((Double_t)(number_of_samples_y))+Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y)))/Magnitude_Spectrum_Year[i]/((Double_t)(number_of_samples_y))<<endl;*/
        
        }
        datay.close();
       
        TCanvas *cy = new TCanvas("ch","Spectrum's Magnitudes Per Year",0,0,600,600);
        gStyle->SetOptStat(0);
        cy->Divide(1,3);
        cy->cd(1);
        h_Magnitude_y->Draw();
        h_Positive_Mean_Magnitude_y->SetLineColor(kBlue);
        h_Negative_Mean_Magnitude_y->SetLineColor(kRed);
        h_Positive_Mean_Magnitude_y->Draw("same");
        h_Negative_Mean_Magnitude_y->Draw("same");
        cy->cd(2);
        h_SD_Magnitude_y->Draw();
        cy->cd(3);
        h_Discriminant_y->SetLineColor(kGreen);
        h_Discriminant_y->Draw();
        cy->Update();
        //gSystem->mkdir("Plots");
        //gSystem->mkdir(("Plots/"+month[m]).c_str());
        cy->SaveAs("Plots/Magnitudes.png");
        cy->Clear();

        TCanvas *c1 = new TCanvas("c1","Discriminant",200,10,200+nbins,600); 
        TGraph *g_Discriminant_y = new TGraph (nbins, X, discriminant);
        g_Discriminant_y->SetLineColor(8);
        g_Discriminant_y->SetLineWidth(3);
        g_Discriminant_y->SetTitle("Discriminant vs Frequency");
        g_Discriminant_y->GetXaxis()->SetTitle("Frecuency(Hz)");
        g_Discriminant_y->GetYaxis()->SetTitle("Discriminant");
        g_Discriminant_y->SetMinimum(-3);
        g_Discriminant_y->SetMaximum(3);
    
        g_Discriminant_y->Draw("");
        c1->SaveAs("Discriminant.png");
        if (test_files) cout<<"Saving Discriminant as Discriminant.png"<<endl;
        c1->Clear();
    

     test_output_Mag_h.close();
     test_output_Mag_m.close();
     test_output.close();
     test_output1.close();
     output.close();

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



         outf->Write();

         // rootapp->Run();
         return 0 ;
}//1 /* main */

