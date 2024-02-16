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


//****************************************************************************************************************************************************
//REQUIRED FUNCTIONS FOR THE CALCULATIONS*************************************************************************************************************
//****************************************************************************************************************************************************


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double deg_torad(double degree)// Function that converts Degrees to Radians
{
	double pi=TMath::Pi();
	double rad= degree*(pi/180);
	return rad;
}

double rad_todeg(double radians)// Function that converts Radians to Degrees 
{
	double pi=TMath::Pi();
	double deg= radians*(180/pi);
	return deg;
}

double algo(double value)   // Algorithm for subtract multiples of 360°
{
    // Start
	int al_one;
	al_one=value/360;
    double al_two[2];
    al_two[0]=(value/360)-al_one;
	al_two[1]=360*al_two[0];
	// Finish
	if(value>360)
	{
	    return al_two[1];
	}
	else
	{
		return value;
	}
	
}

double sp_algo(double value)   // Special algorithm for subtract multiples of 360° for Moon calculation cycles
{
    // Start
	int al_one;
	al_one=value/360;
    double al_two[2];
    al_two[0]=(value/360)-al_one;
	al_two[1]=360*al_two[0];
	// Finish
	if(value>360)
	{
	    return al_two[1];
	}
	if((value>=0)&&(value<=360))
	{
		return value;
	}
	if(value<0)
	{
        if((value<0)&&(value>=(-360)))
		{
			return value+360;
		}
		else
		{
		    al_one=(value/360)*(-1);
            al_two[0]=((value/360)+al_one);
	        al_two[1]=360*al_two[0];
			return al_two[1]+360;
		}
	}
	
	
	
	
	
}

double E_f(double e, double Ms) //Iteration in the calculation of the eccentric anomaly of the Sun in radians
{
	double M_r=deg_torad(Ms);
	double E[50]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double ED;
	double Lm;
	double E_iterada;
	E[0]= Ms+((180/TMath::Pi())*e*sin(M_r)*(1+(e*cos(M_r))));
	E[0]=deg_torad(E[0]);
	Lm=E[0]-(e*sin(E[0]))-M_r;
	ED=Lm/(1-(e*cos(E[0])));
	E_iterada=E[0]-ED;
	int i=0;
    while(Lm>0.000001)
	{
		E[i+1]=E[i]-ED;
	    E_iterada=E[i+1];
		Lm=E[i+1]-(e*sin(E[i+1]))-M_r;
		ED=Lm/(1-(e*cos(E[i+1])));
	    i++;
	}
	
	return E_iterada;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//****************************************************************************************************************************************************
//SUN CALCULATIONS CODE*******************************************************************************************************************************
//****************************************************************************************************************************************************


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double height(int a,int m,int d,double time)
{
	double lon=-73.370721; //Longitude
	double lat=5.865964;   //Latitude
	double JD=0; //Julian Date Number
	double r=0; //The distance of the Earth from the Sun
	double hX=0; //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	double hY=0; //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	double hZ=0; //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	double gX=0; //Rectangular Geocentric Ecliptical Coordinates. X Axis. 
	double gY=0; //Rectangular Geocentric Ecliptical Coordinates. Y Axis.
	double gZ=0; //Rectangular Geocentric Ecliptical Coordinates. Z Axis.
	double hei=0; //The height
    double aziS=0; //The Azimuth from the South
	double aziN=0; //The Azimuth from the North
	double M[2]={0,0}; // Mean Anomaly variable
	double n=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	float e=0.01671; //Eccentricity of the Earth's orbit, an ellipse. 
	double v[2]={0,0}; //True AnoM[1]aly of the Earth. This M[1]easure is a iteration. Check the calculations.
	double O=174.873; //The ecliptic longitude Ohm [large omega] of the ascending node of the orbit.
	double w=288.064; // The argument w[omega] of the perihelion
	double i=0.000; //the inclination i of the orbit
	double delta=0; //The distance large delta of the planet from the Earth
	double lambda=0;  //The ecliptic longitude
	double beta=0; //The ecliptic latitude
	double epsilon=0; //Obliquity of the ecliptic
	double alpha=0; //The right ascension
	double delta_low=0; //The declination
	double Oi[2]={0,0};        //Sidereal time on the prime meridian
	double theta[2]={0,0};
	double H=0; //Hour Angle
	
    JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)-1.5;   // At noon UTC on 1 January 2000, at Julian Day 2451545. 
	n = 0.9856076686 /(pow(1,(3/2)));  
	M[0]= Mo + (n * JD);
	M[1]=algo(M[0]);
	v[0]= deg_torad(M[1]) +(((2*e)-((1/4)*pow(e,3))+((5/96)*pow(e,5))+((107/4608)*pow(e,7)))*sin(deg_torad(M[1]))) 
	+ ((((5/4)*pow(e,2))-((11/24)*pow(e,4))+((17/192)*pow(e,6)))*sin(2*deg_torad(M[1])))
	+((((13/12)*pow(e,3))-((43/64)*pow(e,5))+((95/512)*pow(e,7)))*sin(3*deg_torad(M[1])))
	+((((103/96)*pow(e,4))-((451/480)*pow(e,6)))*sin(4*deg_torad(M[1])))
	+((((1097/960)*pow(e,5))-((5957/4680)*pow(e,7)))*sin(5*deg_torad(M[1])))
	+ (((1223/960)*pow(e,6))*sin(6*deg_torad(M[1])))
	+(((47273/32256)*pow(e,7))*sin(7*deg_torad(M[1])));
	v[1]=rad_todeg(v[0]);
	r=(1*(1-(pow(e,2))))/(1+(e*cos(deg_torad(v[1])))); // The distance between Earth and the Sun, measured in AU
	hX=r*((cos(deg_torad(O))*cos(deg_torad(w+v[1])))-(sin(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v[1])))); //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	hY=r*((sin(deg_torad(O))*cos(deg_torad(w+v[1])))+(cos(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v[1])))); //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	hZ=r*sin(deg_torad(i))*sin(deg_torad(w+v[1])); //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	// This fuction asume that XSun,YSun and ZSun are equal to zero due to the reference frame.
	gX=0-hX;
	gY=0-hY;
	gZ=0-hZ;
	delta=sqrt((pow(gX,2))+(pow(gY,2))+(pow(gZ,2))); //The distance large delta of the planet from the Earth
	lambda=atan2(gY,gX);    //The ecliptic longitude
	beta=asin(gZ/delta);    //The ecliptic latitude
    if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	epsilon=deg_torad(23.4397); //Obliquity of the ecliptic
	alpha=atan2(sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon),cos(lambda)); //The right ascension
	delta_low=asin((sin(beta)*cos(epsilon))+(cos(beta)*sin(epsilon)*sin(lambda))); //The declination
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}
	Oi[0]=M[1]+102.937+(15*(time+(5))); //5 plus hours for Colombia
	Oi[1]=algo(Oi[0]);
	theta[0]=Oi[1]+(lon); //The sideral time
	theta[1]=algo(theta[0]);
	H=theta[1]-rad_todeg(alpha);
	hei=asin(sin(deg_torad(lat))*sin(delta_low)+cos(deg_torad(lat))*cos(delta_low)*cos(deg_torad(H)));
	aziS=atan2(sin(deg_torad(H)),cos(deg_torad(H))*sin(deg_torad(lat))-tan(delta_low)*cos(deg_torad(lat)));
	if(aziS<=0)  /// Correction of the frame reference
	{
		aziN=aziS+deg_torad(180); /// Es la misma mierda pero no pude sacar un mejor condicional xd jjajaja
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}

    return rad_todeg(hei);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double azimuth(int a,int m,int d,double time)
{
	double lon=-73.370721; //Longitude
	double lat=5.865964; //Latitude
	double JD=0; //Julian Date Number
	double r=0; //The distance of the Earth from the Sun
	double hX=0; //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	double hY=0; //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	double hZ=0; //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	double gX=0; //Rectangular Geocentric Ecliptical Coordinates. X Axis. 
	double gY=0; //Rectangular Geocentric Ecliptical Coordinates. Y Axis.
	double gZ=0; //Rectangular Geocentric Ecliptical Coordinates. Z Axis.
	double hei=0; //The height
    double aziS=0; //The Azimuth from the South
	double aziN=0; //The Azimuth from the North
	double M[2]={0,0}; // Mean Anomaly variable
	double n=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	float e=0.01671; //Eccentricity of the Earth's orbit, an ellipse. 
	double v[2]={0,0}; //True AnoM[1]aly of the Earth. This M[1]easure is a iteration. Check the calculations.
	double O=174.873; //The ecliptic longitude Ohm [large omega] of the ascending node of the orbit.
	double w=288.064; // The argument w[omega] of the perihelion
	double i=0.000; //the inclination i of the orbit
	double delta=0; //The distance large delta of the planet from the Earth
	double lambda=0;  //The ecliptic longitude
	double beta=0; //The ecliptic latitude
	double epsilon=0; //Obliquity of the ecliptic
	double alpha=0; //The right ascension
	double delta_low=0; //The declination
	double Oi[2]={0,0}; //Sidereal time on the prime meridian
	double theta[2]={0,0};
	double H=0; //Hour Angle
	
    JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)-1.5;   // At noon UTC on 1 January 2000, at Julian Day 2451545. 
	n = 0.9856076686 /(pow(1,(3/2)));  
	M[0]= Mo + (n * JD);
	M[1]=algo(M[0]);
	v[0]= deg_torad(M[1]) +(((2*e)-((1/4)*pow(e,3))+((5/96)*pow(e,5))+((107/4608)*pow(e,7)))*sin(deg_torad(M[1]))) 
	+ ((((5/4)*pow(e,2))-((11/24)*pow(e,4))+((17/192)*pow(e,6)))*sin(2*deg_torad(M[1])))
	+((((13/12)*pow(e,3))-((43/64)*pow(e,5))+((95/512)*pow(e,7)))*sin(3*deg_torad(M[1])))
	+((((103/96)*pow(e,4))-((451/480)*pow(e,6)))*sin(4*deg_torad(M[1])))
	+((((1097/960)*pow(e,5))-((5957/4680)*pow(e,7)))*sin(5*deg_torad(M[1])))
	+ (((1223/960)*pow(e,6))*sin(6*deg_torad(M[1])))
	+(((47273/32256)*pow(e,7))*sin(7*deg_torad(M[1])));
	v[1]=rad_todeg(v[0]);
	r=(1*(1-(pow(e,2))))/(1+(e*cos(deg_torad(v[1])))); // The distance between Earth and the Sun, measured in AU
	hX=r*((cos(deg_torad(O))*cos(deg_torad(w+v[1])))-(sin(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v[1])))); //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	hY=r*((sin(deg_torad(O))*cos(deg_torad(w+v[1])))+(cos(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v[1])))); //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	hZ=r*sin(deg_torad(i))*sin(deg_torad(w+v[1])); //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	// This fuction asume that XSun,YSun and ZSun are equal to zero due to the reference frame.
	gX=0-hX;
	gY=0-hY;
	gZ=0-hZ;
	delta=sqrt((pow(gX,2))+(pow(gY,2))+(pow(gZ,2))); //The distance large delta of the planet from the Earth
	lambda=atan2(gY,gX);    //The ecliptic longitude
	beta=asin(gZ/delta);    //The ecliptic latitude
    if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	epsilon=deg_torad(23.4397); //Obliquity of the ecliptic
	alpha=atan2(sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon),cos(lambda)); //The right ascension
	delta_low=asin((sin(beta)*cos(epsilon))+(cos(beta)*sin(epsilon)*sin(lambda))); //The declination
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}
	Oi[0]=M[1]+102.937+(15*(time+(5))); //5 plus hours for Colombia
	Oi[1]=algo(Oi[0]);
	theta[0]=Oi[1]+(lon); //The sideral time
	theta[1]=algo(theta[0]);
	H=theta[1]-rad_todeg(alpha);
	hei=asin(sin(deg_torad(lat))*sin(delta_low)+cos(deg_torad(lat))*cos(delta_low)*cos(deg_torad(H)));
	aziS=atan2(sin(deg_torad(H)),cos(deg_torad(H))*sin(deg_torad(lat))-tan(delta_low)*cos(deg_torad(lat)));
	if(aziS<=0)  /// Correction of the frame reference
	{
		aziN=aziS+deg_torad(180); /// Es la misma mierda pero no pude sacar un mejor condicional xd jjajaja
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}

    return rad_todeg(aziN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////


////////

void daily_mode() //Daily mode Canvas generator
{
	
	int year=2019;
	int month[3];
	month[0]=3;
	month[1]=4;
	month[2]=5;
	int day[31];
	for(int a=0;a<31;a++)
	{  
	  day[a]=a+1;
	}

	
	string y="2019";
	string m[3];
	m[0]="03";
	m[1]="04";
	m[2]="05";
	string d[31];
    ostringstream vchange;
	for(int a=0;a<31;a++)
	{
		vchange=NULL;
		vchange<<day[a];
		d[a] = vchange.str(); 
	}
	
	

	
	gSystem->Exec("mkdir EcoacousterSun_Data");
	gSystem->Exec("mkdir EcoacousterSun_Motion");
	gSystem->Exec("mkdir EcoacousterSun_PolarSurvey");
	gSystem->Exec("mkdir EcoacousterSun_RF");
	
	string data_arranged="DailyMode_";
    string base_canvas1="DailyMotion_";
	string base_canvas2="SunSurvey_";
	
	string ref[3]={"Mar","Apr","May"};


    for(int box=0;box<3;box++)
	{

        int t=box;
		

    for(int s=0;s<31;s++)
	{
		
		double azimutal[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    double altitud[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	    ofstream daily;  //Data
	    daily.open(("C://root_v5.34.38/EcoacousterSun_Data/"+m[t]+"-"+ref[t]+"-"+data_arranged+"-"+d[s]+"-"+m[t]+"-"+y+".txt").c_str());
	    daily<<"****Data saved in degrees****"<<endl;
	    daily<<"************HDQS*************"<<endl;
	    daily<<"Date: "<<day[s]<<"/"<<month[t]<<"/"<<year<<endl;
	    daily<<"Hour"<<"	"<<"Azimuth"<<"	"<<"Height"<<endl;
	
	    for(int u=0;u<24;u++)
	    {
		    altitud[u]=height(year,month[t],day[s],u);
		    azimutal[u]=azimuth(year,month[t],day[s],u);
		    daily<<u<<"	"<<azimutal[u]<<"	"<<altitud[u]<<endl;
	    }
	
	    daily.close();
	

	    c1= new TCanvas("c1","Daily mode",1300,800);  //Motion
	    c1->SetTitle("Daily motion of the Sun");
        c1->Divide(3,1);

	
	    c1->cd(1);
	    a_ho= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    a_ho->SetPoint(u,u,azimutal[u]);
	    }
	    a_ho->SetTitle("Azimuth vs Time;Time(hours);Azimuth from the North(degrees)");
	    a_ho->SetLineColor(kRed);
	    a_ho->SetMarkerStyle(8);
	    a_ho->SetMarkerColor(kGreen+2);
	    a_ho->SetLineWidth(3);
	    a_ho->Draw("ACP");
	    a_ho->GetXaxis()->SetTitle("Time(hours)");
        a_ho->GetYaxis()->SetTitle("Azimuth from the North(degrees)");
        a_ho->GetXaxis()->CenterTitle();
        a_ho->GetYaxis()->CenterTitle();
	    a_ho->GetXaxis()->SetLabelSize(0.025);
	    a_ho->GetYaxis()->SetLabelSize(0.025);
	    a_ho->GetXaxis()->SetTitleOffset(1);
	    a_ho->GetYaxis()->SetTitleOffset(1.5);
		
		
	
	    c1->cd(2);
	    h_ho= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    h_ho->SetPoint(u,u,altitud[u]);
	    }
	    h_ho->SetTitle("Height vs Time;Time(hours);Height from the Horizon(degrees)");
	    h_ho->Draw("ACP");
	    h_ho->SetLineWidth(3);
	    h_ho->SetLineColor(kRed);
	    h_ho->SetMarkerStyle(8);
	    h_ho->SetMarkerColor(kGreen+2);
	    h_ho->GetXaxis()->SetTitle("Time(hours)");
        h_ho->GetYaxis()->SetTitle("Height from the Horizon(degrees)");
        h_ho->GetXaxis()->CenterTitle();
        h_ho->GetYaxis()->CenterTitle();
	    h_ho->GetXaxis()->SetLabelSize(0.025);
	    h_ho->GetYaxis()->SetLabelSize(0.025);
	    h_ho->GetXaxis()->SetTitleOffset(1);
	    h_ho->GetYaxis()->SetTitleOffset(1.5);
	
	    c1->cd(3);
	    a_h= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    a_h->SetPoint(u,altitud[u],azimutal[u]);
	    }
	    a_h->SetTitle("Azimuth vs Height;Height(degrees);Azimuth(degrees)");
	    a_h->SetLineColor(kRed);
	    a_h->SetMarkerStyle(8);
        a_h->SetMarkerColor(kGreen+2);
	    a_h->SetLineWidth(3);
	    a_h->Draw("ACP");
	    a_h->GetXaxis()->SetTitle("Height(degrees)");
        a_h->GetYaxis()->SetTitle("Azimuth(degrees)");
        a_h->GetXaxis()->CenterTitle();
        a_h->GetYaxis()->CenterTitle();
	    a_h->GetXaxis()->SetLabelSize(0.025);
	    a_h->GetYaxis()->SetLabelSize(0.025);
	    a_h->GetXaxis()->SetTitleOffset(1);
	    a_h->GetYaxis()->SetTitleOffset(1.5);
		


	
   	    c1->SaveAs(("C://root_v5.34.38/EcoacousterSun_Motion/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		c1->SaveAs(("C://root_v5.34.38/EcoacousterSun_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());




	    c2= new TCanvas("c2","Daily Polar mode",1300,750);  //Survey
	
	    c2->cd(1);
        PL= new TGraphPolar();
        PL->SetTitle("Sun Survey through the Sky");
	    for(int u=0;u<24;u++)
	    {
		    PL->SetPoint(u,deg_torad(360-azimutal[u]),altitud[u]);
	    }  
        PL->SetMarkerStyle(20);
        PL->SetMarkerSize(1.);
        PL->SetMarkerColor(4);
        PL->SetLineColor(4);
        PL->Draw("ACP");
		// Update, otherwise GetPolargram returns 0 and the graph get wrong.
		c2->Update();
        PL->GetPolargram()->SetToDegree();
		PL->GetPolargram()->SetPolarOffset(0.01);
		PL->GetPolargram()->SetPolarLabelSize(0.03);
		
		
		
		//Textual conventions
		c2->cd(2); 
        SubT = new TText(0.8,1,("Date: "+d[s]+"/"+m[t]+"/"+y).c_str());
        SubT->SetTextAlign(22);
        SubT->SetTextFont(43);
        SubT->SetTextSize(20);
        SubT->Draw();
		
		c2->cd(3);
        CPN = new TText(1.07,0,"N");
        CPN->SetTextAlign(22);
        CPN->SetTextFont(23);
        CPN->SetTextSize(20);
        CPN->Draw();
		
		c2->cd(4);
        CPE = new TText(0,-1.1,"E");
        CPE->SetTextAlign(22);
        CPE->SetTextFont(23);
        CPE->SetTextSize(20);
        CPE->Draw();
		
		c2->cd(5);
        CPS = new TText(-1.12,0,"S");
        CPS->SetTextAlign(22);
        CPS->SetTextFont(23);
        CPS->SetTextSize(20);
        CPS->Draw();
		
		c2->cd(6);
        CPW = new TText(0.05,1.07,"W");
        CPW->SetTextAlign(22);
        CPW->SetTextFont(23);
        CPW->SetTextSize(20);
        CPW->Draw();
		
		
		
	
	    c2->SaveAs(("C://root_v5.34.38/EcoacousterSun_PolarSurvey/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		c2->SaveAs(("C://root_v5.34.38/EcoacousterSun_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());
		
		
	}
	}
	

	
}




void wholeday_mode() //Whole day year mode Canvas generator
{
	
	int year=2019;
	int month[12];
	for(int a=0;a<12;a++)
	{  
	  month[a]=a+1;
	}
	int day[31];
	for(int a=0;a<31;a++)
	{  
	  day[a]=a+1;
	}

	
	string y="2019";
	string m[12];
    ostringstream vchangem;
	for(int a=0;a<12;a++)
	{
		vchangem=NULL;
		vchangem<<month[a];
		m[a] = vchangem.str(); 
	}
	string d[31];
    ostringstream vchanged;
	for(int a=0;a<31;a++)
	{
		vchanged=NULL;
		vchanged<<day[a];
		d[a] = vchanged.str(); 
	}
	
	

	
	gSystem->Exec("mkdir EcoacousterSun_Data");
	gSystem->Exec("mkdir EcoacousterSun_Motion");
	gSystem->Exec("mkdir EcoacousterSun_PolarSurvey");
	gSystem->Exec("mkdir EcoacousterSun_RF");
	
	string data_arranged="DailyMode_";
    string base_canvas1="DailyMotion_";
	string base_canvas2="SunSurvey_";
	
	string ref[12]={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};


    for(int box=0;box<12;box++)
	{

        int t=box;
		

    for(int s=0;s<31;s++)
	{
		
		if((t==1 && s==28)||(t==1 && s==29)||(t==1 && s==30)||(t==3 && s==30)||(t==5 && s==30)||(t==8 && s==30)||(t==10 && s==30))// These days don't exist
		{
			ofstream daily;  //Data
	        daily.open(("C://root_v5.34.38/EcoacousterSun_Data/"+m[t]+"-"+ref[t]+"-"+"No Data"+"-"+d[s]+"-"+m[t]+"-"+y+".txt").c_str());
	        daily<<"****Data saved in degrees****"<<endl;
	        daily<<"************HDQS*************"<<endl;
			daily<<"**********No Data***********"<<endl;
			
		}
		else
		{
			double azimutal[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	        double altitud[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	        ofstream daily;  //Data
	        daily.open(("C://root_v5.34.38/EcoacousterSun_Data/"+m[t]+"-"+ref[t]+"-"+data_arranged+"-"+d[s]+"-"+m[t]+"-"+y+".txt").c_str());
	        daily<<"****Data saved in degrees****"<<endl;
	        daily<<"************HDQS*************"<<endl;
	        daily<<"Date: "<<day[s]<<"/"<<month[t]<<"/"<<year<<endl;
	        daily<<"Hour"<<"	"<<"Azimuth"<<"	"<<"Height"<<endl;
	
	        for(int u=0;u<24;u++)
	        {
		        altitud[u]=height(year,month[t],day[s],u);
		        azimutal[u]=azimuth(year,month[t],day[s],u);
		        daily<<u<<"	"<<azimutal[u]<<"	"<<altitud[u]<<endl;
	        }
	
	        daily.close();
	

	        c1= new TCanvas("c1","Daily mode",1300,800);  //Motion
	        c1->SetTitle("Daily motion of the Sun");
            c1->Divide(3,1);

	
	        c1->cd(1);
	        a_ho= new TGraph();
	        for(int u=0;u<24;u++)
	        {
		        a_ho->SetPoint(u,u,azimutal[u]);
	        }
	        a_ho->SetTitle("Azimuth vs Time;Time(hours);Azimuth from the North(degrees)");
	        a_ho->SetLineColor(kRed);
	        a_ho->SetMarkerStyle(8);
	        a_ho->SetMarkerColor(kGreen+2);
	        a_ho->SetLineWidth(3);
	        a_ho->Draw("ACP");
	        a_ho->GetXaxis()->SetTitle("Time(hours)");
            a_ho->GetYaxis()->SetTitle("Azimuth from the North(degrees)");
            a_ho->GetXaxis()->CenterTitle();
            a_ho->GetYaxis()->CenterTitle();
	        a_ho->GetXaxis()->SetLabelSize(0.025);
	        a_ho->GetYaxis()->SetLabelSize(0.025);
	        a_ho->GetXaxis()->SetTitleOffset(1);
	        a_ho->GetYaxis()->SetTitleOffset(1.5);
		
		
	
	        c1->cd(2);
	        h_ho= new TGraph();
	        for(int u=0;u<24;u++)
	        {
		        h_ho->SetPoint(u,u,altitud[u]);
	        }
	        h_ho->SetTitle("Height vs Time;Time(hours);Height from the Horizon(degrees)");
	        h_ho->Draw("ACP");
	        h_ho->SetLineWidth(3);
	        h_ho->SetLineColor(kRed);
	        h_ho->SetMarkerStyle(8);
	        h_ho->SetMarkerColor(kGreen+2);
	        h_ho->GetXaxis()->SetTitle("Time(hours)");
            h_ho->GetYaxis()->SetTitle("Height from the Horizon(degrees)");
            h_ho->GetXaxis()->CenterTitle();
            h_ho->GetYaxis()->CenterTitle();
	        h_ho->GetXaxis()->SetLabelSize(0.025);
	        h_ho->GetYaxis()->SetLabelSize(0.025);
	        h_ho->GetXaxis()->SetTitleOffset(1);
	        h_ho->GetYaxis()->SetTitleOffset(1.5);
	
	        c1->cd(3);
	        a_h= new TGraph();
	        for(int u=0;u<24;u++)
	        {
		        a_h->SetPoint(u,altitud[u],azimutal[u]);
	        }
	        a_h->SetTitle("Azimuth vs Height;Height(degrees);Azimuth(degrees)");
	        a_h->SetLineColor(kRed);
	        a_h->SetMarkerStyle(8);
            a_h->SetMarkerColor(kGreen+2);
	        a_h->SetLineWidth(3);
	        a_h->Draw("ACP");
	        a_h->GetXaxis()->SetTitle("Height(degrees)");
            a_h->GetYaxis()->SetTitle("Azimuth(degrees)");
            a_h->GetXaxis()->CenterTitle();
            a_h->GetYaxis()->CenterTitle();
	        a_h->GetXaxis()->SetLabelSize(0.025);
	        a_h->GetYaxis()->SetLabelSize(0.025);
	        a_h->GetXaxis()->SetTitleOffset(1);
	        a_h->GetYaxis()->SetTitleOffset(1.5);
		


	
   	        c1->SaveAs(("C://root_v5.34.38/EcoacousterSun_Motion/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		    c1->SaveAs(("C://root_v5.34.38/EcoacousterSun_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());




	        c2= new TCanvas("c2","Daily Polar mode",1300,750);  //Survey
	
	        c2->cd(1);
            PL= new TGraphPolar();
            PL->SetTitle("Sun Survey through the Sky");
	        for(int u=0;u<24;u++)
	        {
		        PL->SetPoint(u,deg_torad(360-azimutal[u]),altitud[u]);
	        }  
            PL->SetMarkerStyle(20);
            PL->SetMarkerSize(1.);
            PL->SetMarkerColor(4);
            PL->SetLineColor(4);
            PL->Draw("ACP");
		    // Update, otherwise GetPolargram returns 0 and the graph get wrong.
		    c2->Update();
            PL->GetPolargram()->SetToDegree();
		    PL->GetPolargram()->SetPolarOffset(0.01);
		    PL->GetPolargram()->SetPolarLabelSize(0.03);
		
		
		
		    //Textual conventions
		    c2->cd(2); 
            SubT = new TText(0.8,1,("Date: "+d[s]+"/"+m[t]+"/"+y).c_str());
            SubT->SetTextAlign(22);
            SubT->SetTextFont(43);
            SubT->SetTextSize(20);
            SubT->Draw();
		
		    c2->cd(3);
            CPN = new TText(1.07,0,"N");
            CPN->SetTextAlign(22);
            CPN->SetTextFont(23);
            CPN->SetTextSize(20);
            CPN->Draw();
		
		    c2->cd(4);
            CPE = new TText(0,-1.1,"E");
            CPE->SetTextAlign(22);
            CPE->SetTextFont(23);
            CPE->SetTextSize(20);
            CPE->Draw();
		
		    c2->cd(5);
            CPS = new TText(-1.12,0,"S");
            CPS->SetTextAlign(22);
            CPS->SetTextFont(23);
            CPS->SetTextSize(20);
            CPS->Draw();
		
		    c2->cd(6);
            CPW = new TText(0.05,1.07,"W");
            CPW->SetTextAlign(22);
            CPW->SetTextFont(23);
            CPW->SetTextSize(20);
            CPW->Draw();
		
		
		
	
	        c2->SaveAs(("C://root_v5.34.38/EcoacousterSun_PolarSurvey/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		    c2->SaveAs(("C://root_v5.34.38/EcoacousterSun_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());
			
			
			
			
			
			
			
		}
		
		
	}
	}
	

	
}

void stat_testSUN() //Sampling for statistical test comparision
{
	int y=2019;
	int m[240];
	int d[240];
	double h[240];
	//double min[240];
	
	ifstream sample1;
	sample1.open("Month(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample1>>m[i];
	}
	sample1.close();
	
	ifstream sample2;
	sample2.open("Day(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample2>>d[i];
	}
	sample2.close();
	
	ifstream sample3;
	sample3.open("Hour(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample3>>h[i];
	}
	sample3.close();
	/*
	ifstream sample4;
	sample4.open("Minute(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample4>>min[i];
	}
	sample4.close();
	
	for(int a=0;a<240;a++)
	{
		h[a]=h[a]+(min[a]/60);
	}
	*/
	ofstream TEST1;  //Data output
	TEST1.open("Statistical test Sun Azimuth.txt");
    for(int i=0;i<240;i++)
	{
		TEST1<<azimuth(y,m[i],d[i],h[i])<<endl;
	}
	TEST1.close();
	
	ofstream TEST2;  //Data output
	TEST2.open("Statistical test Sun Height.txt");
    for(int i=0;i<240;i++)
	{
		TEST2<<height(y,m[i],d[i],h[i])<<endl;;
	}
	TEST2.close();
	
	cout<<"DATA SAVED!!!"<<endl;
	
}

//****************************************************************************************************************************************************
//MOON CALCULATIONS CODE******************************************************************************************************************************
//****************************************************************************************************************************************************

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double height_moon(int a,int m,int d,double time)
{
	double lon=-73.370721; //Longitude
	double lat=5.865964; //Latitude
	double JD=0; //Julian Date Number from epoch 2010 January 0
	double RJD=0;//The true Julian Date Number
	double T=0; //Number of Julian Centuries since 2000 January 0.5
	double eg=0;//The Sun’s mean ecliptic longitude measured in degrees
	double wg=0;//The longitude of the Sun at perigee measured in degrees
	double e=0;//The eccentricity of the Sun–Earth orbit
	double Ms=0;//Mean Anomaly of the Sun
	double E[2]={0,0};//The eccentric anomaly of the Sun measured in radians/degrees
	double v=0;//True anomaly of the Sun in degrees
	double Ls=0; //The ecliptic longitude of the Sun in degrees
	double Lm=0;//Moon’s mean longitude in degrees
	double Mm=0;//The Moon’s mean anomaly
	double N=0;//The Moon’s ascending node's mean longitude
	double C=0; //Constant of correction
	double Ev=0; // Correction for evection measured in degrees
	double Ae=0; //Annual equation correction measured in degrees
	double A3=0; // A third correction measured in degrees
	double MmR=0;// Moon's corrected anomaly measured in degrees
	double Ec=0; //Correction fo the equation of centre measured in degrees
	double A4=0; // A fourth correction measured in degrees
	double LmR=0; // Value of the Moon's corrected longitude measured in degrees
	double VC=0; //The final correction to apply to the Moon's longitude, the variation measured in degrees
	double LmR2=0; //The Moon's true orbital longitude measured in degrees
	double NR=0; //Corrected longitude of the node measured in degrees
	double lambda=0;//Ecliptic longitude Lambda
	double beta=0;//Ecliptic latitude Beta
	double epsilon=0;//The obliquity of the ecliptic, the angle between the planes of the equator and the ecliptic measured in degrees
	double alpha=0;//The right ascension of the ecuatorial plane
	double delta=0;//The declination of the ecuatorial plane
	double Me=0; // Mean Anomaly of the Earth
	double n=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Oi[2]={0,0}; //Sidereal time on the prime meridian + 5 plus hours for Colombia
	double theta[2]={0,0};//The sideral time
	double H=0; //Hour Angle
	double D=0; //Change in the Moon's orbit from the line of reference through the month
	double F=0; ////Moon's phase
	double hei=0; //The height of the Moon
    double aziS=0; //The Moon's Azimuth from the South
	double aziN=0; //The Moon's Azimuth from the North
	
	double Lm0=91.929336; //Moon’s mean longitude at the epoch  measured in degrees
	double P0=130.143076; //Moon’s mean longitude of the perigee  at the epoch measured in degrees
	double N0=291.682547; //Moon’s mean longitude of the node  at the epoch measured in degrees
	double i=5.145396; //Inclination of Moon’s orbit measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	
	int dm[3]={0,0,0};//Dummy 1
	double abc[3]={0,0,0};//Dummy 2
	double right=0;//Dummy 3
	double y=0; //Dummy 4
	double x=0; //Dummy 5
	double k=0; //Dummy 6
	double arcs_deg[2]={0,0};//Dummy 7
	double b=0; //Dummy 8
	
	JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)+(time/24);   // At midnight UTC on 31 December of 2010, at Julian Day 2455196.5. Check http://www.onlineconversion.com/julian_date.htm
	JD=JD-3653;
	dm[0]=(m+9)/12;
    dm[1]=(7*(a+dm[0]))/4;
	dm[2]=(275*m)/9;
	RJD=(367*a)-dm[1]+dm[2]+d+1721013.5+(time/24);
	T=(RJD-2451545.0)/36525.0;
	eg=279.6966778+(36000.76892*T)+(0.0003025*pow(T,2));
	if(T==1.09997262149212860) //The Sun’s mean ecliptic longitude at epoch 2010 measured in degrees
	{
		eg=279.557208;
	}
	else
	{
		eg=sp_algo(eg);
	}
	wg=281.2208444+(1.719175*T)+(0.000452778*pow(T,2));
	if(T==1.09997262149212860) //The longitude of the Sun at perigee at epoch 2010 measured in degrees
	{
		wg=283.112438;
	}
	else
	{
		wg=sp_algo(wg);
	}
	e=0.01675104-(0.0000418*T)-(0.000000126*pow(T,2));
	if(T==1.09997262149212860) //The eccentricity of the Sun–Earth orbit at epoch 2010 measured in degrees
	{
		e=0.016705;
	}
	else
	{
		e=sp_algo(e);
	}
	Ms=eg-wg;
	Ms=sp_algo(Ms);
	E[0]=E_f(e,Ms); //Measured in radians
	E[1]=rad_todeg(E[0]);//Measured in degrees
	abc[0]=1+e;
	abc[1]=1-e;
	abc[2]=sqrt(abc[0]/abc[1]);
    right=abc[2]*tan(E[0]/2); // E measured in radians
	v=atan(right)*2;
	v=sp_algo(rad_todeg(v));
	Ls=v+wg;
	Ls=sp_algo(Ls);
    Lm=(13.1763966*JD)+Lm0;
	Lm=sp_algo(Lm);
	Mm=Lm-(0.1114041*JD)-P0;
	Mm=sp_algo(Mm);
	N=N0-(0.0529539*JD);
	N=sp_algo(N);
	C=Lm-Ls; 
	Ev=1.2739*sin(deg_torad((2*C)-Mm));
	Ae=0.1858*sin(deg_torad(Ms)); 
	A3=0.37*sin(deg_torad(Ms));
	MmR=Mm+Ev-Ae-A3;
	Ec=6.2886*sin(deg_torad(MmR)); 
	A4=0.214*sin(deg_torad(2*MmR));
	LmR=Lm+Ev+Ec-Ae+A4; 
	VC=0.6583*sin(deg_torad(2*(LmR-Ls))); 
	LmR2=LmR+VC;
	NR=N-(0.16*(sin(deg_torad(Ms)))); 
	y=cos(deg_torad(i))*sin(deg_torad(LmR2-NR));
	x=cos(deg_torad(LmR2-NR));
	lambda=atan2(y,x)+deg_torad(NR);
	if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	k=sin(deg_torad(i))*sin(deg_torad(LmR2-NR));
	beta=asin(k);
	arcs_deg[0]=(46.815*T)+(0.0006*(T*T))-(0.00181*(T*T*T));
	arcs_deg[1]=arcs_deg[0]/3600;
	epsilon=23.43929167-arcs_deg[1];
    epsilon=epsilon+0.02755;
	b=(sin(lambda)*cos(deg_torad(epsilon)))-(tan(beta)*sin(deg_torad(epsilon)));
	alpha=atan2(b,cos(lambda));
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}
	delta=asin((sin(beta)*cos(deg_torad(epsilon)))+(cos(beta)*sin(deg_torad(epsilon))*sin(lambda)));
	n = 0.9856076686 /(pow(1,(3/2))); 
	Me= Mo + (n * JD);
	Me=sp_algo(Me);
	Oi[0]=Me+102.937+(15*(time+(5))); 
	Oi[1]=sp_algo(Oi[0]);
	theta[0]=Oi[1]+(lon);
	theta[1]=sp_algo(theta[0]); 
	H=theta[1]-rad_todeg(alpha);
	D=LmR2-Ls;
    D=sp_algo(D);
    F=(0.5)*(1-cos(deg_torad(D)));
	hei=asin(sin(deg_torad(lat))*sin(delta)+cos(deg_torad(lat))*cos(delta)*cos(deg_torad(H)));
	aziS=atan2(sin(deg_torad(H)),cos(deg_torad(H))*sin(deg_torad(lat))-tan(delta)*cos(deg_torad(lat)));
	if(aziS<=0)
	{
		aziN=aziS+deg_torad(180);
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}
	

	return rad_todeg(hei);
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double azimuth_moon(int a,int m,int d,double time)
{
	double lon=-73.370721; //Longitude
	double lat=5.865964; //Latitude
	double JD=0; //Julian Date Number from epoch 2010 January 0
	double RJD=0;//The true Julian Date Number
	double T=0; //Number of Julian Centuries since 2000 January 0.5
	double eg=0;//The Sun’s mean ecliptic longitude measured in degrees
	double wg=0;//The longitude of the Sun at perigee measured in degrees
	double e=0;//The eccentricity of the Sun–Earth orbit
	double Ms=0;//Mean Anomaly of the Sun
	double E[2]={0,0};//The eccentric anomaly of the Sun measured in radians/degrees
	double v=0;//True anomaly of the Sun in degrees
	double Ls=0; //The ecliptic longitude of the Sun in degrees
	double Lm=0;//Moon’s mean longitude in degrees
	double Mm=0;//The Moon’s mean anomaly
	double N=0;//The Moon’s ascending node's mean longitude
	double C=0; //Constant of correction
	double Ev=0; // Correction for evection measured in degrees
	double Ae=0; //Annual equation correction measured in degrees
	double A3=0; // A third correction measured in degrees
	double MmR=0;// Moon's corrected anomaly measured in degrees
	double Ec=0; //Correction fo the equation of centre measured in degrees
	double A4=0; // A fourth correction measured in degrees
	double LmR=0; // Value of the Moon's corrected longitude measured in degrees
	double VC=0; //The final correction to apply to the Moon's longitude, the variation measured in degrees
	double LmR2=0; //The Moon's true orbital longitude measured in degrees
	double NR=0; //Corrected longitude of the node measured in degrees
	double lambda=0;//Ecliptic longitude Lambda
	double beta=0;//Ecliptic latitude Beta
	double epsilon=0;//The obliquity of the ecliptic, the angle between the planes of the equator and the ecliptic measured in degrees
	double alpha=0;//The right ascension of the ecuatorial plane
	double delta=0;//The declination of the ecuatorial plane
	double Me=0; // Mean Anomaly of the Earth
	double n=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Oi[2]={0,0}; //Sidereal time on the prime meridian + 5 plus hours for Colombia
	double theta[2]={0,0};//The sideral time
	double H=0; //Hour Angle
	double D=0; //Change in the Moon's orbit from the line of reference through the month
	double F=0; ////Moon's phase
	double hei=0; //The height of the Moon
    double aziS=0; //The Moon's Azimuth from the South
	double aziN=0; //The Moon's Azimuth from the North
	
	double Lm0=91.929336; //Moon’s mean longitude at the epoch  measured in degrees
	double P0=130.143076; //Moon’s mean longitude of the perigee  at the epoch measured in degrees
	double N0=291.682547; //Moon’s mean longitude of the node  at the epoch measured in degrees
	double i=5.145396; //Inclination of Moon’s orbit measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	
	int dm[3]={0,0,0};//Dummy 1
	double abc[3]={0,0,0};//Dummy 2
	double right=0;//Dummy 3
	double y=0; //Dummy 4
	double x=0; //Dummy 5
	double k=0; //Dummy 6
	double arcs_deg[2]={0,0};//Dummy 7
	double b=0; //Dummy 8
	
	JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)+(time/24);   // At midnight UTC on 31 December of 2010, at Julian Day 2455196.5. Check http://www.onlineconversion.com/julian_date.htm
	JD=JD-3653;
	dm[0]=(m+9)/12;
    dm[1]=(7*(a+dm[0]))/4;
	dm[2]=(275*m)/9;
	RJD=(367*a)-dm[1]+dm[2]+d+1721013.5+(time/24);
	T=(RJD-2451545.0)/36525.0;
	eg=279.6966778+(36000.76892*T)+(0.0003025*pow(T,2));
	if(T==1.09997262149212860) //The Sun’s mean ecliptic longitude at epoch 2010 measured in degrees
	{
		eg=279.557208;
	}
	else
	{
		eg=sp_algo(eg);
	}
	wg=281.2208444+(1.719175*T)+(0.000452778*pow(T,2));
	if(T==1.09997262149212860) //The longitude of the Sun at perigee at epoch 2010 measured in degrees
	{
		wg=283.112438;
	}
	else
	{
		wg=sp_algo(wg);
	}
	e=0.01675104-(0.0000418*T)-(0.000000126*pow(T,2));
	if(T==1.09997262149212860) //The eccentricity of the Sun–Earth orbit at epoch 2010 measured in degrees
	{
		e=0.016705;
	}
	else
	{
		e=sp_algo(e);
	}
	Ms=eg-wg;
	Ms=sp_algo(Ms);
	E[0]=E_f(e,Ms); //Measured in radians
	E[1]=rad_todeg(E[0]);//Measured in degrees
	abc[0]=1+e;
	abc[1]=1-e;
	abc[2]=sqrt(abc[0]/abc[1]);
    right=abc[2]*tan(E[0]/2); // E measured in radians
	v=atan(right)*2;
	v=sp_algo(rad_todeg(v));
	Ls=v+wg;
	Ls=sp_algo(Ls);
    Lm=(13.1763966*JD)+Lm0;
	Lm=sp_algo(Lm);
	Mm=Lm-(0.1114041*JD)-P0;
	Mm=sp_algo(Mm);
	N=N0-(0.0529539*JD);
	N=sp_algo(N);
	C=Lm-Ls; 
	Ev=1.2739*sin(deg_torad((2*C)-Mm));
	Ae=0.1858*sin(deg_torad(Ms)); 
	A3=0.37*sin(deg_torad(Ms));
	MmR=Mm+Ev-Ae-A3;
	Ec=6.2886*sin(deg_torad(MmR)); 
	A4=0.214*sin(deg_torad(2*MmR));
	LmR=Lm+Ev+Ec-Ae+A4; 
	VC=0.6583*sin(deg_torad(2*(LmR-Ls))); 
	LmR2=LmR+VC;
	NR=N-(0.16*(sin(deg_torad(Ms)))); 
	y=cos(deg_torad(i))*sin(deg_torad(LmR2-NR));
	x=cos(deg_torad(LmR2-NR));
	lambda=atan2(y,x)+deg_torad(NR);
	if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	k=sin(deg_torad(i))*sin(deg_torad(LmR2-NR));
	beta=asin(k);
	arcs_deg[0]=(46.815*T)+(0.0006*(T*T))-(0.00181*(T*T*T));
	arcs_deg[1]=arcs_deg[0]/3600;
	epsilon=23.43929167-arcs_deg[1];
    epsilon=epsilon+0.02755;
	b=(sin(lambda)*cos(deg_torad(epsilon)))-(tan(beta)*sin(deg_torad(epsilon)));
	alpha=atan2(b,cos(lambda));
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}
	delta=asin((sin(beta)*cos(deg_torad(epsilon)))+(cos(beta)*sin(deg_torad(epsilon))*sin(lambda)));
	n = 0.9856076686 /(pow(1,(3/2))); 
	Me= Mo + (n * JD);
	Me=sp_algo(Me);
	Oi[0]=Me+102.937+(15*(time+(5))); 
	Oi[1]=sp_algo(Oi[0]);
	theta[0]=Oi[1]+(lon);
	theta[1]=sp_algo(theta[0]); 
	H=theta[1]-rad_todeg(alpha);
	D=LmR2-Ls;
    D=sp_algo(D);
    F=(0.5)*(1-cos(deg_torad(D)));
	hei=asin(sin(deg_torad(lat))*sin(delta)+cos(deg_torad(lat))*cos(delta)*cos(deg_torad(H)));
	aziS=atan2(sin(deg_torad(H)),cos(deg_torad(H))*sin(deg_torad(lat))-tan(delta)*cos(deg_torad(lat)));
	if(aziS<=0)
	{
		aziN=aziS+deg_torad(180);
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}
	

	return rad_todeg(aziN);
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double moon_phase(int a,int m,int d,double time)
{
	double lon=-73.370721; //Longitude
	double lat=5.865964; //Latitude
	double JD=0; //Julian Date Number from epoch 2010 January 0
	double RJD=0;//The true Julian Date Number
	double T=0; //Number of Julian Centuries since 2000 January 0.5
	double eg=0;//The Sun’s mean ecliptic longitude measured in degrees
	double wg=0;//The longitude of the Sun at perigee measured in degrees
	double e=0;//The eccentricity of the Sun–Earth orbit
	double Ms=0;//Mean Anomaly of the Sun
	double E[2]={0,0};//The eccentric anomaly of the Sun measured in radians/degrees
	double v=0;//True anomaly of the Sun in degrees
	double Ls=0; //The ecliptic longitude of the Sun in degrees
	double Lm=0;//Moon’s mean longitude in degrees
	double Mm=0;//The Moon’s mean anomaly
	double N=0;//The Moon’s ascending node's mean longitude
	double C=0; //Constant of correction
	double Ev=0; // Correction for evection measured in degrees
	double Ae=0; //Annual equation correction measured in degrees
	double A3=0; // A third correction measured in degrees
	double MmR=0;// Moon's corrected anomaly measured in degrees
	double Ec=0; //Correction fo the equation of centre measured in degrees
	double A4=0; // A fourth correction measured in degrees
	double LmR=0; // Value of the Moon's corrected longitude measured in degrees
	double VC=0; //The final correction to apply to the Moon's longitude, the variation measured in degrees
	double LmR2=0; //The Moon's true orbital longitude measured in degrees
	double NR=0; //Corrected longitude of the node measured in degrees
	double lambda=0;//Ecliptic longitude Lambda
	double beta=0;//Ecliptic latitude Beta
	double epsilon=0;//The obliquity of the ecliptic, the angle between the planes of the equator and the ecliptic measured in degrees
	double alpha=0;//The right ascension of the ecuatorial plane
	double delta=0;//The declination of the ecuatorial plane
	double Me=0; // Mean Anomaly of the Earth
	double n=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Oi[2]={0,0}; //Sidereal time on the prime meridian + 5 plus hours for Colombia
	double theta[2]={0,0};//The sideral time
	double H=0; //Hour Angle
	double D=0; //Change in the Moon's orbit from the line of reference through the month
	double F=0; ////Moon's phase
	double hei=0; //The height of the Moon
    double aziS=0; //The Moon's Azimuth from the South
	double aziN=0; //The Moon's Azimuth from the North
	
	double Lm0=91.929336; //Moon’s mean longitude at the epoch  measured in degrees
	double P0=130.143076; //Moon’s mean longitude of the perigee  at the epoch measured in degrees
	double N0=291.682547; //Moon’s mean longitude of the node  at the epoch measured in degrees
	double i=5.145396; //Inclination of Moon’s orbit measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	
	int dm[3]={0,0,0};//Dummy 1
	double abc[3]={0,0,0};//Dummy 2
	double right=0;//Dummy 3
	double y=0; //Dummy 4
	double x=0; //Dummy 5
	double k=0; //Dummy 6
	double arcs_deg[2]={0,0};//Dummy 7
	double b=0; //Dummy 8
	
	JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)+(time/24);   // At midnight UTC on 31 December of 2010, at Julian Day 2455196.5. Check http://www.onlineconversion.com/julian_date.htm
	JD=JD-3653;
	dm[0]=(m+9)/12;
    dm[1]=(7*(a+dm[0]))/4;
	dm[2]=(275*m)/9;
	RJD=(367*a)-dm[1]+dm[2]+d+1721013.5+(time/24);
	T=(RJD-2451545.0)/36525.0;
	eg=279.6966778+(36000.76892*T)+(0.0003025*pow(T,2));
	if(T==1.09997262149212860) //The Sun’s mean ecliptic longitude at epoch 2010 measured in degrees
	{
		eg=279.557208;
	}
	else
	{
		eg=sp_algo(eg);
	}
	wg=281.2208444+(1.719175*T)+(0.000452778*pow(T,2));
	if(T==1.09997262149212860) //The longitude of the Sun at perigee at epoch 2010 measured in degrees
	{
		wg=283.112438;
	}
	else
	{
		wg=sp_algo(wg);
	}
	e=0.01675104-(0.0000418*T)-(0.000000126*pow(T,2));
	if(T==1.09997262149212860) //The eccentricity of the Sun–Earth orbit at epoch 2010 measured in degrees
	{
		e=0.016705;
	}
	else
	{
		e=sp_algo(e);
	}
	Ms=eg-wg;
	Ms=sp_algo(Ms);
	E[0]=E_f(e,Ms); //Measured in radians
	E[1]=rad_todeg(E[0]);//Measured in degrees
	abc[0]=1+e;
	abc[1]=1-e;
	abc[2]=sqrt(abc[0]/abc[1]);
    right=abc[2]*tan(E[0]/2); // E measured in radians
	v=atan(right)*2;
	v=sp_algo(rad_todeg(v));
	Ls=v+wg;
	Ls=sp_algo(Ls);
    Lm=(13.1763966*JD)+Lm0;
	Lm=sp_algo(Lm);
	Mm=Lm-(0.1114041*JD)-P0;
	Mm=sp_algo(Mm);
	N=N0-(0.0529539*JD);
	N=sp_algo(N);
	C=Lm-Ls; 
	Ev=1.2739*sin(deg_torad((2*C)-Mm));
	Ae=0.1858*sin(deg_torad(Ms)); 
	A3=0.37*sin(deg_torad(Ms));
	MmR=Mm+Ev-Ae-A3;
	Ec=6.2886*sin(deg_torad(MmR)); 
	A4=0.214*sin(deg_torad(2*MmR));
	LmR=Lm+Ev+Ec-Ae+A4; 
	VC=0.6583*sin(deg_torad(2*(LmR-Ls))); 
	LmR2=LmR+VC;
	NR=N-(0.16*(sin(deg_torad(Ms)))); 
	y=cos(deg_torad(i))*sin(deg_torad(LmR2-NR));
	x=cos(deg_torad(LmR2-NR));
	lambda=atan2(y,x)+deg_torad(NR);
	if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	k=sin(deg_torad(i))*sin(deg_torad(LmR2-NR));
	beta=asin(k);
	arcs_deg[0]=(46.815*T)+(0.0006*(T*T))-(0.00181*(T*T*T));
	arcs_deg[1]=arcs_deg[0]/3600;
	epsilon=23.43929167-arcs_deg[1];
    epsilon=epsilon+0.02755;
	b=(sin(lambda)*cos(deg_torad(epsilon)))-(tan(beta)*sin(deg_torad(epsilon)));
	alpha=atan2(b,cos(lambda));
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}
	delta=asin((sin(beta)*cos(deg_torad(epsilon)))+(cos(beta)*sin(deg_torad(epsilon))*sin(lambda)));
	n = 0.9856076686 /(pow(1,(3/2))); 
	Me= Mo + (n * JD);
	Me=sp_algo(Me);
	Oi[0]=Me+102.937+(15*(time+(5))); 
	Oi[1]=sp_algo(Oi[0]);
	theta[0]=Oi[1]+(lon);
	theta[1]=sp_algo(theta[0]); 
	H=theta[1]-rad_todeg(alpha);
	D=LmR2-Ls;
    D=sp_algo(D);
    F=(0.5)*(1-cos(deg_torad(D)));
	hei=asin(sin(deg_torad(lat))*sin(delta)+cos(deg_torad(lat))*cos(delta)*cos(deg_torad(H)));
	aziS=atan2(sin(deg_torad(H)),cos(deg_torad(H))*sin(deg_torad(lat))-tan(delta)*cos(deg_torad(lat)));
	if(aziS<=0)
	{
		aziN=aziS+deg_torad(180);
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}
	

	return F;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void nightly_mode() //Nightly mode Canvas generator
{
    int year=2019;
	int month[3];
	month[0]=3;
	month[1]=4;
	month[2]=5;
	int day[31];
	for(int a=0;a<31;a++)
	{  
	  day[a]=a+1;
	}

	
	string y="2019";
	string m[3];
	m[0]="03";
	m[1]="04";
	m[2]="05";
	string d[31];
    ostringstream vchange;
	for(int a=0;a<31;a++)
	{
		vchange=NULL;
		vchange<<day[a];
		d[a] = vchange.str(); 
	}
	
	

	
	gSystem->Exec("mkdir EcoacousterMoon_Data");
	gSystem->Exec("mkdir EcoacousterMoon_Motion");
	gSystem->Exec("mkdir EcoacousterMoon_PolarSurvey");
	gSystem->Exec("mkdir EcoacousterMoon_RF");
	
	string data_arranged="NightlyMode_";
    string base_canvas1="NightlyMotion_";
	string base_canvas2="MoonSurvey_";
	
	string ref[3]={"Mar","Apr","May"};


    for(int box=0;box<3;box++)
	{

        int t=box;
		

    for(int s=0;s<31;s++)
	{
		
		double azimutal[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    double altitud[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	    ofstream daily;  //Data
	    daily.open(("C://root_v5.34.38/EcoacousterMoon_Data/"+m[t]+"-"+ref[t]+"-"+data_arranged+"-"+d[s]+"-"+m[t]+"-"+y+".txt").c_str());
	    daily<<"****Data saved in degrees****"<<endl;
	    daily<<"************HDQS*************"<<endl;
	    daily<<"Date: "<<day[s]<<"/"<<month[t]<<"/"<<year<<endl;
	    daily<<"Hour"<<"	"<<"Azimuth"<<"	"<<"Height"<<endl;
	
	    for(int u=0;u<24;u++)
	    {
		    altitud[u]=height_moon(year,month[t],day[s],u);
		    azimutal[u]=azimuth_moon(year,month[t],day[s],u);
		    daily<<u<<"	"<<azimutal[u]<<"	"<<altitud[u]<<endl;
	    }
	
	    daily.close();
	

	    c1= new TCanvas("c1","Nightly mode",1300,800);  //Motion
	    c1->SetTitle("Daily motion of the Moon");
        c1->Divide(3,1);

	
	    c1->cd(1);
	    a_ho= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    a_ho->SetPoint(u,u,azimutal[u]);
	    }
	    a_ho->SetTitle("Azimuth vs Time;Time(hours);Azimuth from the North(degrees)");
	    a_ho->SetLineColor(kRed);
	    a_ho->SetMarkerStyle(8);
	    a_ho->SetMarkerColor(kGreen+2);
	    a_ho->SetLineWidth(3);
	    a_ho->Draw("ACP");
	    a_ho->GetXaxis()->SetTitle("Time(hours)");
        a_ho->GetYaxis()->SetTitle("Azimuth from the North(degrees)");
        a_ho->GetXaxis()->CenterTitle();
        a_ho->GetYaxis()->CenterTitle();
	    a_ho->GetXaxis()->SetLabelSize(0.025);
	    a_ho->GetYaxis()->SetLabelSize(0.025);
	    a_ho->GetXaxis()->SetTitleOffset(1);
	    a_ho->GetYaxis()->SetTitleOffset(1.5);
		
		
	
	    c1->cd(2);
	    h_ho= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    h_ho->SetPoint(u,u,altitud[u]);
	    }
	    h_ho->SetTitle("Height vs Time;Time(hours);Height from the Horizon(degrees)");
	    h_ho->Draw("ACP");
	    h_ho->SetLineWidth(3);
	    h_ho->SetLineColor(kRed);
	    h_ho->SetMarkerStyle(8);
	    h_ho->SetMarkerColor(kGreen+2);
	    h_ho->GetXaxis()->SetTitle("Time(hours)");
        h_ho->GetYaxis()->SetTitle("Height from the Horizon(degrees)");
        h_ho->GetXaxis()->CenterTitle();
        h_ho->GetYaxis()->CenterTitle();
	    h_ho->GetXaxis()->SetLabelSize(0.025);
	    h_ho->GetYaxis()->SetLabelSize(0.025);
	    h_ho->GetXaxis()->SetTitleOffset(1);
	    h_ho->GetYaxis()->SetTitleOffset(1.5);
	
	    c1->cd(3);
	    a_h= new TGraph();
	    for(int u=0;u<24;u++)
	    {
		    a_h->SetPoint(u,altitud[u],azimutal[u]);
	    }
	    a_h->SetTitle("Azimuth vs Height;Height(degrees);Azimuth(degrees)");
	    a_h->SetLineColor(kRed);
	    a_h->SetMarkerStyle(8);
        a_h->SetMarkerColor(kGreen+2);
	    a_h->SetLineWidth(3);
	    a_h->Draw("ACP");
	    a_h->GetXaxis()->SetTitle("Height(degrees)");
        a_h->GetYaxis()->SetTitle("Azimuth(degrees)");
        a_h->GetXaxis()->CenterTitle();
        a_h->GetYaxis()->CenterTitle();
	    a_h->GetXaxis()->SetLabelSize(0.025);
	    a_h->GetYaxis()->SetLabelSize(0.025);
	    a_h->GetXaxis()->SetTitleOffset(1);
	    a_h->GetYaxis()->SetTitleOffset(1.5);
		


	
   	    c1->SaveAs(("C://root_v5.34.38/EcoacousterMoon_Motion/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		c1->SaveAs(("C://root_v5.34.38/EcoacousterMoon_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas1+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());




	    c2= new TCanvas("c2","Nightly Polar mode",1300,750);  //Survey
	
	    c2->cd(1);
        PL= new TGraphPolar();
        PL->SetTitle("Moon Survey through the Sky");
	    for(int u=0;u<24;u++)
	    {
		    PL->SetPoint(u,deg_torad(360-azimutal[u]),altitud[u]);
	    }  
        PL->SetMarkerStyle(20);
        PL->SetMarkerSize(1.);
        PL->SetMarkerColor(4);
        PL->SetLineColor(4);
        PL->Draw("ACP");
		// Update, otherwise GetPolargram returns 0 and the graph get wrong.
		c2->Update();
        PL->GetPolargram()->SetToDegree();
		PL->GetPolargram()->SetPolarOffset(0.01);
		PL->GetPolargram()->SetPolarLabelSize(0.03);
		
		
		
		//Textual conventions
		c2->cd(2); 
        SubT = new TText(0.8,1,("Date: "+d[s]+"/"+m[t]+"/"+y).c_str());
        SubT->SetTextAlign(22);
        SubT->SetTextFont(43);
        SubT->SetTextSize(20);
        SubT->Draw();
		
		c2->cd(3);
        CPN = new TText(1.07,0,"N");
        CPN->SetTextAlign(22);
        CPN->SetTextFont(23);
        CPN->SetTextSize(20);
        CPN->Draw();
		
		c2->cd(4);
        CPE = new TText(0,-1.1,"E");
        CPE->SetTextAlign(22);
        CPE->SetTextFont(23);
        CPE->SetTextSize(20);
        CPE->Draw();
		
		c2->cd(5);
        CPS = new TText(-1.12,0,"S");
        CPS->SetTextAlign(22);
        CPS->SetTextFont(23);
        CPS->SetTextSize(20);
        CPS->Draw();
		
		c2->cd(6);
        CPW = new TText(0.05,1.07,"W");
        CPW->SetTextAlign(22);
        CPW->SetTextFont(23);
        CPW->SetTextSize(20);
        CPW->Draw();
		
		
		
	
	    c2->SaveAs(("C://root_v5.34.38/EcoacousterMoon_PolarSurvey/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".png").c_str());
		c2->SaveAs(("C://root_v5.34.38/EcoacousterMoon_RF/"+m[t]+"-"+ref[t]+"-"+base_canvas2+"-"+d[s]+"-"+m[t]+"-"+y+".root").c_str());
		
		
	}
	}
	
	
}

void phase_mode() //Phase mode Canvas generator
{
    int year=2019;
	int hour=0;
	int fase[730]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
	for(int a=0;a<365;a++) //Mostrar
	{
		cout<<fase[a]<<" "<<fase[a+365]<<endl;
	}
	luna = new TGraph();
	for(int a=0;a<365;a++) //Llenar
	{
		
		luna->SetPoint(a,a,moon_phase(year,fase[a],fase[a+365],hour));
		
	}
	
    luna->SetMarkerStyle(7);
	luna->SetLineColor(kBlue-6);
	luna->SetLineWidth(1);
	//luna->SetLineStyle();
	luna->Draw("ACP");
	
	
}



void stat_testMOON() //Sampling for statistical test comparision
{
	int y=2019;
	int m[240];
	int d[240];
	double h[240];
	//double min[240];
	
	ifstream sample1;
	sample1.open("Month(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample1>>m[i];
	}
	sample1.close();
	
	ifstream sample2;
	sample2.open("Day(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample2>>d[i];
	}
	sample2.close();
	
	ifstream sample3;
	sample3.open("Hour(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample3>>h[i];
	}
	sample3.close();
	
	/*
	ifstream sample4;
	sample4.open("Minute(ST).txt");
	for(int i=0;i<240;i++)
	{
		sample4>>min[i];
	}
	sample4.close();
	
	for(int a=0;a<240;a++)
	{
		h[a]=h[a]+(min[a]/60);
	}
	*/
	
	ofstream TEST1;  //Data output
	TEST1.open("Statistical test Moon Azimuth.txt");
    for(int i=0;i<240;i++)
	{
		TEST1<<azimuth_moon(y,m[i],d[i],h[i])<<endl;
	}
	TEST1.close();
	
	ofstream TEST2;  //Data output
	TEST2.open("Statistical test Moon Height.txt");
    for(int i=0;i<240;i++)
	{
		TEST2<<height_moon(y,m[i],d[i],h[i])<<endl;;
	}
	TEST2.close();
	
	ofstream TEST3;  //Data output
	TEST3.open("Statistical test Moon Phase.txt");
    for(int i=0;i<240;i++)
	{
		TEST3<<moon_phase(y,m[i],d[i],h[i])<<endl;;
	}
	TEST3.close();
	
	
	
	
	cout<<"DATA SAVED!!!"<<endl;
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double juliano(int a, int m, int d,double h) //This function calculates the true Julian Date Number from 1901 to 2099
{
	double JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)-1.5;   // At noon UTC on 1 January 2000, at Julian Day 2451545. Check http://www.onlineconversion.com/julian_date.htm 
	return JD+(h/24);

}

void practical()
{
	
	
	int day[730]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};

	ofstream fase;  //Data
	ofstream moonhei;
	fase.open("Moon_phase.txt");
	moonhei.open("Moon_yearly_height.txt");
	
	for(int a=0;a<365;a++) //Llenar
	{
	    for(int u=0;u<24;u++)
		{
			
			//daily<<u<<"	"<<azimuth_moon(2019,day[a],day[a+365],u)<<"	"<<height_moon(2019,day[a],day[a+365],u)<<endl;
            fase<<u<<"	"<<moon_phase(2019,day[a],day[a+365],u)<<endl;
			moonhei<<u<<"	"<<height_moon(2019,day[a],day[a+365],u)<<endl;
		}
	}
	
    fase.close();
	moonhei.close();
    
	int p[2]={-1,-1};
	
	c1 = new TCanvas("Moon behaviour","Moon behaviour",1000,1000);
	c1->Divide(2,1);
	c1->cd(1);
    fase_lunar = new TGraph();
	for(int a=0;a<365;a++) //Llenar
	{
		for(int u=0;u<24;u++)
		{
			p[0]++;
			fase_lunar->SetPoint(p[0],p[0],moon_phase(2019,day[a],day[a+365],u));
		}
		
	}
	
    //fase_lunar->SetMarkerStyle(7);
	fase_lunar->SetLineColor(kBlue-6);
	fase_lunar->SetLineWidth(1);
	//luna->SetLineStyle();
	fase_lunar->Draw("AC");
	c1->cd(2);
	altitud_lunar = new TGraph();
	double tmp[2]={0,0};
	for(int a=0;a<365;a++) //Llenar
	{
		for(int u=0;u<24;u++)
		{
			p[1]++;
			tmp[0]=0;
			tmp[1]=0;
			tmp[0]=moon_phase(2019,day[a],day[a+365],u);
			tmp[1]=height_moon(2019,day[a],day[a+365],u);
			if((tmp[0]>=0.1)&&(tmp[1]>=0)&&((u<=6)||(u>=18)))
			{
				tmp[1]=tmp[1];
			}
			else
			{
				tmp[1]=0;
			}
		    altitud_lunar->SetPoint(p[1],p[1],tmp[1]);

		}
	}
	
    //altitud_lunar->SetMarkerStyle(7);
	altitud_lunar->SetLineColor(kBlue-6);
	altitud_lunar->SetLineWidth(1);
	//luna->SetLineStyle();
	altitud_lunar->Draw("AC");




}

void us()
{
	double data[513];
	
	for(int a=0;a<24;a++)
	{
	string d;
	ostringstream vchange;
	vchange=NULL;
	vchange<<a;
	d= vchange.str(); 
	ifstream abrir;
	abrir.open(("290319_"+d+".txt").c_str());
	for(int i=0;i<513;i++)
	{
		abrir>>data[i];
	}
	abrir.close();
	plot = new TH1D(("H"+d).c_str(),"Histogram",513,0,513);
	for(int i=0;i<513;i++)
	{
		plot->SetBinContent(i,data[i]);
	}
	plot->SetMinimum(20);
	plot->SetMaximum(100);
	cout<<plot->GetRMS()<<endl;
	//plot->Draw("C");
	plot->Draw("same");
	
    }
}

void stat_comparision()
{
	int date[]={2019,3,11};
	
	double azimuthSUN[24]; // Value Azimuth Sun
	double heightSUN[24]; // Value Height Sun
	double azimuthMOON[24]; // Value Azimuth Moon
	double heightMOON[24]; // Value Height Moon
	
	double AeAS[24]; // Average Sun Azimuth Error
	double AeHS[24]; // Average Sun Height Error
	double AeAM[24]; // Average Moon Azimuth Error
	double AeHM[24]; // Average Moon Height Error
	
	ifstream SunAziE;
	ifstream SunHeiE;
	ifstream MoonAziE;
	ifstream MoonHeiE;
	SunAziE.open("stat_AeAS.txt");
	SunHeiE.open("stat_AeHS.txt");
	MoonAziE.open("stat_AeAM.txt");
	MoonHeiE.open("stat_AeHM.txt");
	
	for(int i=0;i<24;i++)
	{
		SunAziE>>AeAS[i];
	    SunHeiE>>AeHS[i];
	    MoonAziE>>AeAM[i];
	    MoonHeiE>>AeHM[i];
	}
	SunAziE.close();
	SunHeiE.close();
	MoonAziE.close();
	MoonHeiE.close();
	
	for(int i=0;i<24;i++)
	{
        azimuthSUN[i]=deg_torad(azimuth(date[0],date[1],date[2],i));
		heightSUN[i]=height(date[0],date[1],date[2],i);
		azimuthMOON[i]=deg_torad(360-azimuth_moon(date[0],date[1],date[2],i));
		heightMOON[i]=height_moon(date[0],date[1],date[2],i);
		AeAS[i]=deg_torad(AeAS[i]);
		AeAM[i]=deg_torad(360-AeAM[i]);
	}
	
	canvas_1 = new TCanvas("Stat SUN Canvas","Stat SUN Comparision ",1300,750);  //Survey
	
	    canvas_1->cd(1);
        SC= new TGraphPolar(24,azimuthSUN,heightSUN,AeAS,AeHS);
        SC->SetTitle("Stat Comparision SUN");
        SC->SetMarkerStyle(7);
        SC->SetMarkerSize(1.);
        SC->SetMarkerColor(4);
        SC->SetLineColor(4);
		SC->SetLineWidth(2);
        SC->Draw("PE");
		// Update, otherwise GetPolargram returns 0 and the graph get wrong.
		//canvas_1->Update();
        //SC->GetPolargram()->SetToDegree();
		//SC->GetPolargram()->SetPolarOffset(0.01);
		//SC->GetPolargram()->SetPolarLabelSize(0.03);

		//Textual conventions
		canvas_1->cd(2);
        CPN = new TText(1.07,0,"N");
        CPN->SetTextAlign(22);
        CPN->SetTextFont(23);
        CPN->SetTextSize(20);
        CPN->Draw();
		
		canvas_1->cd(3);
        CPE = new TText(0,-1.1,"E");
        CPE->SetTextAlign(22);
        CPE->SetTextFont(23);
        CPE->SetTextSize(20);
        CPE->Draw();
		
		canvas_1->cd(4);
        CPS = new TText(-1.12,0,"S");
        CPS->SetTextAlign(22);
        CPS->SetTextFont(23);
        CPS->SetTextSize(20);
        CPS->Draw();
		
		canvas_1->cd(5);
        CPW = new TText(0.05,1.07,"W");
        CPW->SetTextAlign(22);
        CPW->SetTextFont(23);
        CPW->SetTextSize(20);
        CPW->Draw();
	
}

void gavilan()
{
	double data[513];
	

	ifstream abrir;
	abrir.open("gavilan.txt");
	for(int i=0;i<513;i++)
	{
		abrir>>data[i];
	}
	abrir.close();
	plot = new TH1D("Gavilan","Espectro de potencias",513,0,513);
	for(int i=0;i<513;i++)
	{
		plot->SetBinContent(i,data[i]-86);
	}
	plot->SetMinimum(-65);
	plot->SetMaximum(-20);
    plot->GetXaxis()->SetTitle("Frecuencia");
    plot->GetYaxis()->SetTitle("Amplitud");	    
	plot->GetXaxis()->SetLabelSize(0.025);
	plot->GetYaxis()->SetLabelSize(0.025);
	//plot->Draw("C");
	plot->Draw();
	
}

void bucle_proof()
{
   double vari[]={2019,8,25,21};
   for (int i=0;i<61;i++)
   {
	   cout<<"--------------Minuto "<<i<<endl;
	   cout<<"Altitud             "<<height(vari[0],vari[1],vari[2],vari[3]+(i/60))<<endl;
	   cout<<"Azimuth             "<<azimuth(vari[0],vari[1],vari[2],vari[3]+(i/60))<<endl;
   }	   
	
	
	
	
	
}









