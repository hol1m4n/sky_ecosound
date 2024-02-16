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
		aziN=aziS+deg_torad(180); /// Es la misma mierda pero no pude sacar un mejor cond xd jjajaja
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
		aziN=aziS+deg_torad(180); /// Es la misma mierda pero no pude sacar un mejor cond xd jjajaja
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}

    return rad_todeg(aziN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////

//WARNING: Due to I'm a lazy guy you have to run all for loops with int u=0 not i, because in the other functions there are variables with similar names. Att: Hollman :)

////////



void daily_mode() //Daily mode Canvas generator
{
	
	int year=2019;
	int month[3];
	month[0]=3;
	month[1]=4;
	month[2]=5;
	int day[31];
	day[0]=1;
	day[1]=2;
	day[2]=3;
	day[3]=4;
	day[4]=5;
	day[5]=6;
	day[6]=7;
	day[7]=8;
	day[8]=9;
	day[9]=10;
	day[10]=11;
	day[11]=12;
	day[12]=13;
	day[13]=14;
	day[14]=15;
	day[15]=16;
	day[16]=17;
	day[17]=18;
	day[18]=19;
	day[19]=20;
	day[20]=21;
	day[21]=22;
	day[22]=23;
	day[23]=24;
	day[24]=25;
	day[25]=26;
	day[26]=27;
	day[27]=28;
	day[28]=29;
	day[29]=30;
	day[30]=31;

	
	
	string y="2019";
	string m[3];
	m[0]="03";
	m[1]="04";
	m[2]="05";
	string d[31];
	d[0]="1";
	d[1]="2";
	d[2]="3";
	d[3]="4";
	d[4]="5";
	d[5]="6";
	d[6]="7";
	d[7]="8";
	d[8]="9";
	d[9]="10";
	d[10]="11";
	d[11]="12";
	d[12]="13";
	d[13]="14";
	d[14]="15";
	d[15]="16";
	d[16]="17";
	d[17]="18";
	d[18]="19";
	d[19]="20";
	d[20]="21";
	d[21]="22";
	d[22]="23";
	d[23]="24";
	d[24]="25";
	d[25]="26";
	d[26]="27";
	d[27]="28";
	d[28]="29";
	d[29]="30";
	d[30]="31";
	
	
	
	string data_arranged="A.DailyMode_data_";
    string base_canvas1="B.DailyMotion_";
	string base_canvas2="C.SunSurvey_";

    for(int prueba=0;prueba<3;prueba++)
	{

         
        int t=prueba;

		 
    for(int s=0;s<31;s++)
	{
		double azimutal[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    double altitud[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	    ofstream daily;
	    daily.open((m[t]+data_arranged+"_"+d[s]+"_"+m[t]+"_"+y+".txt").c_str());
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
	
	    c1= new TCanvas("c1","Daily mode",1300,800);
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
	
   	    c1->SaveAs((m[t]+base_canvas1+"_"+d[s]+"_"+m[t]+"_"+y+".png").c_str());
		c1->SaveAs((m[t]+base_canvas1+"_"+d[s]+"_"+m[t]+"_"+y+".root").c_str());
	
	
	
	
	
	    c2= new TCanvas("c2","Daily Polar mode",1300,750);
	
	    c2->cd(1);
        PL= new TGraphPolar();
        PL->SetTitle("Sun Survey through the Earth's Atmosphere");
	    for(int u=0;u<24;u++)
	    {
		    PL->SetPoint(u,deg_torad(360-azimutal[u]),altitud[u]);
	    }  
        PL->SetMarkerStyle(20);
        PL->SetMarkerSize(1.);
        PL->SetMarkerColor(4);
        PL->SetLineColor(4);
        PL->Draw("ACP");
	
	    c2->SaveAs((m[t]+base_canvas2+"_"+d[s]+"_"+m[t]+"_"+y+".png").c_str());
		c2->SaveAs((m[t]+base_canvas2+"_"+d[s]+"_"+m[t]+"_"+y+".root").c_str());
		
		
	}
	
	
	}
	
}