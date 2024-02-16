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

double algo(double value)   // Algorithm for subtract multiples of 360º
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

double sp_algo(double value)   // Special algorithm for subtract multiples of 360º for Moon calculation cycles
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
		aziN=aziS+deg_torad(180); /// Es la misma mierda pero no pude sacar un mejor cond xd jjajaja
	}
	else
	{
		aziN=aziS+TMath::Pi();
	}

    return rad_todeg(aziN);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

double juliano(int a, int m, int d,double h) //This function calculates the true Julian Date Number from 1901 to 2099
{
	double JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)-1.5;   // At noon UTC on 1 January 2000, at Julian Day 2451545. Check http://www.onlineconversion.com/julian_date.htm 
	return JD+(h/24);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void first_plotter()
{
	double RMS[]={150.316,149.93,149.866,149.929,149.768,149.586,150.665,149.118,146.009,146.965,149.49,148.63,148.437,148.494,148.901,148.246,148.364,149.108,149.567,149.965,149.826,149.851,150.104,150.44,150.636,145.176,146.621,149.11,149.495,146.7,145.186,146.486,147.196,148.825,148.21,148.064,148.401,148.923,148.723,148.853,148.852,149.129,149.586,150.053,149.644,150.099,150.312,150.037,150.193,150.11,150.136,148.985,149.586,150.14,150.401,149.085,148.134,148.463,148.426,147.97,149.116,147.694,149.211,149.193,149.016,148.224,149.238,149.817,149.891,149.797,150.386,149.962,149.925,150.114,148.38,149.098,150.038,150.859,146.334,147.288,147.029,148.312,148.225,147.724,148.115,149.218,147.865,147.944,147.54,147.862,148.819,149.555,149.783,150.276,150.12,150.698,148.676,149.886,150.093,150.291,150.573,146.036,146.369,148.072,146.722,147.096,148.587,148.411,149.174,148.724,149.333,149.003,149.204,147.294,148.715,150.015,149.715,149.682,146.679,148.333,149.376,148.854,150.151,149.056,150.305,149.993,149.74,148.355,146.389,146.546,147.05,146.371,146.786,147.222,146.137,147.553,148.218,148.04,148.65,149.794,150.087,150.072,150.271,150.118,150.057,150.025,149.971,149.943,149.949,149.544,149.608,148.635,147.284,147.534,146.879,147.18,147.954,148.278,147.973,148.88,149.445,149.472,149.792,149.848,149.969,150.474,150.276,150.14,150.01,150.062,150.038,150.055,149.866,149.816,149.833,148.78,146.553,146.077,146.065,146.336,147.104,147.166,146.948,147.044,147.795,148.471,148.931,149.338,150.146,150.05,150.003,150.202,150.047,150.352,149.902,150.025,149.92,149.82,149.885,148.4,146.845,146.148,145.861,145.067,146.415,147.874,148.971,147.764,148.382,149.426,149.121,149.519,149.201,149.491,150.022,150.124,150.057,150.264,150.071,149.983,149.916,149.943,149.896,148.346,146.913,147.91,147.336,146.686,146.815,147.043,146.617,148.42,148.914,149.247,149.285,149.828,149.347,150.036,150.091,149.955,149.862,149.868,149.893,149.283,149.78,148.004,149.279,149.447,148.357,145.495,147.566,146.202,145.303,147.085,146.81,148.061};
	double HF[]={0.167226,0.174954,0.257291,0.21134,0.336835,0.574433,0.122438,0.233185,0.230117,0.175708,0.730269,0.20438,0.165242,0.141438,0.200775,0.134358,0.151766,0.208897,0.289499,0.159854,0.111749,0.150152,0.13071,0.10026,0.115533,0.402094,0.144474,0.153925,0.144554,0.931975,0.428778,0.395949,0.242494,0.285724,0.171362,0.161971,0.15561,0.174062,0.147232,0.157074,0.145725,0.141051,0.151989,0.131141,0.118848,0.232892,0.153543,0.24159,0.17289,0.190781,0.155618,0.445709,0.207037,0.131186,0.159031,0.212451,0.281039,0.210456,0.201716,0.208105,0.26207,0.133843,0.211418,0.163334,0.158253,0.120156,0.152379,0.111439,0.122964,0.582754,0.145519,0.199668,0.233228,0.166659,0.347929,0.828665,0.251471,0.10469,0.597323,0.238126,0.178126,0.202849,0.232119,0.167174,0.165966,0.291323,0.149355,0.172366,0.161882,0.147923,0.159372,0.11925,0.12355,0.13496,0.18266,0.101343,0.131584,0.125017,0.13979,0.110253,0.0999616,0.311062,0.423837,0.594293,0.36105,0.163554,0.169301,0.176466,0.264706,0.151546,0.175557,0.131644,0.164993,0.160368,0.144492,0.122263,0.113829,0.106226,0.158528,0.22943,0.123454,0.18877,0.208689,0.142925,0.114202,0.217139,0.128825,0.133064,0.18724,0.189113,0.145784,0.146592,0.138858,0.133163,0.165274,0.130206,0.13558,0.130604,0.137159,0.111758,0.125189,0.240253,0.133739,0.151495,0.168372,0.192527,0.198491,0.183086,0.186208,0.219232,0.317153,0.17927,0.161251,0.183008,0.154841,0.150558,0.199401,0.16172,0.123386,0.11385,0.120644,0.121001,0.117307,0.104019,0.164931,0.144739,0.15633,0.170666,0.16649,0.157725,0.173179,0.185145,0.189361,0.158507,0.13193,0.185192,0.152367,0.1464,0.142306,0.157969,0.152248,0.148683,0.130043,0.134847,0.117391,0.127305,0.131701,0.104999,0.108234,0.272753,0.315248,0.197844,0.182874,0.132635,0.150718,0.176595,0.154194,0.20217,0.150967,0.152869,0.20379,0.164812,0.29361,0.308429,0.506151,0.261089,0.196609,0.127266,0.129504,0.217849,0.137763,0.145057,0.106069,0.113685,0.109315,0.158703,0.168405,0.133942,0.185449,0.20634,0.132644,0.202279,0.200246,0.143957,0.202785,0.244415,0.190211,0.138326,0.158378,0.19518,0.137064,0.141105,0.181932,0.232006,0.196799,0.174602,0.107567,0.122803,0.112512,0.162256,0.271128,0.21978,0.413,0.51587,0.293085,0.568453,0.368339,0.256924,0.359866,1.02618,0.209081,0.155063,0.196754,0.165955,0.15465,0.176339};
	double HT[]={9.99994,9.99995,9.99988,9.99992,9.99959,9.999,9.99999,9.99988,9.99996,9.99999,9.99806,9.99991,9.99998,10,9.99994,10,9.99999,9.99992,9.99975,9.99994,9.99999,9.99996,9.99998,10,9.99999,9.99999,10,9.99998,9.99999,9.9992,9.99999,9.99991,9.99998,9.99986,10,10,10,9.99999,10,10,10,10,9.99999,10,10,9.99988,9.99997,9.9999,9.99997,9.99994,9.99997,9.99948,9.99994,9.99999,9.99997,9.99992,9.99988,9.99995,9.99995,9.99994,9.9999,10,9.99991,9.99997,9.99999,10,9.99999,9.99999,9.99999,9.99863,9.99996,9.99994,9.9999,9.99997,9.99969,9.99762,9.99985,10,9.99983,9.99996,9.99999,9.99997,9.99993,9.99999,9.99999,9.99989,10,9.99999,10,10,9.99999,10,9.99999,9.99998,9.99992,9.99999,9.99999,9.99999,9.99997,9.99999,10,9.99997,9.99987,9.99937,9.99988,10,9.99999,9.99999,9.99981,9.99999,9.99997,10,9.99998,10,10,9.99999,10,10,10,9.99994,10,9.99995,9.9999,9.99999,9.99999,9.99988,9.99999,10,9.99999,9.99999,10,10,10,10,10,10,10,10,10,10,9.99998,9.99981,9.99998,9.99998,9.99995,9.99994,9.99992,9.99995,9.99995,9.99992,9.99976,9.99998,10,9.99999,10,10,9.99996,9.99999,10,10,10,10,10,10,9.99991,9.99995,9.99994,9.99995,9.99997,9.99997,9.99996,9.99994,9.99995,9.99996,9.99998,9.99997,10,10,10,10,9.99999,9.99999,10,10,10,9.99999,9.99999,10,9.99999,9.99973,9.99968,9.99988,9.99996,9.99998,9.99996,9.99995,9.99996,9.99991,9.99997,9.99999,9.99997,9.99999,9.99997,9.99999,9.99979,9.99993,9.99994,10,10,9.99989,10,9.99998,10,9.99999,9.99999,9.99997,9.99997,9.99999,9.99996,9.99993,9.99999,9.99984,9.99986,9.99997,9.99998,9.99988,9.99997,10,9.99999,9.99998,10,10,9.99997,9.99981,9.99996,9.99995,10,9.99996,9.99999,9.99996,9.99986,9.99992,9.99923,9.99915,9.99978,9.99944,9.99957,9.99989,9.9997,9.99981,9.99996,10,10,10,10,9.99998};
	
    double SUNa[256];
	double SUNh[256];
	double MOONa[256];
	double MOONh[256];
	double MOONf[256];
	
	for(int a=0;a<256;a++)
	{
		SUNa[a]=0;
	    SUNh[a]=0;
	    MOONa[a]=0;
	    MOONh[a]=0;
	    MOONf[a]=0;
	}
	


	int year=2019;
	int month[256];
	int day[256];
	double hour[256];
	float minute=30;
	
	//Arranca a las 10 am del 28 de marzo (28/03/2019) y termina a la 1 am del 8 de abril (08/04/2019)
	//Arreglo para los meses
	for(int i=0;i<72;i++)
	{
		month[i]=3;
	}
	for(int i=72;i<256;i++)
	{
		month[i]=4;
	}
	//Arreglo para los dias
    for(int i=0;i<14;i++) 
	{
		day[i]=28;
	}
    for(int i=14;i<38;i++) 
	{
		day[i]=29;
	}
	for(int i=38;i<62;i++) 
	{
		day[i]=30;
	}
	for(int i=62;i<86;i++) 
	{
	    day[i]=31;
	}
	for(int i=86;i<110;i++) 
	{
	    day[i]=1;
	}
	for(int i=110;i<134;i++) 
	{
	    day[i]=2;
	}
	for(int i=134;i<158;i++) 
	{
	    day[i]=3;
	}
	for(int i=158;i<182;i++) 
	{
	    day[i]=4;
	}
	for(int i=182;i<206;i++) 
	{
	    day[i]=5;
	}
	for(int i=206;i<230;i++) 
	{
	    day[i]=6;
	}
	for(int i=230;i<254;i++) 
	{
	    day[i]=7;
	}
	for(int i=254;i<256;i++) 
	{
	    day[i]=8;
	}
	//Arreglo para las horas
	int p=9;
	for(int i=0;i<256;i++)
	{
		p++;
		if(p>23)
		{
		    p=0;
		}
		hour[i]=p;
	}
	//Ajuste decimal para las horas
	for(int i=0;i<256;i++)
	{
		hour[i]=hour[i]+(minute/60);
	}

	// Calculadora
	for(int i=0;i<256;i++)
	{
		SUNa[i]=azimuth(year,month[i],day[i],hour[i]);
		SUNh[i]=height(year,month[i],day[i],hour[i]);
		MOONa[i]=azimuth_moon(year,month[i],day[i],hour[i]);
		MOONh[i]=height_moon(year,month[i],day[i],hour[i]);
		MOONf[i]=moon_phase(year,month[i],day[i],hour[i]);
		cout<<SUNa[i]<<"	"<<SUNh[i]<<"	"<<MOONa[i]<<"	"<<MOONh[i]<<"	"<<MOONf[i]<<"	"<<endl;
	}
	
	float JDate[]={7025.92,7025.96,7026,7026.04,7026.08,7026.13,7026.17,7026.21,7026.25,7026.29,7026.33,7026.38,7026.42,7026.46,7026.5,7026.54,7026.58,7026.63,7026.67,7026.71,7026.75,7026.79,7026.83,7026.88,7026.92,7026.96,7027,7027.04,7027.08,7027.13,7027.17,7027.21,7027.25,7027.29,7027.33,7027.38,7027.42,7027.46,7027.5,7027.54,7027.58,7027.63,7027.67,7027.71,7027.75,7027.79,7027.83,7027.88,7027.92,7027.96,7028,7028.04,7028.08,7028.13,7028.17,7028.21,7028.25,7028.29,7028.33,7028.38,7028.42,7028.46,7028.5,7028.54,7028.58,7028.63,7028.67,7028.71,7028.75,7028.79,7028.83,7028.88,7028.92,7028.96,7029,7029.04,7029.08,7029.13,7029.17,7029.21,7029.25,7029.29,7029.33,7029.38,7029.42,7029.46,7029.5,7029.54,7029.58,7029.63,7029.67,7029.71,7029.75,7029.79,7029.83,7029.88,7029.92,7029.96,7030,7030.04,7030.08,7030.13,7030.17,7030.21,7030.25,7030.29,7030.33,7030.38,7030.42,7030.46,7030.5,7030.54,7030.58,7030.63,7030.67,7030.71,7030.75,7030.79,7030.83,7030.88,7030.92,7030.96,7031,7031.04,7031.08,7031.13,7031.17,7031.21,7031.25,7031.29,7031.33,7031.38,7031.42,7031.46,7031.5,7031.54,7031.58,7031.63,7031.67,7031.71,7031.75,7031.79,7031.83,7031.88,7031.92,7031.96,7032,7032.04,7032.08,7032.13,7032.17,7032.21,7032.25,7032.29,7032.33,7032.38,7032.42,7032.46,7032.5,7032.54,7032.58,7032.63,7032.67,7032.71,7032.75,7032.79,7032.83,7032.88,7032.92,7032.96,7033,7033.04,7033.08,7033.13,7033.17,7033.21,7033.25,7033.29,7033.33,7033.38,7033.42,7033.46,7033.5,7033.54,7033.58,7033.63,7033.67,7033.71,7033.75,7033.79,7033.83,7033.88,7033.92,7033.96,7034,7034.04,7034.08,7034.13,7034.17,7034.21,7034.25,7034.29,7034.33,7034.38,7034.42,7034.46,7034.5,7034.54,7034.58,7034.63,7034.67,7034.71,7034.75,7034.79,7034.83,7034.88,7034.92,7034.96,7035,7035.04,7035.08,7035.13,7035.17,7035.21,7035.25,7035.29,7035.33,7035.38,7035.42,7035.46,7035.5,7035.54,7035.58,7035.63,7035.67,7035.71,7035.75,7035.79,7035.83,7035.88,7035.92,7035.96,7036,7036.04,7036.08,7036.13,7036.17,7036.21,7036.25,7036.29,7036.33,7036.38,7036.42,7036.46,7036.5,7036.54};
	//for(int u=0;u<256;u++)
	//{
		//JDate[u]=juliano(year,month[u],day[u],hour[u]-0.5);
	//}
	
	
	//////////////////////RMS CANVAS///////////////////////////
	
	c1= new TCanvas("RMS_ScatSA","Scatter plot",850,850);  //Motion
	c1->SetTitle("RMS_ScatSA Scatter Plot");
    c1->Divide(1,1);
	c1->cd(1);
	SAVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		SAVRMS->SetPoint(u,SUNa[u],RMS[u]);
	}
	SAVRMS->SetTitle("RMS y Azimut solar - Muestra de 256 horas");
	SAVRMS->SetMarkerStyle(8);
	SAVRMS->SetMarkerColor(kOrange+10);
	SAVRMS->Draw("AP");
	SAVRMS->GetXaxis()->SetTitle("Azimut solar ( )");
    SAVRMS->GetYaxis()->SetTitle("RMS");
	SAVRMS->GetYaxis()->SetTitleFont(132);
	SAVRMS->GetYaxis()->SetLabelFont(132);
	SAVRMS->GetYaxis()->CenterTitle();
	SAVRMS->GetXaxis()->SetTitleFont(132);
	SAVRMS->GetXaxis()->SetLabelFont(132);
	SAVRMS->GetXaxis()->CenterTitle();
	c1->Update();
	
	
	c2= new TCanvas("RMS_ScatSH","Scatter plot",850,850);  //Motion
	c2->SetTitle("RMS_ScatSH Scatter Plot");
    c2->Divide(1,1);
	c2->cd(1);
	SHVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		SHVRMS->SetPoint(u,SUNh[u],RMS[u]);
	}
	SHVRMS->SetTitle("RMS y Altitud solar - Muestra de 256 horas");
	SHVRMS->SetMarkerStyle(8);
	SHVRMS->SetMarkerColor(kOrange+10);
	SHVRMS->Draw("AP");
	SHVRMS->GetXaxis()->SetTitle("Altitud solar ( )");
    SHVRMS->GetYaxis()->SetTitle("RMS");
	SHVRMS->GetYaxis()->SetTitleFont(132);
	SHVRMS->GetYaxis()->SetLabelFont(132);
	SHVRMS->GetYaxis()->CenterTitle();
	SHVRMS->GetXaxis()->SetTitleFont(132);
	SHVRMS->GetXaxis()->SetLabelFont(132);
	SHVRMS->GetXaxis()->CenterTitle();
	c2->Update();
	
	
	
	
	c3= new TCanvas("RMS_ScatMA","Scatter plot",850,850);  //Motion
	c3->SetTitle("RMS_ScatMA Scatter Plot");
    c3->Divide(1,1);
	c3->cd(1);
	MAVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		MAVRMS->SetPoint(u,MOONa[u],RMS[u]);
	}
	MAVRMS->SetTitle("RMS y Azimut lunar - Muestra de 256 horas");
	MAVRMS->SetMarkerStyle(8);
	MAVRMS->SetMarkerColor(kGreen-3);
	MAVRMS->Draw("AP");
	MAVRMS->GetXaxis()->SetTitle("Azimut lunar ( )");
    MAVRMS->GetYaxis()->SetTitle("RMS");
	MAVRMS->GetYaxis()->SetTitleFont(132);
	MAVRMS->GetYaxis()->SetLabelFont(132);
	MAVRMS->GetYaxis()->CenterTitle();
	MAVRMS->GetXaxis()->SetTitleFont(132);
	MAVRMS->GetXaxis()->SetLabelFont(132);
	MAVRMS->GetXaxis()->CenterTitle();
	c3->Update();
	
	
	
	
	c4= new TCanvas("RMS_ScatMH","Scatter plot",850,850);  //Motion
	c4->SetTitle("RMS_ScatMH Scatter Plot");
    c4->Divide(1,1);
	c4->cd(1);
	MHVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		MHVRMS->SetPoint(u,MOONh[u],RMS[u]);
	}
	MHVRMS->SetTitle("RMS y Altitud lunar - Muestra de 256 horas");
	MHVRMS->SetMarkerStyle(8);
	MHVRMS->SetMarkerColor(kGreen-3);
	MHVRMS->Draw("AP");
	MHVRMS->GetXaxis()->SetTitle("Altitud lunar ( )");
    MHVRMS->GetYaxis()->SetTitle("RMS");
	MHVRMS->GetYaxis()->SetTitleFont(132);
	MHVRMS->GetYaxis()->SetLabelFont(132);
	MHVRMS->GetYaxis()->CenterTitle();
	MHVRMS->GetXaxis()->SetTitleFont(132);
	MHVRMS->GetXaxis()->SetLabelFont(132);
	MHVRMS->GetXaxis()->CenterTitle();
	c4->Update();
	
	
	
	c5= new TCanvas("RMS_ScatMF","Scatter plot",850,850);  //Motion
	c5->SetTitle("RMS_ScatMF Scatter Plot");
    c5->Divide(1,1);
	c5->cd(1);
	MFVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		MFVRMS->SetPoint(u,MOONf[u]*100,RMS[u]);
	}
	MFVRMS->SetTitle("RMS y Fase lunar - Muestra de 256 horas");
	MFVRMS->SetMarkerStyle(8);
	MFVRMS->SetMarkerColor(kGreen-3);
	MFVRMS->Draw("AP");
	MFVRMS->GetXaxis()->SetTitle("Fase lunar (%)");
    MFVRMS->GetYaxis()->SetTitle("RMS");
	MFVRMS->GetYaxis()->SetTitleFont(132);
	MFVRMS->GetYaxis()->SetLabelFont(132);
	MFVRMS->GetYaxis()->CenterTitle();
	MFVRMS->GetXaxis()->SetTitleFont(132);
	MFVRMS->GetXaxis()->SetLabelFont(132);
	MFVRMS->GetXaxis()->CenterTitle();
	c5->Update();
	
	
	
	
	
	
	//////////////////////HF CANVAS///////////////////////////
	c6= new TCanvas("HF_ScatSA","Scatter plot",850,850);  //Motion
	c6->SetTitle("HF_ScatSA Scatter Plot");
    c6->Divide(1,1);
	c6->cd(1);
	SAVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		SAVHF->SetPoint(u,SUNa[u],HF[u]);
	}
	SAVHF->SetTitle("Entropia frecuencial (Hf) y Azimut solar - Muestra de 256 horas");
	SAVHF->SetMarkerStyle(8);
	SAVHF->SetMarkerColor(kMagenta+2);
	SAVHF->Draw("AP");
	SAVHF->GetXaxis()->SetTitle("Azimut solar ( )");
    SAVHF->GetYaxis()->SetTitle("Hf");
	SAVHF->GetYaxis()->SetTitleFont(132);
	SAVHF->GetYaxis()->SetLabelFont(132);
	SAVHF->GetYaxis()->CenterTitle();
	SAVHF->GetXaxis()->SetTitleFont(132);
	SAVHF->GetXaxis()->SetLabelFont(132);
	SAVHF->GetXaxis()->CenterTitle();
	c6->Update();
	
	
	
	
	c7= new TCanvas("HF_ScatSH","Scatter plot",850,850);  //Motion
	c7->SetTitle("HF_ScatSH Scatter Plot");
    c7->Divide(1,1);
	c7->cd(1);
	SHVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		SHVHF->SetPoint(u,SUNh[u],HF[u]);
	}
	SHVHF->SetTitle("Entropia frecuencial (Hf) y Altitud solar - Muestra de 256 horas");
	SHVHF->SetMarkerStyle(8);
	SHVHF->SetMarkerColor(kMagenta+2);
	SHVHF->Draw("AP");
	SHVHF->GetXaxis()->SetTitle("Altitud solar ( )");
    SHVHF->GetYaxis()->SetTitle("Hf");
	SHVHF->GetYaxis()->SetTitleFont(132);
	SHVHF->GetYaxis()->SetLabelFont(132);
	SHVHF->GetYaxis()->CenterTitle();
	SHVHF->GetXaxis()->SetTitleFont(132);
	SHVHF->GetXaxis()->SetLabelFont(132);
	SHVHF->GetXaxis()->CenterTitle();
	c7->Update();
	
	
	
	
	
	
	c8= new TCanvas("HF_ScatMA","Scatter plot",850,850);  //Motion
	c8->SetTitle("HF_ScatMA Scatter Plot");
    c8->Divide(1,1);
	c8->cd(1);
	MAVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		MAVHF->SetPoint(u,MOONa[u],HF[u]);
	}
	MAVHF->SetTitle("Entropia frecuencial (Hf) y Azimut lunar - Muestra de 256 horas");
	MAVHF->SetMarkerStyle(8);
	MAVHF->SetMarkerColor(kCyan+3);
	MAVHF->Draw("AP");
	MAVHF->GetXaxis()->SetTitle("Azimut lunar ( )");
    MAVHF->GetYaxis()->SetTitle("Hf");
	MAVHF->GetYaxis()->SetTitleFont(132);
	MAVHF->GetYaxis()->SetLabelFont(132);
	MAVHF->GetYaxis()->CenterTitle();
	MAVHF->GetXaxis()->SetTitleFont(132);
	MAVHF->GetXaxis()->SetLabelFont(132);
	MAVHF->GetXaxis()->CenterTitle();
    c8->Update();
	
	
	
	
    c9= new TCanvas("HF_ScatMH","Scatter plot",850,850);  //Motion
	c9->SetTitle("HF_ScatMH Scatter Plot");
    c9->Divide(1,1);
	c9->cd(1);
	MHVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		MHVHF->SetPoint(u,MOONh[u],HF[u]);
	}
	MHVHF->SetTitle("Entropia frecuencial (Hf) y Altitud lunar - Muestra de 256 horas");
	MHVHF->SetMarkerStyle(8);
	MHVHF->SetMarkerColor(kCyan+3);
	MHVHF->Draw("AP");
	MHVHF->GetXaxis()->SetTitle("Altitud lunar ( )");
    MHVHF->GetYaxis()->SetTitle("Hf");
	MHVHF->GetYaxis()->SetTitleFont(132);
	MHVHF->GetYaxis()->SetLabelFont(132);
	MHVHF->GetYaxis()->CenterTitle();
	MHVHF->GetXaxis()->SetTitleFont(132);
	MHVHF->GetXaxis()->SetLabelFont(132);
	MHVHF->GetXaxis()->CenterTitle();

	c9->Update();
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	c10= new TCanvas("HF_ScatMF","Scatter plot",850,850);  //Motion
	c10->SetTitle("HF_ScatMF Scatter Plot");
    c10->Divide(1,1);
	c10->cd(1);
	MFVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		MFVHF->SetPoint(u,MOONf[u]*100,HF[u]);
	}
	MFVHF->SetTitle("Entropia frecuencial (Hf) y Fase lunar - Muestra de 256 horas");
	//MFVHF->SetTitleFont(132);
	MFVHF->GetXaxis()->SetTitle("Fase lunar (%)");
	MFVHF->GetXaxis()->SetTitleFont(132);
	MFVHF->GetXaxis()->SetLabelFont(132);
	MFVHF->GetXaxis()->CenterTitle();
    
	MFVHF->GetYaxis()->SetTitle("Hf");
	MFVHF->GetYaxis()->SetTitleFont(132);
	MFVHF->GetYaxis()->SetLabelFont(132);
	MFVHF->GetYaxis()->CenterTitle();
	
    MFVHF->SetMarkerStyle(8);
	MFVHF->SetMarkerColor(kCyan+3);
	MFVHF->Draw("AP");
	
	c10->Update();
	
    /*
	///////////////////////////////JULIAN//////////////////////////////////////////
	c11= new TCanvas("J_ScatRMS","Scatter plot",850,850);  //Motion
	c11->SetTitle("J_ScatRMS Scatter Plot");
    c11->Divide(1,1);
	c11->cd(1);
	JVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		JVRMS->SetPoint(u,JDate[u],RMS[u]);
	}
	JVRMS->SetTitle("RMS vs Julian date");
	JVRMS->SetMaximum(152);
	JVRMS->SetMinimum(144);
	//JVRMS->SetMarkerStyle(6);
	JVRMS->Draw("AC");
	JVRMS->GetXaxis()->SetTitle("Julian number");
    JVRMS->GetYaxis()->SetTitle("RMS");
    JVRMS->GetXaxis()->CenterTitle();
    JVRMS->GetYaxis()->CenterTitle();
	

	c12= new TCanvas("J_ScatHF","Scatter plot",850,850);  //Motion
	c12->SetTitle("J_ScatHF Scatter Plot");
    c12->Divide(1,1);
	c12->cd(1);
	JVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		JVHF->SetPoint(u,JDate[u],HF[u]);
	}
	JVHF->SetTitle("HF vs Julian date");
	JVHF->SetMarkerStyle(6);
	JVHF->Draw("AC");
	JVHF->GetXaxis()->SetTitle("Julian number");
    JVHF->GetYaxis()->SetTitle("HF");
    JVHF->GetXaxis()->CenterTitle();
    JVHF->GetYaxis()->CenterTitle();
	////////////////////////////3D/////////////////////////////////////
	c13= new TCanvas("TRI_RMS_SUN","RMS and SUN",850,850);
	c13->Divide(1,1);
	c13->cd(1);
	SaRMS= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		SaRMS->SetPoint(a,SUNa[a],SUNh[a],RMS[a]);
	}
	SaRMS->Draw("pcol");
	SaRMS->SetMarkerStyle(20);
	SaRMS->SetTitle("SUN and RMS");
	SaRMS->GetXaxis()->SetTitle("Sun Azimuth");
    SaRMS->GetYaxis()->SetTitle("Sun Height");
	SaRMS->GetZaxis()->SetTitle("RMS");
    SaRMS->GetXaxis()->CenterTitle();
    SaRMS->GetYaxis()->CenterTitle();
	SaRMS->GetZaxis()->CenterTitle();
	SaRMS->GetXaxis()->SetTitleOffset(1.5);
	SaRMS->GetYaxis()->SetTitleOffset(1.5);
	SaRMS->GetZaxis()->SetTitleOffset(1.5);
	
	c14= new TCanvas("TRI_HF_SUN","HF and SUN",850,850);
	c14->Divide(1,1);
	c14->cd(1);
	SaHF= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		SaHF->SetPoint(a,SUNa[a],SUNh[a],HF[a]);
	}
	SaHF->Draw("pcol");
	SaHF->SetMarkerStyle(20);
	SaHF->SetTitle("SUN and HF");
	SaHF->GetXaxis()->SetTitle("Sun Azimuth");
    SaHF->GetYaxis()->SetTitle("Sun Height");
	SaHF->GetZaxis()->SetTitle("HF");
    SaHF->GetXaxis()->CenterTitle();
    SaHF->GetYaxis()->CenterTitle();
	SaHF->GetZaxis()->CenterTitle();
	SaHF->GetXaxis()->SetTitleOffset(1.5);
	SaHF->GetYaxis()->SetTitleOffset(1.5);
	SaHF->GetZaxis()->SetTitleOffset(1.5);
	
	c15= new TCanvas("TRI_RMS_MOON","RMS and MOON",850,850);
	c15->Divide(1,1);
	c15->cd(1);
	MaRMS= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		MaRMS->SetPoint(a,MOONa[a],MOONh[a],RMS[a]);
	}
	MaRMS->Draw("pcol");
	MaRMS->SetMarkerStyle(20);
	MaRMS->SetTitle("MOON and RMS");
	MaRMS->GetXaxis()->SetTitle("Moon Azimuth");
    MaRMS->GetYaxis()->SetTitle("Moon Height");
	MaRMS->GetZaxis()->SetTitle("RMS");
    MaRMS->GetXaxis()->CenterTitle();
    MaRMS->GetYaxis()->CenterTitle();
	MaRMS->GetZaxis()->CenterTitle();
	MaRMS->GetXaxis()->SetTitleOffset(1.5);
	MaRMS->GetYaxis()->SetTitleOffset(1.5);
	MaRMS->GetZaxis()->SetTitleOffset(1.5);
	
	c16= new TCanvas("TRI_HF_MOON","HF and MOON",850,850);
	c16->Divide(1,1);
	c16->cd(1);
	MaHF= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		MaHF->SetPoint(a,MOONa[a],MOONh[a],HF[a]);
	}
	MaHF->Draw("pcol");
	MaHF->SetMarkerStyle(20);
	MaHF->SetTitle("MOON and HF");
	MaHF->GetXaxis()->SetTitle("Moon Azimuth");
    MaHF->GetYaxis()->SetTitle("Moon Height");
	MaHF->GetZaxis()->SetTitle("HF");
    MaHF->GetXaxis()->CenterTitle();
    MaHF->GetYaxis()->CenterTitle();
	MaHF->GetZaxis()->CenterTitle();
	MaHF->GetXaxis()->SetTitleOffset(1.5);
	MaHF->GetYaxis()->SetTitleOffset(1.5);
	MaHF->GetZaxis()->SetTitleOffset(1.5);
	
	*/
	//////////////////////////////////////////////OTROS//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	c17= new TCanvas("RMS_Scat_HF","Scatter plot",850,850);  //Motion
	c17->SetTitle("RMS_Scat_HF Scatter Plot");
    c17->Divide(1,1);
	c17->cd(1);
	RMSVHF= new TGraph();
	for(int u=0;u<256;u++)
	{
		RMSVHF->SetPoint(u,RMS[u],HF[u]);
	}
	RMSVHF->SetTitle("RMS y Entropia frecuencial (Hf)");
	RMSVHF->SetMarkerStyle(8);
	RMSVHF->SetMarkerColor(kBlack);
	RMSVHF->Draw("AP");
	RMSVHF->GetXaxis()->SetTitle("RMS");
    RMSVHF->GetYaxis()->SetTitle("Hf");
	RMSVHF->GetYaxis()->SetTitleFont(132);
	RMSVHF->GetYaxis()->SetLabelFont(132);
	RMSVHF->GetYaxis()->CenterTitle();
	RMSVHF->GetXaxis()->SetTitleFont(132);
	RMSVHF->GetXaxis()->SetLabelFont(132);
	RMSVHF->GetXaxis()->CenterTitle();
	c17->Update();
	/*
	
	c18= new TCanvas("TRI_HF_RMS_JUL","RMS,HF and JUL",850,850);
	c18->Divide(1,1);
	c18->cd(1);
	RHJ= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		RHJ->SetPoint(a,JDate[a],HF[a],RMS[a]);
	}
	RHJ->SetTitle("RMS HF JUL");
	RHJ->GetXaxis()->SetTitle("Julian Number");
    RHJ->GetYaxis()->SetTitle("RMS");
	RHJ->GetZaxis()->SetTitle("HF");
    //RHJ->GetXaxis()->CenterTitle();
    //RHJ->GetYaxis()->CenterTitle();
	//RHJ->GetZaxis()->CenterTitle();
	RHJ->Draw("cont4z");
	c18->Update();

	//RHJ->GetXaxis()->SetTitleOffset(1);
	//RHJ->GetYaxis()->SetTitleOffset(1);
	//RHJ->GetZaxis()->SetTitleOffset(1);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	c19 = new TCanvas("TRI_RMS_Heights","Heights and RMS",850,850);
	c19->Divide(1,1);
	c19->cd(1);
	RMSvHSM= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		RMSvHSM->SetPoint(a,SUNh[a],MOONh[a],RMS[a]);
	}
	RMSvHSM->Draw("pcol");
	RMSvHSM->SetMarkerStyle(20);
	RMSvHSM->SetTitle("TRI_RMS_Heights");
	RMSvHSM->GetXaxis()->SetTitle("Sun Height");
    RMSvHSM->GetYaxis()->SetTitle("Moon Height");
	RMSvHSM->GetZaxis()->SetTitle("RMS");
	//RMSvHSM->SetMaximum(152);
	//RMSvHSM->SetMinimum(144);
    RMSvHSM->GetXaxis()->CenterTitle();
    RMSvHSM->GetYaxis()->CenterTitle();
	RMSvHSM->GetZaxis()->CenterTitle();
	RMSvHSM->GetXaxis()->SetTitleOffset(1.5);
	RMSvHSM->GetYaxis()->SetTitleOffset(1.8);
	RMSvHSM->GetZaxis()->SetTitleOffset(1.8);
	c19->Update();
	
	c20 = new TCanvas("TRI_HF_Heights","Heights and HF",850,850);
	c20->Divide(1,1);
	c20->cd(1);
	HFvHSM= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		HFvHSM->SetPoint(a,SUNh[a],MOONh[a],HF[a]);
	}
	HFvHSM->Draw("pcol");
	HFvHSM->SetMarkerStyle(20);
	HFvHSM->SetTitle("TRI_HF_Heights");
	HFvHSM->GetXaxis()->SetTitle("Sun Height");
    HFvHSM->GetYaxis()->SetTitle("Moon Height");
	HFvHSM->GetZaxis()->SetTitle("HF");
    HFvHSM->GetXaxis()->CenterTitle();
    HFvHSM->GetYaxis()->CenterTitle();
	HFvHSM->GetZaxis()->CenterTitle();
	HFvHSM->GetXaxis()->SetTitleOffset(1.5);
	HFvHSM->GetYaxis()->SetTitleOffset(1.8);
	HFvHSM->GetZaxis()->SetTitleOffset(1.8);
	
	c21 = new TCanvas("TRI_RMS_Azimuths","Azimuths and RMS",850,850);
	c21->Divide(1,1);
	c21->cd(1);
	RMSvASM= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		RMSvASM->SetPoint(a,SUNa[a],MOONa[a],RMS[a]);
	}
	RMSvASM->Draw("pcol");
	RMSvASM->SetMarkerStyle(20);
	RMSvASM->SetTitle("TRI_RMS_Azimuths");
	RMSvASM->GetXaxis()->SetTitle("Sun Azimuth");
    RMSvASM->GetYaxis()->SetTitle("Moon Azimuth");
	RMSvASM->GetZaxis()->SetTitle("RMS");
	//RMSvASM->SetMaximum(152);
	//RMSvASM->SetMinimum(144);
    RMSvASM->GetXaxis()->CenterTitle();
    RMSvASM->GetYaxis()->CenterTitle();
	RMSvASM->GetZaxis()->CenterTitle();
	RMSvASM->GetXaxis()->SetTitleOffset(1.5);
	RMSvASM->GetYaxis()->SetTitleOffset(1.8);
	RMSvASM->GetZaxis()->SetTitleOffset(1.8);
	c21->Update();
	
	c20 = new TCanvas("TRI_HF_Azimuths","Azimuths and HF",850,850);
	c20->Divide(1,1);
	c20->cd(1);
	HFvASM= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		HFvASM->SetPoint(a,SUNa[a],MOONa[a],HF[a]);
	}
	HFvASM->Draw("pcol");
	HFvASM->SetMarkerStyle(20);
	HFvASM->SetTitle("TRI_HF_Azimuths");
	HFvASM->GetXaxis()->SetTitle("Sun Azimuth");
    HFvASM->GetYaxis()->SetTitle("Moon Azimuth");
	HFvASM->GetZaxis()->SetTitle("HF");
    HFvASM->GetXaxis()->CenterTitle();
    HFvASM->GetYaxis()->CenterTitle();
	HFvASM->GetZaxis()->CenterTitle();
	HFvASM->GetXaxis()->SetTitleOffset(1.5);
	HFvASM->GetYaxis()->SetTitleOffset(1.8);
	HFvASM->GetZaxis()->SetTitleOffset(1.8);
	
	*/
	
	gSystem->Exec("mkdir First_Scatter_Plots");
	
	string name[18];
	name[0]="RMS_ScatSA";
	name[1]="RMS_ScatSH";
	name[2]="RMS_ScatMA";
	name[3]="RMS_ScatMH";
	name[4]="RMS_ScatMF";
    name[5]="HF_ScatSA";
	name[6]="HF_ScatSH";
	name[7]="HF_ScatMA";
	name[8]="HF_ScatMH";
	name[9]="HF_ScatMF";
	name[10]="J_ScatRMS";
	name[11]="J_ScatHF";
	name[12]="TRI_RMS_SUN";
	name[13]="TRI_HF_SUN";
	name[14]="TRI_RMS_MOON";
	name[15]="TRI_HF_MOON";
	name[16]="RMS_Scat_HF";
	name[17]="TRI_HF_RMS_JUL";
	
	c1->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[0]+".root").c_str());
    c2->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[1]+".root").c_str());
	c3->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[2]+".root").c_str());
	c4->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[3]+".root").c_str());
	c5->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[4]+".root").c_str());
	c6->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[5]+".root").c_str());
	c7->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[6]+".root").c_str());
	c8->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[7]+".root").c_str());
	c9->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[8]+".root").c_str());
	c10->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[9]+".root").c_str());
	//c11->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[10]+".root").c_str());
	//c12->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[11]+".root").c_str());
	//c13->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[12]+".root").c_str());
	//c14->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[13]+".root").c_str());
	//c15->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[14]+".root").c_str());
	//c16->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[15]+".root").c_str());
	c17->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[16]+".root").c_str());
	//c18->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[17]+".root").c_str());
	c1->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[0]+".png").c_str());
    c2->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[1]+".png").c_str());
	c3->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[2]+".png").c_str());
	c4->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[3]+".png").c_str());
	c5->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[4]+".png").c_str());
	c6->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[5]+".png").c_str());
	c7->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[6]+".png").c_str());
	c8->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[7]+".png").c_str());
	c9->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[8]+".png").c_str());
	c10->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[9]+".png").c_str());
	//c11->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[10]+".png").c_str());
	//c12->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[11]+".png").c_str());
	//c13->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[12]+".png").c_str());
	//c14->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[13]+".png").c_str());
	//c15->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[14]+".png").c_str());
	//c16->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[15]+".png").c_str());
	c17->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[16]+".png").c_str());
	//c18->SaveAs(("C://root_v5.34.38/First_Scatter_Plots/"+name[17]+".png").c_str());
	
	
	
}

void fixed_plotter()
{
	double RMS[]={150.316,149.93,149.866,149.929,149.768,149.586,150.665,149.118,146.009,146.965,149.49,148.63,148.437,148.494,148.901,148.246,148.364,149.108,149.567,149.965,149.826,149.851,150.104,150.44,150.636,145.176,146.621,149.11,149.495,146.7,145.186,146.486,147.196,148.825,148.21,148.064,148.401,148.923,148.723,148.853,148.852,149.129,149.586,150.053,149.644,150.099,150.312,150.037,150.193,150.11,150.136,148.985,149.586,150.14,150.401,149.085,148.134,148.463,148.426,147.97,149.116,147.694,149.211,149.193,149.016,148.224,149.238,149.817,149.891,149.797,150.386,149.962,149.925,150.114,148.38,149.098,150.038,150.859,146.334,147.288,147.029,148.312,148.225,147.724,148.115,149.218,147.865,147.944,147.54,147.862,148.819,149.555,149.783,150.276,150.12,150.698,148.676,149.886,150.093,150.291,150.573,146.036,146.369,148.072,146.722,147.096,148.587,148.411,149.174,148.724,149.333,149.003,149.204,147.294,148.715,150.015,149.715,149.682,146.679,148.333,149.376,148.854,150.151,149.056,150.305,149.993,149.74,148.355,146.389,146.546,147.05,146.371,146.786,147.222,146.137,147.553,148.218,148.04,148.65,149.794,150.087,150.072,150.271,150.118,150.057,150.025,149.971,149.943,149.949,149.544,149.608,148.635,147.284,147.534,146.879,147.18,147.954,148.278,147.973,148.88,149.445,149.472,149.792,149.848,149.969,150.474,150.276,150.14,150.01,150.062,150.038,150.055,149.866,149.816,149.833,148.78,146.553,146.077,146.065,146.336,147.104,147.166,146.948,147.044,147.795,148.471,148.931,149.338,150.146,150.05,150.003,150.202,150.047,150.352,149.902,150.025,149.92,149.82,149.885,148.4,146.845,146.148,145.861,145.067,146.415,147.874,148.971,147.764,148.382,149.426,149.121,149.519,149.201,149.491,150.022,150.124,150.057,150.264,150.071,149.983,149.916,149.943,149.896,148.346,146.913,147.91,147.336,146.686,146.815,147.043,146.617,148.42,148.914,149.247,149.285,149.828,149.347,150.036,150.091,149.955,149.862,149.868,149.893,149.283,149.78,148.004,149.279,149.447,148.357,145.495,147.566,146.202,145.303,147.085,146.81,148.061};
	double HF[]={0.167226,0.174954,0.257291,0.21134,0.336835,0.574433,0.122438,0.233185,0.230117,0.175708,0.730269,0.20438,0.165242,0.141438,0.200775,0.134358,0.151766,0.208897,0.289499,0.159854,0.111749,0.150152,0.13071,0.10026,0.115533,0.402094,0.144474,0.153925,0.144554,0.931975,0.428778,0.395949,0.242494,0.285724,0.171362,0.161971,0.15561,0.174062,0.147232,0.157074,0.145725,0.141051,0.151989,0.131141,0.118848,0.232892,0.153543,0.24159,0.17289,0.190781,0.155618,0.445709,0.207037,0.131186,0.159031,0.212451,0.281039,0.210456,0.201716,0.208105,0.26207,0.133843,0.211418,0.163334,0.158253,0.120156,0.152379,0.111439,0.122964,0.582754,0.145519,0.199668,0.233228,0.166659,0.347929,0.828665,0.251471,0.10469,0.597323,0.238126,0.178126,0.202849,0.232119,0.167174,0.165966,0.291323,0.149355,0.172366,0.161882,0.147923,0.159372,0.11925,0.12355,0.13496,0.18266,0.101343,0.131584,0.125017,0.13979,0.110253,0.0999616,0.311062,0.423837,0.594293,0.36105,0.163554,0.169301,0.176466,0.264706,0.151546,0.175557,0.131644,0.164993,0.160368,0.144492,0.122263,0.113829,0.106226,0.158528,0.22943,0.123454,0.18877,0.208689,0.142925,0.114202,0.217139,0.128825,0.133064,0.18724,0.189113,0.145784,0.146592,0.138858,0.133163,0.165274,0.130206,0.13558,0.130604,0.137159,0.111758,0.125189,0.240253,0.133739,0.151495,0.168372,0.192527,0.198491,0.183086,0.186208,0.219232,0.317153,0.17927,0.161251,0.183008,0.154841,0.150558,0.199401,0.16172,0.123386,0.11385,0.120644,0.121001,0.117307,0.104019,0.164931,0.144739,0.15633,0.170666,0.16649,0.157725,0.173179,0.185145,0.189361,0.158507,0.13193,0.185192,0.152367,0.1464,0.142306,0.157969,0.152248,0.148683,0.130043,0.134847,0.117391,0.127305,0.131701,0.104999,0.108234,0.272753,0.315248,0.197844,0.182874,0.132635,0.150718,0.176595,0.154194,0.20217,0.150967,0.152869,0.20379,0.164812,0.29361,0.308429,0.506151,0.261089,0.196609,0.127266,0.129504,0.217849,0.137763,0.145057,0.106069,0.113685,0.109315,0.158703,0.168405,0.133942,0.185449,0.20634,0.132644,0.202279,0.200246,0.143957,0.202785,0.244415,0.190211,0.138326,0.158378,0.19518,0.137064,0.141105,0.181932,0.232006,0.196799,0.174602,0.107567,0.122803,0.112512,0.162256,0.271128,0.21978,0.413,0.51587,0.293085,0.568453,0.368339,0.256924,0.359866,1.02618,0.209081,0.155063,0.196754,0.165955,0.15465,0.176339};
	double HT[]={9.99994,9.99995,9.99988,9.99992,9.99959,9.999,9.99999,9.99988,9.99996,9.99999,9.99806,9.99991,9.99998,10,9.99994,10,9.99999,9.99992,9.99975,9.99994,9.99999,9.99996,9.99998,10,9.99999,9.99999,10,9.99998,9.99999,9.9992,9.99999,9.99991,9.99998,9.99986,10,10,10,9.99999,10,10,10,10,9.99999,10,10,9.99988,9.99997,9.9999,9.99997,9.99994,9.99997,9.99948,9.99994,9.99999,9.99997,9.99992,9.99988,9.99995,9.99995,9.99994,9.9999,10,9.99991,9.99997,9.99999,10,9.99999,9.99999,9.99999,9.99863,9.99996,9.99994,9.9999,9.99997,9.99969,9.99762,9.99985,10,9.99983,9.99996,9.99999,9.99997,9.99993,9.99999,9.99999,9.99989,10,9.99999,10,10,9.99999,10,9.99999,9.99998,9.99992,9.99999,9.99999,9.99999,9.99997,9.99999,10,9.99997,9.99987,9.99937,9.99988,10,9.99999,9.99999,9.99981,9.99999,9.99997,10,9.99998,10,10,9.99999,10,10,10,9.99994,10,9.99995,9.9999,9.99999,9.99999,9.99988,9.99999,10,9.99999,9.99999,10,10,10,10,10,10,10,10,10,10,9.99998,9.99981,9.99998,9.99998,9.99995,9.99994,9.99992,9.99995,9.99995,9.99992,9.99976,9.99998,10,9.99999,10,10,9.99996,9.99999,10,10,10,10,10,10,9.99991,9.99995,9.99994,9.99995,9.99997,9.99997,9.99996,9.99994,9.99995,9.99996,9.99998,9.99997,10,10,10,10,9.99999,9.99999,10,10,10,9.99999,9.99999,10,9.99999,9.99973,9.99968,9.99988,9.99996,9.99998,9.99996,9.99995,9.99996,9.99991,9.99997,9.99999,9.99997,9.99999,9.99997,9.99999,9.99979,9.99993,9.99994,10,10,9.99989,10,9.99998,10,9.99999,9.99999,9.99997,9.99997,9.99999,9.99996,9.99993,9.99999,9.99984,9.99986,9.99997,9.99998,9.99988,9.99997,10,9.99999,9.99998,10,10,9.99997,9.99981,9.99996,9.99995,10,9.99996,9.99999,9.99996,9.99986,9.99992,9.99923,9.99915,9.99978,9.99944,9.99957,9.99989,9.9997,9.99981,9.99996,10,10,10,10,9.99998};
	
    double SUNa[256];
	double SUNh[256];
	double MOONa[256];
	double MOONh[256];
	double MOONf[256];
	
	for(int a=0;a<256;a++)
	{
		SUNa[a]=0;
	    SUNh[a]=0;
	    MOONa[a]=0;
	    MOONh[a]=0;
	    MOONf[a]=0;
	}
	


	int year=2019;
	int month[256];
	int day[256];
	double hour[256];
	float minute=30;
	
	//Arranca a las 10 am del 28 de marzo (28/03/2019) y termina a la 1 am del 8 de abril (08/04/2019)
	//Arreglo para los meses
	for(int i=0;i<72;i++)
	{
		month[i]=3;
	}
	for(int i=72;i<256;i++)
	{
		month[i]=4;
	}
	//Arreglo para los dias
    for(int i=0;i<14;i++) 
	{
		day[i]=28;
	}
    for(int i=14;i<38;i++) 
	{
		day[i]=29;
	}
	for(int i=38;i<62;i++) 
	{
		day[i]=30;
	}
	for(int i=62;i<86;i++) 
	{
	    day[i]=31;
	}
	for(int i=86;i<110;i++) 
	{
	    day[i]=1;
	}
	for(int i=110;i<134;i++) 
	{
	    day[i]=2;
	}
	for(int i=134;i<158;i++) 
	{
	    day[i]=3;
	}
	for(int i=158;i<182;i++) 
	{
	    day[i]=4;
	}
	for(int i=182;i<206;i++) 
	{
	    day[i]=5;
	}
	for(int i=206;i<230;i++) 
	{
	    day[i]=6;
	}
	for(int i=230;i<254;i++) 
	{
	    day[i]=7;
	}
	for(int i=254;i<256;i++) 
	{
	    day[i]=8;
	}
	//Arreglo para las horas
	int p=9;
	for(int i=0;i<256;i++)
	{
		p++;
		if(p>23)
		{
		    p=0;
		}
		hour[i]=p;
	}
	//Ajuste decimal para las horas
	for(int i=0;i<256;i++)
	{
		hour[i]=hour[i]+(minute/60);
	}
    
	// Calculadora
	for(int i=0;i<256;i++)
	{
		SUNa[i]=azimuth(year,month[i],day[i],hour[i]);
		SUNh[i]=height(year,month[i],day[i],hour[i]);
		MOONa[i]=azimuth_moon(year,month[i],day[i],hour[i]);
		MOONh[i]=height_moon(year,month[i],day[i],hour[i]);
		MOONf[i]=moon_phase(year,month[i],day[i],hour[i]);
	}
	
	float JDate[]={7025.92,7025.96,7026,7026.04,7026.08,7026.13,7026.17,7026.21,7026.25,7026.29,7026.33,7026.38,7026.42,7026.46,7026.5,7026.54,7026.58,7026.63,7026.67,7026.71,7026.75,7026.79,7026.83,7026.88,7026.92,7026.96,7027,7027.04,7027.08,7027.13,7027.17,7027.21,7027.25,7027.29,7027.33,7027.38,7027.42,7027.46,7027.5,7027.54,7027.58,7027.63,7027.67,7027.71,7027.75,7027.79,7027.83,7027.88,7027.92,7027.96,7028,7028.04,7028.08,7028.13,7028.17,7028.21,7028.25,7028.29,7028.33,7028.38,7028.42,7028.46,7028.5,7028.54,7028.58,7028.63,7028.67,7028.71,7028.75,7028.79,7028.83,7028.88,7028.92,7028.96,7029,7029.04,7029.08,7029.13,7029.17,7029.21,7029.25,7029.29,7029.33,7029.38,7029.42,7029.46,7029.5,7029.54,7029.58,7029.63,7029.67,7029.71,7029.75,7029.79,7029.83,7029.88,7029.92,7029.96,7030,7030.04,7030.08,7030.13,7030.17,7030.21,7030.25,7030.29,7030.33,7030.38,7030.42,7030.46,7030.5,7030.54,7030.58,7030.63,7030.67,7030.71,7030.75,7030.79,7030.83,7030.88,7030.92,7030.96,7031,7031.04,7031.08,7031.13,7031.17,7031.21,7031.25,7031.29,7031.33,7031.38,7031.42,7031.46,7031.5,7031.54,7031.58,7031.63,7031.67,7031.71,7031.75,7031.79,7031.83,7031.88,7031.92,7031.96,7032,7032.04,7032.08,7032.13,7032.17,7032.21,7032.25,7032.29,7032.33,7032.38,7032.42,7032.46,7032.5,7032.54,7032.58,7032.63,7032.67,7032.71,7032.75,7032.79,7032.83,7032.88,7032.92,7032.96,7033,7033.04,7033.08,7033.13,7033.17,7033.21,7033.25,7033.29,7033.33,7033.38,7033.42,7033.46,7033.5,7033.54,7033.58,7033.63,7033.67,7033.71,7033.75,7033.79,7033.83,7033.88,7033.92,7033.96,7034,7034.04,7034.08,7034.13,7034.17,7034.21,7034.25,7034.29,7034.33,7034.38,7034.42,7034.46,7034.5,7034.54,7034.58,7034.63,7034.67,7034.71,7034.75,7034.79,7034.83,7034.88,7034.92,7034.96,7035,7035.04,7035.08,7035.13,7035.17,7035.21,7035.25,7035.29,7035.33,7035.38,7035.42,7035.46,7035.5,7035.54,7035.58,7035.63,7035.67,7035.71,7035.75,7035.79,7035.83,7035.88,7035.92,7035.96,7036,7036.04,7036.08,7036.13,7036.17,7036.21,7036.25,7036.29,7036.33,7036.38,7036.42,7036.46,7036.5,7036.54};
    
	
	double rain[]={435,8,140,205,10,0.1,0,135,0,218};
	double tmp[]={0,0,0,0,0,0,0,0};
	double promHF[]={0.240745,0.19601,0.25456,0.20713,0.156872,0.172596,0.151341,0.194776};
	double promRMS[]={147.625,148.841,148.441,148.049,148.177,148.463,148.502,147.943};
	double promJDate[]={7026.5,7027.5,7028.5,7029.5,7030.5,7031.5,7032.5,7033.5};
	double promMP[10];
    p=-1;
	for(int i=0;i<10;i++)
	{
        for(int u=0;u<24;u++)
		{
			p++;
            //tmp[0]=HF[p+14];
			//tmp[1]=tmp[0]+tmp[1];
			//tmp[2]=RMS[p+14];
			//tmp[3]=tmp[2]+tmp[3];
			//tmp[4]=JDate[p+14];
			//tmp[5]=tmp[4]+tmp[5];
			tmp[6]=MOONf[p+14];
			tmp[7]=tmp[6]+tmp[7];
			if(u==23)
			{
				//promHF[i]=tmp[1]/24;
				//promRMS[i]=tmp[3]/24;
				//promJDate[i]=tmp[5]/24;
				//promMP[i]=tmp[7]/24;
				//tmp[1]=0;
				//tmp[3]=0;
				//tmp[5]=0;
				tmp[7]=0;
			}
			
		}
	}
	
    /*
    c18= new TCanvas("TRI_RAIN_RMS_JUL","RMS,RAIN and JUL",850,850);
	c18->Divide(1,1);
	c18->cd(1);
	RHJ= new TGraph2D();
	for(int a=0;a<8;a++)
	{
		RHJ->SetPoint(a,promJDate[a],rain[a],promRMS[a]);
	}
	RHJ->SetTitle("RMS RAIN JUL");
	RHJ->GetXaxis()->SetTitle("Julian Number");
    RHJ->GetYaxis()->SetTitle("PP");
	RHJ->GetZaxis()->SetTitle("RMS");
    RHJ->GetXaxis()->CenterTitle();
    RHJ->GetYaxis()->CenterTitle();
	RHJ->GetZaxis()->CenterTitle();
	RHJ->Draw("cont4z");
    RHJ->SetMarkerStyle(20);
	
	
	c19= new TCanvas("TRI_RAIN_HF_JUL","HF,RAIN and JUL",850,850);
	c19->Divide(1,1);
	c19->cd(1);
	RRJ= new TGraph2D();
	for(int a=0;a<8;a++)
	{
		RRJ->SetPoint(a,promJDate[a],rain[a],promHF[a]);
	}
	RRJ->SetTitle("HF RAIN JUL");
	RRJ->GetXaxis()->SetTitle("Julian Number");
    RRJ->GetYaxis()->SetTitle("PP");
	RRJ->GetZaxis()->SetTitle("HF");
    RRJ->GetXaxis()->CenterTitle();
    RRJ->GetYaxis()->CenterTitle();
	RRJ->GetZaxis()->CenterTitle();
	RRJ->Draw("cont4z");
    RRJ->SetMarkerStyle(20);
	
    

	//////////////////////RMS CANVAS///////////////////////////
    
	c1= new TCanvas("RMS && PP","Scatter plot",850,850);  //Motion
	c1->SetTitle("RMS && PP");
    c1->Divide(1,1);
	c1->cd(1);
	RMSPP= new TGraph();
	for(int u=0;u<8;u++)
	{
		RMSPP->SetPoint(u,rain[u],promRMS[u]);
	}
	//MBVRMS->SetTitle("RMS vs Moon brightness");
	RMSPP->SetMarkerStyle(8);
	RMSPP->Draw("ACP");
	RMSPP->GetXaxis()->SetTitle("pp (mm) ");
    RMSPP->GetYaxis()->SetTitle("RMS");
    RMSPP->GetXaxis()->CenterTitle();
    RMSPP->GetYaxis()->CenterTitle();
	
	c2= new TCanvas("HF && PP","Scatter plot",850,850);  //Motion
	c2->SetTitle("HF && PP");
    c2->Divide(1,1);
	c2->cd(1);
	HFPP= new TGraph();
	for(int u=0;u<8;u++)
	{
		HFPP->SetPoint(u,rain[u],promHF[u]);
	}
	//MBVRMS->SetTitle("RMS vs Moon brightness");
	HFPP->SetMarkerStyle(8);
	HFPP->Draw("ACP");
	HFPP->GetXaxis()->SetTitle("pp (mm) ");
    HFPP->GetYaxis()->SetTitle("HF");
    HFPP->GetXaxis()->CenterTitle();
    HFPP->GetYaxis()->CenterTitle();
	
    c3= new TCanvas("JD && PP","Scatter plot",850,850);  //Motion
	c3->SetTitle("JD && PP");
    c3->Divide(1,1);
	c3->cd(1);
	JDPP= new TGraph();
	for(int u=0;u<8;u++)
	{
		JDPP->SetPoint(u,promJDate[u],rain[u]);
	}
	//MBVRMS->SetTitle("RMS vs Moon brightness");
	JDPP->SetMarkerStyle(8);
	JDPP->Draw("ACP");
	JDPP->GetXaxis()->SetTitle("Julian Date");
    JDPP->GetYaxis()->SetTitle("pp (mm) ");
    JDPP->GetXaxis()->CenterTitle();
    JDPP->GetYaxis()->CenterTitle();
	
	*/
	
	
	
	
	
	/*
	c11= new TCanvas("J_ScatRMS","Scatter plot",850,850);  //Motion
	c11->SetTitle("Sun H vs HF");
    c11->Divide(1,1);
	c11->cd(1);
	JVRMS= new TGraph();
	for(int u=0;u<256;u++)
	{
		JVRMS->SetPoint(u,JDate[u],HT[u]);
	}
	JVRMS->SetTitle("HF vs Sun Height");
	//JVRMS->SetLineColor(kBlue+2);
	//JVRMS->SetMarkerStyle();
	JVRMS->Draw("AC");
	f = new TF1("Fitting","[0]*(sin(x^[1])*sin(x^[2])*sin(x^[3]))+(cos(x^[1])*cos(x^[2])*cos(x^[3]))");
    JVRMS->Fit(f);
     */
	 
	 
	/*
	string name;
	name="first_fit";
	ofstream first_fit;
	first_fit.open(("C://root_v5.34.38/First_Scatter_Plots/"+name+".txt").c_str());
     */
	 
	 string name="tridimensinal";
	c21 = new TCanvas("TRI_RMS_Heights","Heights and RMS",2000,2000);
	c21->Divide(1,1);
	c21->cd(1);
	RMSvHSM= new TGraph2D();
	for(int a=0;a<256;a++)
	{
		RMSvHSM->SetPoint(a,MOONh[a],MOONf[a]*100,RMS[a]);
		//first_fit<<a<<"	"<<SUNh[a]<<"	"<<MOONh[a]<<"	"<<RMS[a]<<endl;
	}
	//first_fit.close();
	RMSvHSM->SetMarkerStyle(20);
	RMSvHSM->SetTitle("RMS vs Altitud lunar vs Fase lunar- Muestra de 256 puntos");
	RMSvHSM->GetXaxis()->SetTitle("ALTURA DE LA LUNA (o)");
    RMSvHSM->GetYaxis()->SetTitle("FASE DE LA LUNA (%)");
	RMSvHSM->GetZaxis()->SetTitle("RMS");
    RMSvHSM->GetXaxis()->CenterTitle();
    RMSvHSM->GetYaxis()->CenterTitle();
	RMSvHSM->GetZaxis()->CenterTitle();
	RMSvHSM->GetXaxis()->SetTitleOffset(2);
	RMSvHSM->GetYaxis()->SetTitleOffset(2);
	RMSvHSM->GetZaxis()->SetTitleOffset(1.3);
	
	RMSvHSM->GetYaxis()->SetTitleFont(132);
	RMSvHSM->GetYaxis()->SetLabelFont(132);
	RMSvHSM->GetYaxis()->CenterTitle();
	RMSvHSM->GetXaxis()->SetTitleFont(132);
	RMSvHSM->GetXaxis()->SetLabelFont(132);
	RMSvHSM->GetXaxis()->CenterTitle();
	RMSvHSM->GetZaxis()->SetTitleFont(132);
	RMSvHSM->GetZaxis()->SetLabelFont(132);
	RMSvHSM->GetZaxis()->CenterTitle();
	
	RMSvHSM->Draw("pcol & line");
	
	c21->Update();
	c21->SaveAs(("C://root_v5.34.38/variables/"+name+".png").c_str());
	//f = new TF2("Fitting","[0]*sin(x^[2])*sin(y^[1])");
    //RMSvHSM->Fit(f);
	
	
	
}



