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

void main()
{
	cout<<"Astronomic calculations code"<<endl;
	cout<<"--------------------------"<<endl;
	cout<<"These calculations only apply for Time zone +5 Colombia regions"<<endl;
	cout<<"--------------------------"<<endl;
    cout<<"--------------------------"<<endl;
	cout<<"--------------------------"<<endl;
	cout<<"**************************"<<endl;
	cout<<"Please type the longitude:"<<endl;
	double lon; //Longitude
	cin>>lon;
	cout<<"Please type the latitude:"<<endl;
	double lat;   //Latitude
	cin>>lat;
	cout<<"Please type the year:"<<endl;
	int a;   //Year
	cin>>a;
	cout<<"Please type the month:"<<endl;
	int m;   //Month
	cin>>m;
	cout<<"Please type the day:"<<endl;
	int d;   //Day
	cin>>d;
	cout<<"Please type the hour(time):"<<endl;
	double time;   //Hour
	cin>>time;
	
	//////////////////////////SUN part/////////////////////////////////////////////
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
	
	//////////////////////////MOON part/////////////////////////////////////////////
	double JD_m=0; //Julian Date Number from epoch 2010 January 0
	double RJD_m=0;//The true Julian Date Number
	double T=0; //Number of Julian Centuries since 2000 January 0.5
	double eg=0;//The Sun’s mean ecliptic longitude measured in degrees
	double wg=0;//The longitude of the Sun at perigee measured in degrees
	double e_xm=0;//The eccentricity of the Sun–Earth orbit
	double Ms=0;//Mean Anomaly of the Sun
	double E[2]={0,0};//The eccentric anomaly of the Sun measured in radians/degrees
	double v_m=0;//True anomaly of the Sun in degrees
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
	double lambda_m=0;//Ecliptic longitude Lambda
	double beta_m=0;//Ecliptic latitude Beta
	double epsilon_m=0;//The obliquity of the ecliptic, the angle between the planes of the equator and the ecliptic measured in degrees
	double alpha_m=0;//The right ascension of the ecuatorial plane
	double delta_m=0;//The declination of the ecuatorial plane
	double Me=0; // Mean Anomaly of the Earth
	double n_m=0;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	double Oi_m[2]={0,0}; //Sidereal time on the prime meridian + 5 plus hours for Colombia
	double theta_m[2]={0,0};//The sideral time
	double H_m=0; //Hour Angle
	double D=0; //Change in the Moon's orbit from the line of reference through the month
	double F=0; ////Moon's phase
	double hei_m=0; //The height of the Moon
    double aziS_m=0; //The Moon's Azimuth from the South
	double aziN_m=0; //The Moon's Azimuth from the North
	
	double Lm0=91.929336; //Moon’s mean longitude at the epoch  measured in degrees
	double P0=130.143076; //Moon’s mean longitude of the perigee  at the epoch measured in degrees
	double N0=291.682547; //Moon’s mean longitude of the node  at the epoch measured in degrees
	double i_m=5.145396; //Inclination of Moon’s orbit measured in degrees
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	
	int dm[3]={0,0,0};//Dummy 1
	double abc[3]={0,0,0};//Dummy 2
	double right=0;//Dummy 3
	double y=0; //Dummy 4
	double x=0; //Dummy 5
	double k=0; //Dummy 6
	double arcs_deg[2]={0,0};//Dummy 7
	double b=0; //Dummy 8
	
	JD_m=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)+(time/24);   // At midnight UTC on 31 December of 2010, at Julian Day 2455196.5. Check http://www.onlineconversion.com/julian_date.htm
	JD_m=JD_m-3653;
	dm[0]=(m+9)/12;
    dm[1]=(7*(a+dm[0]))/4;
	dm[2]=(275*m)/9;
	RJD_m=(367*a)-dm[1]+dm[2]+d+1721013.5+(time/24);
	T=(RJD_m-2451545.0)/36525.0;
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
	e_xm=0.01675104-(0.0000418*T)-(0.000000126*pow(T,2));
	if(T==1.09997262149212860) //The eccentricity of the Sun–Earth orbit at epoch 2010 measured in degrees
	{
		e_xm=0.016705;
	}
	else
	{
		e_xm=sp_algo(e_xm);
	}
	Ms=eg-wg;
	Ms=sp_algo(Ms);
	E[0]=E_f(e_xm,Ms); //Measured in radians
	E[1]=rad_todeg(E[0]);//Measured in degrees
	abc[0]=1+e_xm;
	abc[1]=1-e_xm;
	abc[2]=sqrt(abc[0]/abc[1]);
    right=abc[2]*tan(E[0]/2); // E measured in radians
	v_m=atan(right)*2;
	v_m=sp_algo(rad_todeg(v_m));
	Ls=v_m+wg;
	Ls=sp_algo(Ls);
    Lm=(13.1763966*JD_m)+Lm0;
	Lm=sp_algo(Lm);
	Mm=Lm-(0.1114041*JD_m)-P0;
	Mm=sp_algo(Mm);
	N=N0-(0.0529539*JD_m);
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
	y=cos(deg_torad(i_m))*sin(deg_torad(LmR2-NR));
	x=cos(deg_torad(LmR2-NR));
	lambda_m=atan2(y,x)+deg_torad(NR);
	if( lambda_m<0)
	{
		lambda_m=lambda_m+(2*TMath::Pi());
	}
	k=sin(deg_torad(i_m))*sin(deg_torad(LmR2-NR));
	beta_m=asin(k);
	arcs_deg[0]=(46.815*T)+(0.0006*(T*T))-(0.00181*(T*T*T));
	arcs_deg[1]=arcs_deg[0]/3600;
	epsilon_m=23.43929167-arcs_deg[1];
    epsilon_m=epsilon_m+0.02755;
	b=(sin(lambda_m)*cos(deg_torad(epsilon_m)))-(tan(beta_m)*sin(deg_torad(epsilon_m)));
	alpha_m=atan2(b,cos(lambda_m));
	if( alpha_m<0)
	{
		alpha_m=alpha_m+(2*TMath::Pi());
	}
	delta_m=asin((sin(beta_m)*cos(deg_torad(epsilon_m)))+(cos(beta_m)*sin(deg_torad(epsilon_m))*sin(lambda_m)));
	n_m = 0.9856076686 /(pow(1,(3/2))); 
	Me= Mo + (n_m * JD_m);
	Me=sp_algo(Me);
	Oi_m[0]=Me+102.937+(15*(time+(5))); 
	Oi_m[1]=sp_algo(Oi_m[0]);
	theta_m[0]=Oi_m[1]+(lon);
	theta_m[1]=sp_algo(theta_m[0]); 
	H_m=theta_m[1]-rad_todeg(alpha_m);
	D=LmR2-Ls;
    D=sp_algo(D);
    F=(0.5)*(1-cos(deg_torad(D)));
	hei_m=asin(sin(deg_torad(lat))*sin(delta_m)+cos(deg_torad(lat))*cos(delta_m)*cos(deg_torad(H_m)));
	aziS_m=atan2(sin(deg_torad(H_m)),cos(deg_torad(H_m))*sin(deg_torad(lat))-tan(delta_m)*cos(deg_torad(lat)));
	if(aziS_m<=0)
	{
		aziN_m=aziS_m+deg_torad(180);
	}
	else
	{
		aziN_m=aziS_m+TMath::Pi();
	}

/////////////////////////////////////////////////////////////////////////////////
    cout<<"Coordinates:"<<endl;
	cout<<"Longitude: "<<lon<<" degrees"<<endl;
	cout<<"Latitude: "<<lat<<" degrees"<<endl;
	cout<<"Date:"<<a<<"/"<<m<<"/"<<d<<endl;
	cout<<"Time: "<<time<<" horas"<<endl;
	cout<<"Astronomic coordinates:"<<endl;
	cout<<"---- Azimuth ---- Height ---- Phase"<<endl;
	cout<<"Sun:"<<rad_todeg(aziN)<<"|"<<rad_todeg(hei)<<"|"<<"-"<<endl;
	cout<<"Moon:"<<rad_todeg(aziN_m)<<"|"<<rad_todeg(hei_m)<<"|"<<F<<endl;
	
}


