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




void test()
{

	double JD=0; //Julian Date Number
	double r=0; //The distance of the Earth from the Sun
	double v=0; //True Anomaly of the Earth
	double hX=0; //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	double hY=0; //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	double hZ=0; //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	double gX=0; //Rectangular Geocentric Ecliptical Coordinates. X Axis. 
	double gY=0; //Rectangular Geocentric Ecliptical Coordinates. Y Axis.
	double gZ=0; //Rectangular Geocentric Ecliptical Coordinates. Z Axis.
	double delta=0; //The distance of the Sun from the Earth
	double lambda=0; //The ecliptic longitude
	double beta=0; //The ecliptic latitude
	double alpha=0; //The right ascension
	double delta_low=0; //The declination
	double hei=0; //The height
    double aziS=0; //The Azimuth from the South
	double aziN=0; //The Azimuth from the North
	
   
	int prof[3];
	double c_t[4];
	//cout<<"Enter the date in Gregorian mode"<<endl;
	//cout<<"Year:"<<endl;
	//cin>>prof[0];
	prof[0]=2019;
	//cout<<"Month:"<<endl;
	prof[1]=03;
	//cout<<"Day:"<<endl;
	prof[2]=28;
	//cout<<"Enter coordinates from the Greenwich meridian and Equator frame"<<endl;
	//cout<<"Longitude:"<<endl;
	//cin>>c_t[1];
	c_t[1]=-73.370721;
	//cout<<"Latitude:"<<endl;
	//cin>>c_t[2];
	c_t[2]=5.865964;
    //cout<<"Enter the actual time in hours"<<endl;
	cout<<"Time(hours):"<<endl;
    cin>>c_t[0];
	//c_t[0]=18; 
	cout<<"Time(minutes):"<<endl;
	cin>>c_t[3];
	//c_t[3]=29;
	cout<<"Time format: "<<c_t[0]<<":"<<c_t[3]<<endl;
	c_t[0]=c_t[0]+(c_t[3]/60);   // Convert sexagesimal hour record to decimal record
	cout<<"Decimal format: "<<c_t[0]<<endl;
	
	
	
	JD = date_conversion(prof[0],prof[1],prof[2]);
	
	r=distance_tosun(true_anomaly(mean_anomaly(JD)));
	
	v=true_anomaly(mean_anomaly(JD));
	
	helio_coor(r,v,hX,hY,hZ);
	
	geo_coor(hX,hY,hZ,gX,gY,gZ);
	
	GElonlat(gX,gY,gZ,delta,lambda,beta);
	
	updown(lambda,beta,alpha,delta_low);
	
	H_A(hour_angle(side_time(c_t[0],c_t[1],JD),rad_todeg(alpha)),c_t[2],delta_low,hei,aziS,aziN);
	
    /*
	//Showing results
	cout<<"The date in Julian Day Number: "<<JD<<"days"<<endl;
	cout<<"Mean Anomaly of the Earth: "<<mean_anomaly(JD)<<" degrees"<<endl;
	cout<<"True Anomaly of the Earth: "<<true_anomaly(mean_anomaly(JD))<<" degrees"<<endl;
	cout<<"The distance between Earth and the Sun: "<<distance_tosun(true_anomaly(mean_anomaly(JD)))<<" AU"<<endl;
	cout<<"The Rectangular Heliocentric Ecliptical Coordinate in the X Axis: "<<hX<<" AU"<<endl;
	cout<<"The Rectangular Heliocentric Ecliptical Coordinate in the Y Axis: "<<hY<<" AU"<<endl;
	cout<<"The Rectangular Heliocentric Ecliptical Coordinate in the Z Axis: "<<hZ<<" AU"<<endl;
	cout<<"The Rectangular Geocentric Ecliptical Coordinate in the X Axis: "<<gX<<" AU"<<endl;
	cout<<"The Rectangular Geocentric Ecliptical Coordinate in the Y Axis: "<<gY<<" AU"<<endl;
	cout<<"The Rectangular Geocentric Ecliptical Coordinate in the Z Axis: "<<gZ<<" AU"<<endl;
	cout<<"The distance of the Sun from the Earth Delta: "<<delta<<" AU"<<endl;
	cout<<"The ecliptic latitude Beta: "<<rad_todeg(beta)<<" degrees"<<endl;
	cout<<"The ecliptic longitude Lambda: "<<rad_todeg(lambda)<<" degrees"<<endl;
	cout<<"The right ascension Alpha: "<<rad_todeg(alpha)<<" degrees"<<endl;
	cout<<"The The declination Delta(lowercase): "<<rad_todeg(delta_low)<<" degrees"<<endl;
	cout<<"The The Sidereal time Theta: "<<side_time(c_t[0],c_t[1],JD)<<" degrees"<<endl;
	cout<<"The hour angle: "<<hour_angle(side_time(c_t[0],c_t[1],JD),rad_todeg(alpha))<<" degrees"<<endl;
	*/
	cout<<"The height above the horizon: "<<rad_todeg(hei)<<" degrees"<<endl;
	cout<<"The azimuth from the North: "<<rad_todeg(aziN)<<" degrees"<<endl;
	cout<<"The date in Julian Day Number: "<<JD<<"days"<<endl;
	//cout<<"The azimuth from the South: "<<rad_todeg(aziS)<<" degrees"<<endl;
	
	/*
	
	
	cout<<"Made out by..."<<endl;
	cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"**   "<<"     "<<"   **"<<"*****"<<"**   "<<"     "<<"*************"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"***  "<<"     "<<"  ***"<<"*****"<<"***  "<<"     "<<"*************"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"**** "<<"     "<<" ****"<<"*****"<<"**** "<<"     "<<"*****        "<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"*****"<<"*****"<<"     "<<"******       "<<endl;
    cout<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"********     "<<endl;
    cout<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*************"<<endl;
    cout<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"  ***"<<"*****"<<"     "<<"*************"<<endl;
    cout<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"  ***"<<"*****"<<"     "<<"        *****"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"*****"<<"     "<<"*****"<<"*****"<<"***  "<<"     "<<"        *****"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"**** "<<"     "<<" ****"<<"*****"<<"**** "<<"     "<<"     ********"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"***  "<<"     "<<"  ***"<<"*****"<<"  ***"<<"     "<<"*************"<<endl;
    cout<<"*****"<<"     "<<"*****"<<"     "<<"*****"<<"*****"<<"**   "<<"     "<<"   **"<<"*****"<<"  ***"<<"     "<<"*************"<<endl;

    */
	
	//gSystem->Exec("mkdir new_folder"); //This command allows you to make a new directory. Like mkdir.
	
	
	
	
	///// Example of how to convert numbers to strings
	
	
	/*
	int Number = 123;       // number to be converted to a string

    string Result;          // string which will contain the result

    ostringstream convert;   // stream used for the conversion

    convert << Number;      // insert the textual representation of 'Number' in the characters in the stream

    Result = convert.str(); // set 'Result' to the contents of the stream

    // 'Result' now is equal to "123" 
	
	*/
	
	
	
	
	
    
}
	

double date_conversion(int a, int m, int d)  // This function transforms Gregorian Date to Julian Date     Only Date,not time
{
	double JD=( 367*a- 7 * ( a + (m+9)/12 ) / 4 - 3 * ( ( a + (m-9)/7 ) / 100 + 1 ) / 4 + 275*m/9 + d - 730515)-1.5;   // At noon UTC on 1 January 2000, at Julian Day 2451545. Check http://www.onlineconversion.com/julian_date.htm 
	return JD;
}

double mean_anomaly(double JD)  //This function calculates the Mean Anomaly of the Earth's orbit. All data in degrees
{
	double M; // Mean Anomaly variable
	double n;// Angle (n) of the Earth traverses on average per day. Measured in degrees
	n = 0.9856076686 /(pow(1,(3/2)));  
	double Mo=357.529; // Mean Anomaly at specific date. At noon UTC on 1 January 2000, at Julian Day 2451545 Mean Anomaly_0 of the Earth
	M= Mo + (n * JD);
	
	return algo(M);
}

double true_anomaly(double M) //This function calculates the True Anomaly of the Earth's orbit. All data in degrees
{
	float e=0.01671; //Eccentricity of the Earth's orbit, an ellipse. 
	double v; //True Anomaly of the Earth. This measure is a iteration. Check the calculations.
	
	v= deg_torad(M) +(((2*e)-((1/4)*pow(e,3))+((5/96)*pow(e,5))+((107/4608)*pow(e,7)))*sin(deg_torad(M))) 
	+ ((((5/4)*pow(e,2))-((11/24)*pow(e,4))+((17/192)*pow(e,6)))*sin(2*deg_torad(M)))
	+((((13/12)*pow(e,3))-((43/64)*pow(e,5))+((95/512)*pow(e,7)))*sin(3*deg_torad(M)))
	+((((103/96)*pow(e,4))-((451/480)*pow(e,6)))*sin(4*deg_torad(M)))
	+((((1097/960)*pow(e,5))-((5957/4680)*pow(e,7)))*sin(5*deg_torad(M)))
	+ (((1223/960)*pow(e,6))*sin(6*deg_torad(M)))
	+(((47273/32256)*pow(e,7))*sin(7*deg_torad(M)));
	
	v=rad_todeg(v);
	
	
	return v;
}
	
double distance_tosun(double v) //This function calculates the distance between the Sun and the Earth. All data in AU (Astronomical Units).
{
	double r; //The distance of the Earth from the Sun
	float e=0.01671 ; //Eccentricity of the Earth's orbit, an ellipse. 
	r=(1*(1-(pow(e,2))))/(1+(e*cos(deg_torad(v)))); // The distance between Earth and the Sun, measured in AU
	
	return r;
}

void  helio_coor(double ratio,double v,double& hX,double& hY,double& hZ)    //The Rectangular Heliocentric Ecliptical Coordinates calculation.
{
	double O=174.873; //The ecliptic longitude Ohm [large omega] of the ascending node of the orbit.
	double w=288.064; // The argument w[omega] of the perihelion
	double i=0.000; //the inclination i of the orbit
	
	hX=ratio*((cos(deg_torad(O))*cos(deg_torad(w+v)))-(sin(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v)))); //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	hY=ratio*((sin(deg_torad(O))*cos(deg_torad(w+v)))+(cos(deg_torad(O))*cos(deg_torad(i))*sin(deg_torad(w+v)))); //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	hZ=ratio*sin(deg_torad(i))*sin(deg_torad(w+v)); //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
}

void geo_coor(double hX,double hY,double hZ,double& gX,double& gY,double& gZ)         //    The Rectangular Geocentric Ecliptical Coordinates calculation
{
	// This fuction asume that XSun,YSun and ZSun are equal to zero due to the reference frame.
	gX=0-hX;
	gY=0-hY;
	gZ=0-hZ;
}
	
void GElonlat(double gX,double gY, double gZ,double& delta,double& lambda,double& beta) //Geocentric Ecliptical Longitude and Latitude calculation	
{
	double delta=sqrt((pow(gX,2))+(pow(gY,2))+(pow(gZ,2))); //The distance large delta of the planet from the Earth
	double lambda=atan2(gY,gX);    //The ecliptic longitude
	double beta=asin(gZ/delta);    //The ecliptic latitude
	
    if( lambda<0)
	{
		lambda=lambda+(2*TMath::Pi());
	}
	
}

void updown(double lambda,double beta,double& alpha,double& delta_low)// Right Ascension and Declination calculation
{
	double epsilon=deg_torad(23.4397); //Obliquity of the ecliptic
	double alpha=atan2(sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon),cos(lambda)); //The right ascension
	double delta_low=asin((sin(beta)*cos(epsilon))+(cos(beta)*sin(epsilon)*sin(lambda))); //The declination
	
	if( alpha<0)
	{
		alpha=alpha+(2*TMath::Pi());
	}

}

double side_time(double time,double lon,double JD) // The Sidereal Time calculation
{
	double Oi[2];        //Sidereal time on the prime meridian
	Oi[0]=mean_anomaly(JD)+102.937+(15*(time+(5))); //5 plus hours for Colombia
	Oi[1]=algo(Oi[0]);
	double theta[2];
	theta[0]=Oi[1]+(lon); //The sideral time
	theta[1]=algo(theta[0]); 
	
	return theta[1];
}

double hour_angle(double theta, double alpha) //Hour Angle calculation
{
	double H; //Hour Angle
	H=theta-alpha;
	return H;
}

void H_A(double H,double lat,double delta_low,double& hei,double& aziS,double& aziN)  //Height and Azimuth calculation
{
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
}

void calendar() //Testing and testing...
{
	
	int year;
	year=2019;
    int march[31];
	int april[30];
	int may[31];
	int hour[24];
	double coor[2]; // Latitude, Longitude
	coor[0]=5.865964;
	coor[1]=-73.370721;
	
	
	
	for(int i=0;i<31;i++) // Loop to fill the array (March and May)
	{
		march[i]=i+1;
		may[i]=i+1;
		//cout<<"Marzo: "<<march[i]<<" and "<<"Mayo: "<<may[i]<<endl;
	}
	
    for(int i=0;i<30;i++) // Loop to fill the array (April)
	{
	    april[i]=i+1;
		//cout<<"April "<<april[i]<<endl;
	}
	
	for(int i=0;i<24;i++) // Loop to fill the array (Hour)
	{
	    hour[i]=i;
		//cout<<"Hora "<<hour[i]<<endl;
	}

	
	
	double azimutal[24];
	double altitud[24];
	
	

	for(int i=0;i<1;i++)
	{
		
		double JD=0; //Julian Date Number
	    double r=0; //The distance of the Earth from the Sun
	    double v=0; //True Anomaly of the Earth
	    double hX=0; //Rectangular Heliocentric Ecliptical Coordinates. X Axis.
	    double hY=0; //Rectangular Heliocentric Ecliptical Coordinates. Y Axis.
	    double hZ=0; //Rectangular Heliocentric Ecliptical Coordinates. Z Axis.
	    double gX=0; //Rectangular Geocentric Ecliptical Coordinates. X Axis. 
	    double gY=0; //Rectangular Geocentric Ecliptical Coordinates. Y Axis.
	    double gZ=0; //Rectangular Geocentric Ecliptical Coordinates. Z Axis.
	    double delta=0; //The distance of the Sun from the Earth
	    double lambda=0; //The ecliptic longitude
	    double beta=0; //The ecliptic latitude
	    double alpha=0; //The right ascension
	    double delta_low=0; //The declination
	    double hei=0; //The height
        double aziS=0; //The Azimuth from the South
	    double aziN=0; //The Azimuth from the North
		
		
		
		
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		JD = date_conversion(year,4,april[16]);
		
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	    r=distance_tosun(true_anomaly(mean_anomaly(JD)));
	
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    v=true_anomaly(mean_anomaly(JD));
		
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	    helio_coor(r,v,hX,hY,hZ); 
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	    geo_coor(hX,hY,hZ,gX,gY,gZ);
		
		cout<<gX<<endl;
		cout<<gY<<endl;
		cout<<gZ<<endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //GElonlat(gX,gY,gZ,delta,lambda,beta);
		
		
	    delta=sqrt((pow(gX,2))+(pow(gY,2))+(pow(gZ,2))); //The distance large delta of the planet from the Earth
        lambda=atan2(gY,gX);    //The ecliptic longitude
	    beta=asin(gZ/delta);   //The ecliptic latitude
		
		if( lambda<0)
	    {
		  //  lambda=lambda+(2*TMath::Pi());
	    }
	    
		cout<<delta<<endl;
		cout<<lambda<<endl;
		cout<<beta<<endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //updown(lambda,beta,alpha,delta_low);
	    
		double epsilon=deg_torad(23.4397); //Obliquity of the ecliptic
	    alpha=atan2(sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon),cos(lambda)); //The right ascension
	    delta_low=asin((sin(beta)*cos(epsilon))+(cos(beta)*sin(epsilon)*sin(lambda))); //The declination
	
	    if( alpha<0)
	    {
		    alpha=alpha+(2*TMath::Pi());
	    }
 
	    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	    H_A(hour_angle(side_time(hour[i],coor[1],JD),rad_todeg(alpha)),coor[0],delta_low,hei,aziS,aziN);
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		azimutal[i]=aziN;
		altitud[i]=hei;
		
	}
	
	ofstream salida;
	salida.open("un dia X.txt");
	salida<<hour[0]<<" "<<azimutal[0]<<" "<<altitud[0]<<endl;
	salida<<hour[1]<<" "<<azimutal[1]<<" "<<altitud[1]<<endl;
	salida<<hour[2]<<" "<<azimutal[2]<<" "<<altitud[2]<<endl;
	salida<<hour[3]<<" "<<azimutal[3]<<" "<<altitud[3]<<endl;
	salida<<hour[4]<<" "<<azimutal[4]<<" "<<altitud[4]<<endl;
	salida<<hour[5]<<" "<<azimutal[5]<<" "<<altitud[5]<<endl;
	salida<<hour[6]<<" "<<azimutal[6]<<" "<<altitud[6]<<endl;
	salida<<hour[7]<<" "<<azimutal[7]<<" "<<altitud[7]<<endl;
	salida<<hour[8]<<" "<<azimutal[8]<<" "<<altitud[8]<<endl;
	salida<<hour[9]<<" "<<azimutal[9]<<" "<<altitud[9]<<endl;
	salida<<hour[10]<<" "<<azimutal[10]<<" "<<altitud[10]<<endl;
	salida<<hour[11]<<" "<<azimutal[11]<<" "<<altitud[11]<<endl;
	salida<<hour[12]<<" "<<azimutal[12]<<" "<<altitud[12]<<endl;
	salida<<hour[13]<<" "<<azimutal[13]<<" "<<altitud[13]<<endl;
	salida<<hour[14]<<" "<<azimutal[14]<<" "<<altitud[14]<<endl;
	salida<<hour[15]<<" "<<azimutal[15]<<" "<<altitud[15]<<endl;
	salida<<hour[16]<<" "<<azimutal[16]<<" "<<altitud[16]<<endl;
	salida<<hour[17]<<" "<<azimutal[17]<<" "<<altitud[17]<<endl;
	salida<<hour[18]<<" "<<azimutal[18]<<" "<<altitud[18]<<endl;
	salida<<hour[19]<<" "<<azimutal[19]<<" "<<altitud[19]<<endl;
	salida<<hour[20]<<" "<<azimutal[20]<<" "<<altitud[20]<<endl;
	salida<<hour[21]<<" "<<azimutal[21]<<" "<<altitud[21]<<endl;
	salida<<hour[22]<<" "<<azimutal[22]<<" "<<altitud[22]<<endl;
	salida<<hour[23]<<" "<<azimutal[23]<<" "<<altitud[23]<<endl;
	
	
	c1= new TCanvas("c1","First testing",1200,700);
    c1->Divide(3,1);
	
	c1->cd(1);
	a_ho= new TGraph();
	for(int i=0;i<24;i++)
	{
		
		a_ho->SetPoint(i,hour[i],azimutal[i]);
		
	}
	a_ho->SetTitle("Azimuth vs Tiempo;Tiempo(h);Azimuth Norte(grados)");
	a_ho->SetLineColor(kRed);
	a_ho->SetMarkerStyle(7);
	a_ho->Draw("ACP");

	
	
	
	c1->cd(2);
	h_ho= new TGraph();
	for(int i=0;i<24;i++)
	{
		
		h_ho->SetPoint(i,hour[i],altitud[i]);
		
	}
	h_ho->SetTitle("Altitud vs Tiempo;Tiempo(h);Altitud desde el horizonte(grados)");
	h_ho->SetLineColor(kRed);
	h_ho->SetMarkerStyle(7);
	h_ho->Draw("ACP");
	
	c1->cd(3);
	a_h= new TGraph();
	for(int i=0;i<24;i++)
	{
		
		a_h->SetPoint(i,altitud[i],azimutal[i]);
		
	}
	a_h->SetTitle("Azimuth vs Altitud;Altitud(grados);Azimuth(grados)");
	a_h->SetLineColor(kRed);
	a_h->SetMarkerStyle(7);
	a_h->Draw("ACP");
	
	string rolo;
	rolo="fuck";
	
   	c1->SaveAs("Dia.png");
	
	
	
	
	
	
	
}