//#include "stdafx.h"
//#include "TRandom3.h"
//#include "TGraph2D.h"
#ifdef __CINT__
int rod_trap()  //put here the name of the file
{gSystem->Load("bemsolver.dll");return main(-1,NULL);}
class D3world;
class D3electrode;
#else
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>

#include "ExpLog.h"
#include "bem.h"
#include "bem2d.h"
#include <string.h>
#include <time.h>
#endif

D3world *wr;
#define ELCNT 5

D3electrode *eldca[ELCNT],*eldcb[ELCNT],*elrfa,*elrfb;


int main(int argc, char* argv[]){
#ifndef __CINT__
	char *str=argv[0];TApplication theApp("App", &argc, argv);
#endif
	
logfile.SetLogFileName();
logfile.Append(false); //funzt!!!!!!!!!
logfile.Overwrite(true); //nicht auf true setzen sondern per hand loeschen
logfile.Logging(false);

TString savedworldfname;
logfile.AbsoluteFileName(TString("./")+logfile.GetLogFileName()+".cache",savedworldfname);
gSystem->mkdir(savedworldfname);
logfile.AbsoluteFileName(TString("./")+logfile.GetLogFileName()+".cache/savedworld.data",savedworldfname);
wr=new D3world(savedworldfname,0.00001,32,6,6);



const int NUM_ELECTRODES=2; //specify number of non-ground electrodes

D3ImportedElectrodes *impel=new D3ImportedElectrodes();
TString importfilename;
logfile.AbsoluteFileName("./rod-trap-drawing_real.dxf",importfilename);
if(!impel->Import(importfilename)) return 0;

ofstream datafile;
ofstream voltagefile;
datafile.open ("./rodtrap_field.txt"); //select file name for data output
voltagefile.open("rod_voltages.txt");

stringstream stm; //string object for int to string conversion.



//inserts all the electrodes into the world, electrodes have to be on layers 0,1,2...
for(int i = 0; i <  NUM_ELECTRODES; i++){
	stm << i;
	string ELname = stm.str();
	wr->insert(& (impel->FindElectrode(ELname.c_str()))); //FindElectrode accepts only c strings.
	stm.str(""); //clears the stm object
}

wr->insert(& (impel->FindElectrode("GROUND")));


wr->refine(150); //sets the resolution
wr->correctNorm(0,0,0);
cout << "Started Solving\n";
wr->solve(); //solves the problem
cout << "Done Solving\n";

//precaches calc_slow(x,y,z) for the specified region of interest
//seems no speedup comes from this!



double zeit = clock();


int Voltages[NUM_ELECTRODES]; //defines voltage array and initializes it to 0
for(int i = 0; i < NUM_ELECTRODES; i++) Voltages[i] = 0;

//specify 1V on one electode at a time and outputs data to a text file for each configuration
//for(int curEl = 0; curEl <= 1; curEl = curEl + 1){  //in this version only first solve for DC electrode, then for all 

for(int curEl = 0; curEl < NUM_ELECTRODES; curEl++){ 
	cout << "\nCurrentElectrode is " << curEl <<"\n";
	for(int i = 0; i < NUM_ELECTRODES; i++) Voltages[i] = 0; //clears the voltage array

	Voltages[curEl] = 1;		
	for(int i = 0; i< NUM_ELECTRODES; i++){
		stm << i;
		string ELname = stm.str();
		impel->FindElectrode(ELname.c_str()).SetVoltage(Voltages[i]);
		stm.str("");
	}

	impel->FindElectrode("GROUND").SetVoltage(0);

	const int num=31; //how many data points per axis,makes this odd to be centered at 0
	double xstart=-50; //defines boundaries
	double xstop=50; 
	//center at -30

	double ystart=-50; 
	double ystop =50;
	//center at 3865

	double zstart =-50; 
	double zstop = 50;
	//center at 208

	double phi[pow(num,3)];
	double phip[pow(num,3)];
	double Ex[pow(num,3)];
	double Ey[pow(num,3)];
	double Ez[pow(num,3)];	
	double coord[3*pow(num,3)]; //coordinates in form [x0 y0 z0 x1 y1 z1...]
	//fills in the coordinate array in the right format
	cout << "FILLING OUT COORDINATE ARRAY\n";
	int xnum=0;
	int ynum=1;
	int znum=2;

	for(int i = 0; i< num; i++){
		for(int j = 0; j< num; j++){
			for(int k = 0; k< num; k++){
			coord[xnum] = xstart+(-xstart+xstop)/double(num-1)*i;
			coord[ynum] = ystart+(-ystart+ystop)/double(num-1)*j;
			coord[znum] = zstart+(-zstart+zstop)/double(num-1)*k;
			xnum+=3;
			ynum+=3;
			znum+=3;
			}
		}
	}
	cout << "DONE FILLING OUT COORDINATE ARRAY\n";


	//performs the calculation and fills in the phi and E arrays
	for(int n = 0; n < pow(num,3);n++){
		phi[n] = wr->calc_slow(coord[3*n],coord[3*n+1],coord[3*n+2]);
		if (n % 300 == 0){
			cout <<"processing point " << n << "out of "<< pow(num,3)<<"\n";
		}
	}


	for(int i = 0; i<NUM_ELECTRODES;i++){
		voltagefile << Voltages[i]<<",";	
	}
	voltagefile <<"\n";
	for(int i =0; i< pow(num,3); i++){
	//datafile << coord[3*i] <<","<< coord[3*i+1] <<","<< coord[3*i+2] <<","<< phi[i] << "," <<Ex[i]<< "," <<Ey[i]<< "," <<Ez[i]<< "\n";
		datafile << coord[3*i] <<","<< coord[3*i+1] <<","<< coord[3*i+2] <<","<< phi[i] << "\n";

	}
}



//		//performs the calculation and fills in the phi and E arrays
//		wr.calc_slow(pow(num,3), coord,phi, Ex, Ey, Ez);
//
//		//print results to a file in the format 
//		//Voltages 1 0 ... 0
//		//x,y,z,phi,Ex,Ey,Ez:
//		//where x is axial
//		//      y is radial
//		//      z is height
//		for(int i = 0; i<NUM_ELECTRODES;i++){
//			voltagefile << Voltages[i]<<",";	
//		}
//		voltagefile <<"\n";
//		for(int i =0; i< pow(num,3); i++){
//			datafile << coord[3*i] <<","<< coord[3*i+1] <<","<< coord[3*i+2] <<","<< phi[i] <<","<<Ex[i]<<","<<Ey[i]<<","<<Ez[i]<<"\n";
//		}
//
//	}
//	




cout << "FINISHED. Took in sec: " << (clock() - zeit)/1000 << "\n";
datafile.close();
voltagefile.close();


logfile.Finish();

#ifndef __CINT__
   while(false){
		gSystem->ProcessEvents();
		gSystem->Sleep(10);
	}
	theApp.Terminate();//root memory error workaround call before destructor of ExpJob can be called!
#endif
}

