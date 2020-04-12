#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_erf.h>


using namespace std;
using std::string;


/*###############################################################*/
//#
//# Distyfit.py - for your Hsp:Client fitting needs
//#
//# This version is rigged for fitting sHSPs with no clients
//#
//# Note on modes:
//# mode=0 - 1D gaussian. gridsearch performed.
//# mode=7 - single oligomeric species #xinit0 is species, no gridsearch
//#
//# OLDER VERSION:
//#mode=0, and the distribution is a 2D gaussian specified by 4 params
//#mode=1 and the distribution is a 2D gaussian with a skew in the HSP dirction, 5 params
//#mode=2 twisted, slanted Gaussian, 6 params
//#mode=3 Gaussian skewed in both Client and HSP dimensions
//#mode=4 2d skewed Gaussian but with seperate skewed gaussian for Client=0 trace
//#mode=5 Seperate gaussian for each client number trace
//#mode=6 Seperate skewed gaussian for each client number trace
//#mode=7 AlphaB model (rat and mult)  //HERE THIS HAS BEEN RIGGED TO BE SINGLE COMPLEX ONLY
//#mode>10 Read in last fitted distribution and then run completely free fit 

  
#include "disty_aux.c"
#include "disty_fit.c"


int main(int argc, char *argv[])
{
  if (argc !=2) //Make sure we will get the arugments we need
    {
      cout << "Usage: " << argv[0] << " inputfile" << endl;
      exit(100);
    }

  char *inputfile= argv[1];
  int fileno=countlines(inputfile)-2; //number of datasets in dataset file

  //declare relevant arrays
  string raw[fileno],ident[fileno];
  double shspMass[fileno],subM[fileno],minMZ[fileno],maxMZ[fileno],xinit0[fileno],xinit1[fileno];
  int mode[fileno],smooth[fileno];
  double Zwidth[fileno],MRes[fileno],Zfudge[fileno],max0[fileno],min0[fileno];

  fileno=ReadInput(inputfile,fileno,raw,ident,smooth,shspMass,minMZ,maxMZ,mode,xinit0,xinit1,Zwidth,MRes,Zfudge,max0,min0);//read in dataset summary

  //Formerly read-in parser for values.txt//char *valuefile= argv[2];//double adduction=ReadValueDouble(valuefile,"adduction");

  //these guys are currently not under individual user control
  double adduction = 1E-3;
  double thresh    = 0.00001;   // Abundance threshold for inclusion of a charge state.
  double ResFudge  = 1E-12;   // linear term for m/z peak width
  int minZ         = 1;       // Lowest possible charge state.
  int maxZ         = 200;     // Maximum possible charge state...provide a comfortable margin.
  //Parameters for fitting (this massively effects calculation time, so keep this numbers small)
  int limHSP= 30;    //highest HSP subunit to input
  int limSub= 35;     //largest substract in complex
  int testmax=50000;  //maximum number of trial complexes above threshold (error message comes out if this is too small)
  //Details on the system of interest
  double tw = 0.1;           //relative proportion of residual (12,0) mer (if required)

  initfile("figs/spectrafit.gp");   //initialise gnuplot file
  initfile("figs/make_summary.com");//initialise summary file

  

  FILE*fp;
  fp=fopen("output.txt","w");
  fclose(fp);
  for(int i=0;i<fileno;i++)
    {
      make_summary_init("figs/make_summary.com");  //write first line of summary file
      //for(int j=0;j<2;j++)
		int j=0;
	{
	  if(j==0)
	    mode[i]=0; //run in gaussian mode
	  if(j==2) //Poisson
	      mode[i]=2;
	  if(j==1) //run in single species mode
	    mode[i]=7;


	  int trial=0; //initialise run number [ultimately bring this into the struct]
	  //for gridsearch
	  double max1=0.0; //max for gaussian width grid
	  double min1=0.0; //min for gaussian width grid
	  char *GridType="lin"; //either linear (lin) or logarithmic (log) scaling of second dimension

	  if(mode[i]==0){//gaussian mode
	    max1=5.0;
	    min1=0.1;
	    trial=0;
	    GridType="log";}
	  if(mode[i]==2){//alphaB Poisson mode
	    max1=20.0;
	    min1=0.0;
	    trial=1;
	    GridType="lin";}
	  if(mode[i]==7){//single species mode
	    max1=20.0;
	    min1=0.0;
	    trial=1;
	    GridType="lin";}



	    
	  double x_init[limHSP*limSub];//Set initial conditions
	  int lines=analfile(raw[i],minMZ[i],maxMZ[i]);  //count number of relevant datalines
	  double data[lines*3/smooth[i]];init_array(data,3*lines/smooth[i]); //declare data array
	  lines=readfile(data,lines,raw[i],smooth[i],minMZ[i],maxMZ[i]); //load in input data and smooth
	  normdata(data,lines);                     //normalise input data
	  prin_spec(raw[i],".out",4,data,lines,1);  //print file with adjusted input spectrum
	  
	  x_init[0]=xinit0[i]; //first distribution parameter
	  x_init[1]=xinit1[i]; //second distribution parameter

	  //subM[i] = 0.0 ;  //client is irrelevant here
	  
	  subM[i]=22782;
	  shspMass[i]=20159;
	  mode[i]=10;

	  //note for Andy- changed distributions, gridsearch function and which grads are set to zero to make 1D program
	  struct mass massy={thresh,Zwidth[i],minMZ[i],maxMZ[i],MRes[i],adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass[i],subM[i],Zfudge[i],ResFudge,tw};


	  struct control params={mode[i],0,0,1,1,1,0};//chap, distparam, tw_flg,adduct_flg,Zfudge_flg,Mres_flg,ResFudge_flg
	  //struct control params={mode[i],1,0,1,1,1,0};//run type,(ignore) then: 12mer flag,adduction,Zfudge,Mres,ResFudge
	  //so x,x,0,1,1,1,0  will fit adduction,zfudge and Mres, and not ResFudge nor free 12mer
	  
	    
	  
	  if(mode[i]==0 || mode[i]==2)//do the grid search if doing gaussian fit
	    gridsearch1D(i,mode[i],trial,lines,data,x_init,massy,params,raw,ident,max0[i],min0[i],max1,min1,GridType); //run gridsearch and update x_init
	  if(mode[i]==7)
	   gridsearchSingle(i,mode[i],trial,lines,data,x_init,massy,params,raw,ident,max0[i],min0[i]); //run gridsearch and update x_init
	  
	  
	  //should be good to get here: need to read in input distribution to initialise distribution

	  double chi2=runchap("sim",trial,lines,data,ident[i],params,x_init,massy,1);trial++;  //run minimisation (j specifies skew or not)

	  chi2=runchap("fitty",trial,lines,data,ident[i],params,x_init,massy,1);trial++;  //run minimisation (j specifies skew or not)


	  //x_init[0]=xinit0[i]; //first distribution parameter
	  //x_init[1]=xinit1[i]; //second distribution parameter	  
	  //x_init[2]=1.0;
	  //params.chap=1;
	  //chi2=runchap(trial,lines,data,ident[i],params,x_init,massy,1);trial++;  //run minimisation (j specifies skew or not)
	  //params.chap=11;
	  //chi2=runchap(trial,lines,data,ident[i],params,x_init,massy,1);trial++;  //run minimisation (j specifies skew or not)


	  //cout << "Running with second species... " << endl;
	  //cout << trial << endl;

	  
	  //adding a second species...  //need to grid search on this.
	 	  
	  fp=fopen("output.txt","a");
	  fprintf(fp,"%s\t%s\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",raw[i].c_str(),ident[i].c_str(),mode[i],chi2,massy.MRes,massy.Zfudge,massy.adduction,shspMass[i],x_init[0],x_init[1] );
	  fclose(fp);
	  
	}
      //FINISH RUN FOR GIVEN FILE
      make_summary_end("figs/make_summary.com",ident[i]);
    }//FINISHED
  return 0;
}
