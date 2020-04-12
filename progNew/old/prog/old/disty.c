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

/* Globals  */
double Rgas=8.3145;  // gas constant J K-1 mol-1
double mon=24E-6;    // total subunit concentration
double T0=273.19;    // 0 K in oC

/*###############################################################*/
//#
//# Distyfit.py - for your Hsp:Client fitting needs
//#
//# AJB 040510    Added stuff to do things.
//# AJB 020410    Added distribution function where each slice has free distribution to win bet with Justin
//# AJB 100316.02 Skewed Gaussian (chap=1) and regular distribution (chap=0) added
//# AJB 100316.01 Levenberg-Marquet least-squares minimiser added
//# AJB 100314.01 Ported into C++
//# AJB 100313.01 Automated the mo'fo'.
//# MFB 090526.01.SimMS.2D Modifying to read 2D files
//# MFB 090526.02 m/z increments now depend on m/z and resolution
//# MFB 090526.01 Fixed last version, further improved speed
//# MFB 090520.02 Working on speed, but there is an error in code.
//# MFB 090520.01 Cleaning up code before sharing
//# MFB 090424.01 Beginning work to simulate spectra


#include "disty_aux.c"
#include "disty_fit.c"


int main()
{
  //Mass spec parameters
  double thresh    = 0.001;   // Abundance threshold for inclusion of a charge state.
  double Zwidth    = 1.3;     // Width of the charge state distribution, value from Justin = 1.1352
  double minMZ     = 7500.0;  // Minimum m/z to plot.
  double maxMZ     = 13000;   // Maximum m/z to plot.
  double MRes      = 800;     // Mass spectral resolving power (FWHM/(m/z))
  double adduction = 0.0;     // Adduction
  double Zfudge    = 0.0;     // constant term added to avgZ
  double ResFudge  = 1E-12;   // linear term for m/z peak width
  int minZ         = 1;       // Lowest possible charge state.
  int maxZ         = 200;     // Maximum possible charge state...provide a comfortable margin.

  //Parameters for fitting (this massively effects calculation time, so keep this numbers small)
  int limHSP= 60;    //highest HSP subunit to input
  int limSub= 7;     //largest substract in complex
  int testmax=5000;  //maximum number of trial complexes above threshold (error message comes out if this is too small)

  //Details on the system of interest
  double shspMass = 16167.4;   //mass of sHSP
  double tw = 0.0;           //relative proportion of residual (12,0) mer

  //chap=0, and the distribution is a 2D gaussian specified by 4 params
  //chap=1 and the distribution is a 2D gaussian with a skew in the HSP dirction, 5 params
  //chap=2 twisted, slanted Gaussian, 6 params
  //chap=3 Gaussian skewed in both Client and HSP dimensions
  //chap=4 2d skewed Gaussian but with seperate skewed gaussian for Client=0 trace
  //chap=5 Seperate gaussian for each client number trace
  //chap=6 Seperate skewed gaussian for each client number trace
  //chap>10 Read in last fitted distribution and then run completely free fit 

  //runchap param syntax - distribution type (0-6 +10 as needed) then:
  //fit dist params(1/0),adduction(1/0),Zfudge(1/0),Mres(1/0),MRes_Fudge(1/0)
  //initial guess parameters, struct containing everything useful and client mass

  int fileno=2;   //number of files to be taken in
  string raw[fileno];double subM[fileno];string ident[fileno];
  //raw[0] = "raw/test.txt";       
subM[0] = 60803 ;  ident[0]="Yp1a";
raw[0] = "raw/Luc_1_01_sm1lin1.tst";

  raw[1] = "raw/Luc_1_01_sm1lin1.tst2";      subM[0] = 60803 ;  ident[0]="Luc1";// 0 - 1:0.1 Luc
  //raw[1] = "raw/Luc_1_05_sm1lin1.txt.fix";   subM[1] = 60803 ;  ident[1]="Luc2";// 1 - 1:0.5 Luc
  //raw[2] = "raw/Luc_1_1_sm1lin1.txt.fix";    subM[2] = 60803 ;  ident[2]="Luc3";// 2 - 1:1   Luc
  //raw[3] = "raw/CS_1_01_sm1lin1.txt.fix";    subM[3] = 49072 ;  ident[3]="CSS1";// 3 - 1:0.1 CSS
  //raw[4] = "raw/CS_1_1_sm1lin1.txt.fix";     subM[4] = 49072 ;  ident[4]="CSS2";// 4 - 1:1 CSS
  //raw[5] = "raw/MDH_1_01_sm1lin1.txt.fix";   subM[5] = 33100 ;  ident[5]="MDH1";// 5 - 1:0.1 MDH
  //raw[6] = "raw/MDH_1_1_sm1lin1.txt.fix";    subM[6] = 33100 ;  ident[6]="MDH2";// 6 - 1:1 MDH

  initfile("figs/spectrafit.gp");   //initialise gnuplot file
  initfile("figs/make_summary.com");//initialise summary file
  cout << "Hello! I am here!" << endl;
  printf("shit\n");


  FILE*fp;
  //   for (int i=0;i<fileno;i++)//select which file
  for (int ix=0;ix<1;ix++)//select which file
    {
         int i=0;
      if(ix==0)
      	i=0;
      

      int trial=0;
      int tw_flg;
      make_summary_init("figs/make_summary.com");  //write first line of summary file
      
      if(i==5 || i==6)  //decide whether or not to add free 12mer (luc 1:0.1, MDH 1:0.1 and MDH 1:1)
	tw_flg=1;
      else
	tw_flg=0;
      
      
      for (int j=0;j<1;j++)//j=0 for Gaussian, j=1 for Skewed Gaussian
	{
	  j=1;//skewed Gaussian
	  double x_init[limHSP*limSub];//Set initial conditions
	  
	  {//START THE MEAT
	    
	    struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge,tw};
	    gridsearch(i,j,trial,x_init,massy,raw,ident,subM,tw_flg); //run gridsearch and update x_init


	    // values for Luc1
	    //x_init[0]=23.5;
	    //x_init[2]=0.75;
	    // values for MDH2
	    // x_init[0]=23.5;
	    // x_init[2]=2.25;
	   
	    //x_init[1]=2.0;
	    //x_init[3]=0.5;
	    //x_init[4]=0.01;

 
	    
	    double chi2=runchap(trial,raw[i],ident[i],j,1,tw_flg,1,1,1,0,x_init,massy,subM[i]);trial++;  //run minimisation (j specifies skew or not)
	    chi2=runchap(trial,raw[i],ident[i],j+5 ,1,tw_flg,1,1,1,0,x_init,massy,subM[i]);trial++;  // sliced fitting, taking in previous fits values
	    //run completely free minimisation
	    chi2=runchap(trial,raw[i],ident[i],j+10 ,1,tw_flg,1,0,1,0,x_init,massy,subM[i]);trial++;  // skewed gaussian with client=0 seperate + all params
	    
	  }
	}//FINISH RUN FOR GIVEN FILE
      make_summary_end("figs/make_summary.com",ident[i]);
      
      
    }//FINISH RUN FOR GIVEN FILE
  
  
  
  
  
  return 0;
}
