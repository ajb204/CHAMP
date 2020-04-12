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
  double adduction = 0.0;// Adduction
  double Zfudge    = 0.0;  // constant term added to avgZ
  double ResFudge  = 1E-12;   // linear term for m/z peak width
  int minZ         = 1;            // Lowest possible charge state.
  int maxZ         = 200;          // Maximum possible charge state...provide a comfortable margin.


  //Parameters for fitting (this massively effects calculation time, so keep this numbers small)
  int limHSP= 60;     //highest HSP subunit to input
  int limSub= 5;     //largest substract in complex
  int testmax=5000;   //maximum number of trial complexes above threshold (error message comes out if this is too small)


  //Details on the system of interest
  double shspMass = 17985;    //mass of sHSP

  struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};

  //chap=0, and the distribution is a 2D gaussian specified by 4 params
  //chap=1 and the distribution is a 2D gaussian with a skew in the HSP dirction, 5 params
  //chap=2 twisted, slanted Gaussian, 6 params
  //chap=3 Gaussian skewed in both Client and HSP dimensions
  //chap=4 2d skewed Gaussian but with seperate skewed gaussian for Client=0 trace
  //chap=5 Seperate gaussian for each client number trace
  //chap=6 Seperate skewed gaussian for each client number trace


  //chap>10 Read in last fitted distribution and then run completely free fit 

  int chap=2;

  
  // 0 - 1:1 Luc
  // 1 - 1:1 CS
  // 2 - 1:1 MDH
  // 3 - 1:0.1 MDH
  // 4 - 1:0.1 MDH
  // 5 - 1:1 Luc
  // 6 - 1:1 CCS
  // 7 - 1:1 CCS


  int FIT=3;

  //runchap param syntax - distribution type (0-6 +10 as needed) then:
  //fit dist params(1/0),adduction(1/0),Zfudge(1/0),Mres(1/0),MRes_Fudge(1/0)
  //initial guess parameters, struct containing everything useful and client mass

  string raw[10];
  raw[0] = "raw/Q2_FS_020908_6_Luc1_1.fix";
  raw[1] = "raw/Q2_FS_060409_10_CS1_1.fix";
  raw[2] = "raw/Q2_FS_260808_8_MDH1_1.fix3";
  raw[3] = "raw/Q2_FS_270808_5_MDH1_0.1.fix3";
  raw[4] = "raw/Q2_FS_020908_16.fix3";
  raw[5] = 





  

  if(FIT==0)
    {//fit Luc1:1 to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Q2_FS_020908_6_Luc1_1.fix";
      double subMass = 60803; //mass of subunit (Luciferase)
      //double x_init[10]= { 25.0, 5.0, 1.0,0.5,-1.0,0.0001 };  
      double x_init[10]= { 25.0,5.0,1.0,0.5,0.01,   18.0,1.0,1.0,0.01 };  
      runchap(infile,"Luc1_a",0,1,0,0,0,0,x_init,massy,subMass);  //run fit with 2D gauss
      //runchap(infile,"Luc1_b",0,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      //runchap(infile,"Luc1_c",5,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params
      runchap(infile,"Luc1_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"Luc1_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"Luc1_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params
      
      //    runchap(infile,"Luc1_b",4,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate skewed gauss for client=0
      //      runchap(infile,"goLuc1",14,1,1,1,1,1,x_init,massy,subMass);   //run fit with all params free
    }
  
  //    FIT++;
  if(FIT==1)
    {//fit CS1:1 to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Q2_FS_060409_10_CS1_1.fix";
      double subMass = 49072;  //mass of subunit (CS)
      //double x_init[5]= { 25.0, 5.0, 2.0,0.5,10.0 };  
      double x_init[10]= { 25.0,5.0,2.0,0.5,1.0,   18.0,1.0,1.0,0.01 };  
      //      runchap(infile,"CSS1_a",0,1,0,0,0,0,x_init,massy,subMass);
      runchap(infile,"CSS1_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"CSS1_b",1,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"CSS1_b",4,1,1,1,1,1,x_init,massy,subMass);
      runchap(infile,"CSS1_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"CSS1_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"CSS1_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params

      //runchap(infile,"goCSS1",11,1,1,1,1,1,x_init,massy,subMass);


    }
  // FIT++;
  if(FIT==2)
    {//fit MDH:1 to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Q2_FS_260808_8_MDH1_1.fix3";
      double subMass = 33100; //mass of subunit (MDH)
      double x_init[10]= { 25.0,5.0,2.0,0.5,0.01,   18.0,1.0,1.0,0.01};  
      massy.Zfudge= -2.0;
      runchap(infile,"MDH1_a",0,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"MDH1_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"MDH1_c",1,1,1,1,1,1,x_init,massy,subMass);

      runchap(infile,"MDH1_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"MDH1_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"MDH1_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params


      //runchap(infile,"MDH1_b",4,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goMDH1",14,1,1,1,1,1,x_init,massy,subMass);
    }
  //FIT++;
  if(FIT==3)
    {//fit MDH:1 to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Q2_FS_270808_5_MDH1_0.1.fix3";
      double subMass = 33100; //mass of subunit (MDH)
      double x_init[10]= { 25.0,5.0,2.0,0.5,0.01,   18.0,1.0,1.0,0.01 };  
      //double x_init[5]= { 25.0, 5.0, 2.0,0.5,10.0 };  
      runchap(infile,"MDH2_a",0,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"MDH2_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"MDH2_c",1,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"MDH2_b",4,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goMDH2",14,1,1,1,1,1,x_init,massy,subMass);


      //runchap(infile,"MDH2_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"MDH2_b",0,1,1,1,1,1,x_init,massy,subMass);  // run fit with twisted skewed gauss + all params
      runchap(infile,"MDH2_e",7,1,1,1,1,1,x_init,massy,subMass);  // run fit with twisted skewed gauss + all params
      //runchap(infile,"MDH2_e",2,1,1,1,1,1,x_init,massy,subMass);  // run fit with twisted skewed gauss + all params
 
     //runchap(infile,"MDH2_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params


    }

  //FIT=4;
  //  FIT++;
  if(FIT==4)
    {//fit Hsp:Luc to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Q2_FS_020908_16.fix3";
      double subMass = 60803; //mass of subunit (Luc)
      massy.Zfudge=-2.0;
      //double x_init[5]= { 25.0, 5.0,1.0,1.0,0.5 };  
      double x_init[10]= { 25.0,5.0,1.0,1.0,0.01,   18.0,1.0,1.0,0.01 };  
      runchap(infile,"Luc2_a",0,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"Luc2_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"Luc2_b",1,1,1,1,1,1,x_init,massy,subMass);
      runchap(infile,"Luc2_b",1,1,1,1,1,1,x_init,massy,subMass);
      runchap(infile,"Luc2_e",2,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc2",11,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc2",1,1,1,1,1,1,x_init,massy,subMass);

      //runchap(infile,"Luc2_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      //runchap(infile,"Luc2_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      //runchap(infile,"Luc2_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params


    }


  //FIT=6;
  //FIT++;
  if(FIT==5)
    {//fit Hsp:Luc to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Hsp_Luc_1_05_Q2_FS_020908_24_sm12lin1-1.fix";
      double subMass = 60803; //mass of subunit (Luc)
      massy.Zfudge=0.0;
      //double x_init[5]= { 25.0, 5.0,1.0,1.0,0.5 };  
      double x_init[10]= { 25.0,5.0,1.0,1.0,0.01,   18.0,1.0,1.0,0.01 };  
      runchap(infile,"Luc3_a",0,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"Luc3_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"Luc3_c",1,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"Luc3_b",4,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc2",11,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc3",11,1,1,1,1,1,x_init,massy,subMass);


      runchap(infile,"Luc3_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"Luc3_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"Luc3_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params

    }

  //FIT++;
  if(FIT==6)
    {//fit Hsp:Luc to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Hsp_CS_01_better_spectrum.fix";
      double subMass = 49072;  //mass of subunit (CS)
      massy.Zfudge=0.0;
      //double x_init[5]= { 25.0, 5.0,1.0,1.0,0.5 };  
      //double x_init[10]= { 25.0,5.0,1.0,1.0,0.01,   18.0,1.0,1.0,0.01 };  
      double x_init[10]= { 25.0,5.0,2.0,0.5,0.01,   18.0,1.0,1.0,0.01 };  
      runchap(infile,"CSS2_a",0,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"CSS2_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"CSS2_c",1,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"CSS2_b",4,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc2",11,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goCSS2",14,1,1,1,1,1,x_init,massy,subMass);

      //runchap(infile,"CSS2_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"CSS2_b",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"CSS2_e",2,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      //runchap(infile,"CSS2_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params


    }
  //FIT++;
  if(FIT==7)
    {//fit Hsp:Luc to skewed gaussian
      
      struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass,Zfudge,ResFudge};
      string infile="raw/Hsp_CS_01_Used_for_fig.fix";
      double subMass = 49072;  //mass of subunit (CS)
      massy.Zfudge=0.0;
      //double x_init[5]= { 25.0, 5.0,1.0,1.0,0.5 };  
      //      double x_init[10]= { 25.0,5.0,1.0,1.0,0.01,   18.0,1.0,1.0,0.01 };  
      double x_init[10]= { 25.0,5.0,1.0,0.5,0.01,   18.0,1.0,1.0,0.01 };  
      //runchap(infile,"CSS3_a",0,1,0,0,0,0,x_init,massy,subMass);
      runchap(infile,"CSS3_a",1,1,0,0,0,0,x_init,massy,subMass);
      //runchap(infile,"CSS3_c",1,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"CSS3_b",4,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goLuc2",11,1,1,1,1,1,x_init,massy,subMass);
      //runchap(infile,"goCSS3",14,1,1,1,1,1,x_init,massy,subMass);

      runchap(infile,"CSS3_d",1,1,0,0,0,0,x_init,massy,subMass);  // run fit with skewed gauss
      runchap(infile,"CSS3_e",1,1,1,1,1,1,x_init,massy,subMass);  // run fit with skewed gauss + all params
      runchap(infile,"CSS3_f",6,1,1,1,1,1,x_init,massy,subMass);  // run fit with seperate distos + all params


    }





  /*

  // calc distribution and chi^2 for one condition
   {
     string infile="raw/Q2_FS_020908_16.fix3";
     int lines=countlines_string(infile);
     double data[(lines+1)*3];init_array(data,3*(lines+1));    
     readfile(data,lines,infile);            //load in input data
     double subMass = 60803; //mass of subunit (Luciferase)
     double Hx0 = 27.6363;
     double Hsig = 10.666;
     double Sx0   = 1.75336;  
     double Ssig  = 2;     //Sub width
     double alpha= 4.7372;
     double skew = 0.02;
     double par[10];par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;par[4]=alpha;par[5]=skew;       //parameters for chap=0    
     double chi2=run_calc(1,2,".go",3,massy.limHSP,massy.limSub,par,massy.adduction,massy.Zwidth,massy.shspMass,subMass,massy.testmax,massy.minZ,massy.maxZ,massy.thresh,data,lines,massy.MRes,massy.Zfudge,massy.ResFudge);
     cout << "chi2/dof " << chi2 << endl;
   }
   
  */

  
  return 0;
}
