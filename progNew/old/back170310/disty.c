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
//# AJB 100314.01 Ported into C++
//# AJB 100313.01 Automated the mo'fo'.
//# MFB 090526.01.SimMS.2D Modifying to read 2D files
//# MFB 090526.02 m/z increments now depend on m/z and resolution
//# MFB 090526.01 Fixed last version, further improved speed
//# MFB 090520.02 Working on speed, but there is an error in code.
//# MFB 090520.01 Cleaning up code before sharing

//# MFB 090424.01 Beginning work to simulate spectra


#include "disty_aux.c"



double run_calc(int flg,int chap,char* lab,int tag,int limHSP,int limSub,double *par,double adduction,double Zwidth,double shspMass,double subMass,int testmax,int minZ,int maxZ,double thresh,double* data,int lines,double MRes)
{
  double input[limHSP*limSub];init_array(input,limHSP*limSub);
  double Complexes[limHSP*limSub*5];init_array(Complexes,limHSP*limSub*5);
  double Complex_test[testmax*3];init_array(Complex_test,testmax*3);


  reset_sim(data,lines);   //reset calculated spectrum
  make_input(input,limHSP,limSub,chap,par);                                           //Generate input matrix
  complex_anal(input,Complexes,limHSP,limSub,adduction,Zwidth,shspMass,subMass);      //analyse input species
  int cnt=complex_def(Complexes,Complex_test,limHSP,limSub,testmax,minZ,maxZ,thresh); //find all charge states for input distribution
  double chi2=eval_spec(Complex_test,cnt,testmax,data,lines,MRes);                    //evaluate spectrum based on Complexes
  if(flg==1)
    {
      prin_input("out/test.inp",lab,tag,input,limHSP,limSub);   //print input matrix
      prin_spec("out/test.out",lab,tag,data,lines,2);           //print output spectrum
    }
  return chi2/lines;
}




double gridHSPonly(string chiout,char* tempy,int num,double* params,int limHSP,int limSub,double Hx0min,double Hx0max,int Hx0pts,double Hsigmin,double Hsigmax,int Hsigpts,double Sx0,double Ssig,double adduction,double Zwidth,double shspMass,double subMass,int testmax,int minZ,int maxZ,double thresh,double* data,int lines,double MRes)
{
  double chi2run[Hx0pts*Hsigpts*5];
  chiout.append(tempy,num);
  FILE *fp;
  fp=fopen(chiout.c_str(),"w");
  fclose(fp);
  for (int i=0;i<Hx0pts;i++)
    {
      for (int j=0;j<Hsigpts;j++)
	{	

	  double Hx0=ceil(Hx0min*1.0+((i*1.0)/(Hx0pts*1.0-1.0))*(Hx0max*1.0-Hx0min*1.0))*1.0;      //#Hsp x0
	  double Hsig=(Hsigmin*1.0+((j*1.0)/(Hsigpts*1.0-1.0))*(Hsigmax*1.0-Hsigmin*1.0))*1.0;  //#Hsp width    
	  double par[10];par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;                     //parameters for chap=0
 	  double chi2=run_calc(0,0,tempy,num,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
	  //cout << "Run with: Hx0 " << Hx0 << " Hsig " << Hsig << " Sx0 " << Sx0 << " Ssig " << Ssig <<  " chi2/dof " << chi2 << endl;
	  
	  chi2run[i+j*Hx0pts + Hx0pts*Hsigpts*0]=Hx0;
	  chi2run[i+j*Hx0pts + Hx0pts*Hsigpts*1]=Hsig;
	  chi2run[i+j*Hx0pts + Hx0pts*Hsigpts*2]=Sx0;
	  chi2run[i+j*Hx0pts + Hx0pts*Hsigpts*3]=Ssig;
	  chi2run[i+j*Hx0pts + Hx0pts*Hsigpts*4]=chi2;
	  fp=fopen(chiout.c_str(),"a");
	  fprintf(fp,"%f\t%f\t%f\t%f\t%e\n",Hx0,Hsig,Sx0,Ssig,chi2);
	  fclose(fp);
	}
      fp=fopen(chiout.c_str(),"a");
      fprintf(fp,"\n");
      fclose(fp);
    }
  
  double chi2=0;
  {
    //Find the lowest chi^2 and print the solution
    int imin=findmin(chi2run,Hx0pts*Hsigpts,4);
    double Hx0 =chi2run[imin+ Hx0pts*Hsigpts*0];
    double Hsig=chi2run[imin+ Hx0pts*Hsigpts*1];
    Sx0        =chi2run[imin+ Hx0pts*Hsigpts*2];
    Ssig       =chi2run[imin+ Hx0pts*Hsigpts*3];
    double par[10];par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;                     //parameters for chap=0
    chi2=run_calc(1,0,tempy,num,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
    params[0]=Hx0;params[1]=Hsig;params[2]=Sx0;params[3]=Ssig;params[4]=chi2;
  }

  return chi2;
}









#include "disty_fit.c"


void runchap(string infile,char* outfile,int chap,double* x_init,struct mass massy,double subMass)
{

  int lines=countlines_string(infile);
  double data[(lines+1)*3];init_array(data,3*(lines+1));    
  struct data spec={chap,0,data,lines,massy.limHSP,massy.limSub,massy.adduction,massy.Zwidth,massy.shspMass,subMass,massy.testmax,massy.minZ,massy.maxZ,massy.thresh,massy.MRes};
  readfile(data,lines,infile);            //load in input data
  prin_spec(infile,".out",4,data,lines,1);//print out input spectrum
  
  if(chap==0)
       fitty(outfile,spec,x_init,4);

  if(chap==1)    
      fitty(outfile,spec,x_init,5);

  return;
}





int main()
{

  //Mass spec parameters
  double thresh = 0.001;   // Abundance threshold for inclusion of a charge state.
  double Zwidth = 1.3;     // Width of the charge state distribution, value from Justin = 1.1352
  double minMZ  = 7500.0;  // Minimum m/z to plot.
  double maxMZ  = 13000;   // Maximum m/z to plot.
  double MRes = 400;       // Mass spectral resolving power (FWHM/(m/z))
  double adduction = 1.001;// Adduction
  int minZ = 1;            // Lowest possible charge state.
  int maxZ = 200;          // Maximum possible charge state...provide a comfortable margin.

  //Parameters for fitting
  int limHSP= 60;     //highest HSP subunit to input
  int limSub= 10;     //largest substract in complex
  int testmax=5000;   //maximum number of trial complexes above threshold

  //Details on the system of interest
  double shspMass = 17985;    //mass of sHSP
  struct mass massy={thresh,Zwidth,minMZ,maxMZ,MRes,adduction,minZ,maxZ,limHSP,limSub,testmax,shspMass};


  //chap=0, and the distribution is a 2D gaussian specified by 4 params
  //chap=1 and the distribution is a 2D gaussian with a skew in the HSP dirction, 5 params
  int chap=1;

  
  // 0 - 1:1 Luciferase
  // 1 - 1:1 CS
  // 2 - 1:1 MDH
  // 3 - 1:0.1 MDH
  int FIT=3;
  
  
  if(FIT==0)
    {//fit Luc1:1 to skewed gaussian
      string infile="raw/Q2_FS_020908_6_Luc1_1.fix";
      double subMass = 60803; //mass of subunit (Luciferase)
      double x_init[5]= { 25.0, 5.0, 1.0,0.5,10.0 };  
      runchap(infile,"goLuc1",1,x_init,massy,subMass);
    }
  if(FIT==1)
    {//fit CS1:1 to skewed gaussian
      string infile="raw/Q2_FS_060409_10_CS1_1.fix";
      double subMass = 49072;  //mass of subunit (CS)
      double x_init[5]= { 25.0, 5.0, 2.0,0.5,10.0 };  
      runchap(infile,"goCSS1",1,x_init,massy,subMass);
    }
  if(FIT==2)
    {//fit MDH:1 to skewed gaussian
      string infile="raw/Q2_FS_260808_8_MDH1_1.fix";
      double subMass = 33100; //mass of subunit (MDH)
      double x_init[5]= { 25.0, 5.0, 2.0,0.5,10.0 };  
      runchap(infile,"goMDH1",1,x_init,massy,subMass);
    }
  if(FIT==3)
    {//fit MDH:1 to skewed gaussian
      string infile="raw/Q2_FS_270808_5_MDH1_0.1.fix";
      double subMass = 33100; //mass of subunit (MDH)
      double x_init[5]= { 25.0, 5.0, 2.0,0.5,10.0 };  
      runchap(infile,"goMDH2",1,x_init,massy,subMass);
    }
    





    
   //     double x_init[4]= { 25.0, 5.0, 1.0,0.5 };  
  /*  {
    int lines=countlines_string(infile4);
    double data[(lines+1)*3];init_array(data,3*(lines+1));    
    struct data spec={chap,0,data,lines,limHSP,limSub,adduction,Zwidth,shspMass,subMass_MDH,testmax,minZ,maxZ,thresh,MRes};
    readfile(data,lines,infile4);            //load in input data
    prin_spec(infile4,".out",4,data,lines,1);//print out input spectrum
    runchap(infile4,"goMDH2",1,spec);

    }*/
  
  



  // calc distribution and chi^2 for one condition
  /*  {
    double Hx0 = 27.6363;
    double Hsig = 6.666;
    double Sx0   = 0.75336;  
    double Ssig  = 0.36841;     //Sub width
    double alpha= -0.7372;
    double par[10];par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;par[4]=alpha;       //parameters for chap=0    
    double chi2=run_calc(1,1,".go",3,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
    cout << "chi2/dof " << chi2 << endl;
  }
*/  
  
  /*
   {//Run 2d - fixed substrate but variable Hsps
    double Hx0max= 55.0;
    double Hx0min= 15.0;
    int Hx0pts= 10;
    double Hsigmax= 10.1;
    double Hsigmin= 0.1;
    int Hsigpts= 10;
    
    double Sx0   = 1.0;  
    double Ssig  = 0.001;     //Sub width

    
    double params[5];
    string chiout = "out/chi.out";
    char* tempy="_go"
    gridHSPonly(chiout,tempy,3,params,limHSP,limSub,Hx0min,Hx0max,Hx0pts,Hsigmin,Hsigmax,Hsigpts,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
  }
  */  

  /*
  {//Run 2D fit but over a range of substract conditions
    double Hx0max= 45.0;
    double Hx0min= 15.0;
    int Hx0pts= 20;
    double Hsigmax= 10.1;
    double Hsigmin= 0.1;
    int Hsigpts= 20;

    for (int i=0;i<5;i++)
      {
	cout << "Running " << i+1 << " of 5" <<endl; 
	double Sx0   = i*1.0;  
	double Ssig  = 0.001;     //Sub width
	double params[5];
	string chiout = "out/chi.out";
	char tempy[10];sprintf(tempy,"%i",i);
	gridHSPonly(chiout,tempy,3,params,limHSP,limSub,Hx0min,Hx0max,Hx0pts,Hsigmin,Hsigmax,Hsigpts,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
      }
      }*/




  /*  
  {//Run brute force fit with four free parameters
    double Hx0max= 45.0;
    double Hx0min= 15.0;
    int Hx0pts= 20;

    double Hsigmax= 10.1;
    double Hsigmin= 0.1;
    int Hsigpts= 20;

    double Sx0max= 3.0;
    double Sx0min= 0.0;
    int Sx0pts= 20;

    double Ssigmax= 3.0;
    double Ssigmin= 0.01;
    int Ssigpts= 20;

    
    FILE *fp;
    string chioutL="out/chiLong.out";
    double chi2runL[Sx0pts*Ssigpts*5];
    for (int i=0;i<Sx0pts;i++)
      {
	for (int j=0;j<Ssigpts;j++)
	  {
	    cout << "Running i " << i+1 << " of " << Sx0pts << " and j " << j+1 << " of " << Ssigpts <<endl; 
	    double Sx0   = Sx0min+(Sx0max-Sx0min)*i/(Sx0pts-1.0);
	    double Ssig  = Ssigmin+(Ssigmax-Ssigmin)*j/(Ssigpts-1.0);

	    double params[5];
	    string chiout = "out/chi.out";
	    char tempy[10];sprintf(tempy,"%i-%i",i,j);
	    double chi2=gridHSPonly(chiout,tempy,6,params,limHSP,limSub,Hx0min,Hx0max,Hx0pts,Hsigmin,Hsigmax,Hsigpts,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
	    
	    chi2runL[i+j*Sx0pts + Sx0pts*Ssigpts*0]=params[0];
	    chi2runL[i+j*Sx0pts + Sx0pts*Ssigpts*1]=params[1];
	    chi2runL[i+j*Sx0pts + Sx0pts*Ssigpts*2]=Sx0;
	    chi2runL[i+j*Sx0pts + Sx0pts*Ssigpts*3]=Ssig;
	    chi2runL[i+j*Sx0pts + Sx0pts*Ssigpts*4]=chi2;
	    
	    fp=fopen(chioutL.c_str(),"a");
	    fprintf(fp,"%f\t%f\t%f\t%f\t%e\n",params[0],params[1],Sx0,Ssig,chi2);
	    fclose(fp);
	    
	  }
	fp=fopen(chioutL.c_str(),"a");
	fprintf(fp,"\n");
	fclose(fp);
      }
  
    {
      //Find the lowest chi^2 and print the solution
      int imin=findmin(chi2runL,Sx0pts*Ssigpts,4);
      double Hx0 =chi2runL[imin+ Sx0pts*Ssigpts*0];
      double Hsig=chi2runL[imin+ Sx0pts*Ssigpts*1];
      double Sx0 =chi2runL[imin+ Sx0pts*Ssigpts*2];
      double Ssig=chi2runL[imin+ Sx0pts*Ssigpts*3];
      double chi2=run_calc(1,"global",6,limHSP,limSub,Hx0,Hsig,Sx0,Ssig,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,data,lines,MRes);
      //params[0]=Hx0;params[1]=Hsig;params[2]=Sx0;params[3]=Ssig;params[4]=chi2;
    }
    }*/
  
  


  return 0;
}
