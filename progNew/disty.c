
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
 

 string inputfile= argv[1];

 initfile("figs/spectrafit.gp");   //initialise gnuplot file
 initfile("figs/make_summary.com");//initialise summary file

 anal inst;

 inst.ParseInpFile(inputfile);
 for(int i=0;i<inst.files.size();++i)
   {
     inst.files[i].make_summary_init("figs/make_summary.com");  //write first line of summary file
     inst.files[i].readfile();  //read datafile
     if(inst.files[i].input_read)
       inst.files[i].read_input();
     else
       inst.files[i].make_input();
     inst.RunProtocol(i);       //execute protocol
     inst.files[i].make_summary_end("figs/make_summary.com");
   }
 return 0;
/*


 inst.ReadInput(inputfile); //parse the input file

 FILE*fp;
 fp=fopen("output.txt","w");
 fclose(fp);
 for(int i=0;i<inst.files.size();++i)
   {
     inst.files[i].make_summary_init("figs/make_summary.com");  //write first line of summary file
     
     inst.files[i].adduction=1E-3;
     inst.files[i].thresh   =0.00001;  // Abundance threshold for inclusion of a charge state.
     inst.files[i].ResFudge =1E-12; // linear term for m/z peak width
     inst.files[i].minZ     =1;  // Lowest possible charge state.
     inst.files[i].maxZ     =200;   // Maximum possible charge state...provide a comfortable margin.
     inst.files[i].limHSP   =30; //highest HSP subunit to input
     inst.files[i].limSub   =35;  //largest substract in complex
     inst.files[i].testmax  =50000; //maximum number of trial complexes above threshold (error message comes out if this is too small)
     
     inst.files[i].subMass    =22782;
     inst.files[i].shspMass   =20159;
     //inst.files[i].mode       =10;
     
     inst.files[i].type="Orbi";
     inst.files[i].readfile();  //read datafile

     /*
     inst.files[i].fitpar.adduction_flg=1;
     inst.files[i].fitpar.MRes_flg=1;
     inst.files[i].fitpar.ResFudge_flg=1;
     inst.files[i].fitpar.Zfudge_flg=1;
     inst.files[i].fitpar.Zwidth_flg=1;
     inst.files[i].fitpar.dist_flg=0;
     
     inst.files[i].readfile();  //read datafile
     inst.files[i].read_input("raw/equil.txt"); //read a distribution from text file
     inst.runchap("sim",i,1);   //run minimisation (j specifies skew or not)
     inst.runchap("fitty",i,1); //run minimisation (j specifies skew or not)
     inst.runchap("jiggle",i,1); //run minimisation (j specifies skew or not)

     */     /*

     //inst.files[i].chap=10;
     
     

     
     inst.files[i].fitpar.adduction_flg=1;
     inst.files[i].fitpar.MRes_flg=1;
     inst.files[i].fitpar.ResFudge_flg=1;
     inst.files[i].fitpar.Zfudge_flg=1;
     inst.files[i].fitpar.Zwidth_flg=1;
     inst.files[i].fitpar.dist_flg=1;

     inst.files[i].limHSP   =120; //highest HSP subunit to input
     inst.files[i].limSub   =1;  //largest substract in complex
     inst.files[i].dim=1;

     inst.files[i].AddDist("alphaB");
     inst.files[i].distList[0].pars[0]=40;
     inst.files[i].distList[0].pars[1]=1;
     inst.files[i].distList[0].fit_flg=0;

     inst.files[i].fitpar.dist_flg=1;

     inst.runchap("sim",i,1);   //run minimisation (j specifies skew or not)
     inst.runchap("fitty",i,1); //run minimisation (j specifies skew or not)


     inst.files[i].distList[0].fit_flg=1;
     inst.runchap("fitty",i,1); //run minimisation (j specifies skew or not)

     inst.runchap("jiggle",i,1); //run minimisation (j specifies skew or not)




     inst.files[i].make_summary_end("figs/make_summary.com");

     //for gridsearch
   */  /*
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
     *//*
     
     
    */ /*fp=fopen("output.txt","a");
       fprintf(fp,"%s\t%s\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",raw[i].c_str(),ident[i].c_str(),mode[i],chi2,massy.MRes,massy.Zfudge,massy.adduction,shspMass[i],x_init[0],x_init[1] );
       fclose(fp);*//*
     
     //FINISH RUN FOR GIVEN FILE
     
   }//FINISHED
*/
 return 0;
}


