/*################################################################*/
//# Functions for disty.c
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <complex>
#include <gsl/gsl_sf_erf.h>

using namespace std;

//initialise a new file
void initfile(char *infile)
{
  FILE *fp;
  fp=fopen(infile,"w");
  fclose(fp); 
}
//initialise a new file
void initfile(string infile)
{
  FILE *fp;
  fp=fopen(infile.c_str(),"w");
  fclose(fp); 
}


//stick a file in a vector of vectors
vector<vector<string> > MakeFileVec(string file_name)
{
  
  ifstream in(file_name.c_str());
  vector<vector<string> > infile;
  string line;
  while(getline(in,line)){
    istringstream iss(line);
    vector<string> tokens;
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter(tokens));
    infile.push_back(tokens);
  }
  in.close();
  return infile;
}  




class mass
{
 public:

  FILE *fp;

  string raw,ident,identTag; //data file, identifier, progress identifier
  double Zwidth,MRes,Zfudge,adduction,ResFudge; //mass spec parameters
  double shspMass,subMass; //massess of species i and j
  double minMZ,maxMZ;  //max and min for raw data
  int smooth;   //smoothing factor for raw data
  string mode;  //sim fitty or jiggler
  double chi2; //chi2
  double cWid; //MZ width cutoff when evalulating spectra (MZrange: centroid+/-cSTD*cWid)

  string type;  //spectrometer type (QTof, Orbi)
  string outfile;  //output file
  string input_file; //input distribution to be read in
  int input_read=0;  //read in a distribution from a text file

  double thresh;  //lower limit on ion weight to be considered
  int minZ;      //lowest charge state
  int maxZ;     //highest charge state
  int limHSP;   //biggest species i
  int limSub;   //biggest species j

  int trial=0;  //progress indicator (each sim/fitty/jiggle, this goes up by 1)
  int lines;  //number of data points in smoothed data

  double *err;      //fitting errors
  double *datatmp;  //full length raw data
  double *data;     //smoothed raw data
  double *par;      //fitting parameters
  double *input;    //input ij array
  double *inputprop;   //slice of ij array
  double *input_comp;  //subset of ij array
  double *candidateZ;  // Array for all possible charge states.
    
  double rawmax;  //maximum value of raw intensity
  int dim=1; //1D or 2D for client distributions

  //parameters used in the various models
  //WITH CLIENTS
  //0  2D Gaussian  Hx0 Hsig Sx0 Ssig
  //8  2D Gaussian (grid) H0 Hsig Sx0 Ssig Alpha
  //1  Skew Gaussian
  //9  Skew Gaussian (grid)
  //2  AlphaB  (rat,mult)
  //3  Skewed in both dimensions  Hx0 Hsig alpha Sx0 Ssig alpha2
  //4  Skewed Gaussian, client 0 free  Hx01 Hsig1 alpha1 rel  Hx0 Hsig Sx0 Ssig alpha	    
  //5  2D Gaussian (separate slices) relA Hx0A HsigA	    
  //6  Skew Gaussian (separate slices) relA Hx0A HsigA alphaA	    
  //7  Poisson (no clients)
  //10 Free.
  //NO CLIENTS

  class free //container for free species
  {
  public:
    int i;     //size of species 
    int j=0;   //size of species
    double conc,err=0;  //conc and err
    int fit_flg=0;      //fitting?
    void CalcDist(double *input,int limHSP,int limSub) //put concentratio into input
    {
      input[i+limHSP*j]=conc;
    }
    void ShowPars() //show the current state of this species
    {
      if(fit_flg)
	{
	  printf("Oligo: i %5i j %5i conc %.3f  +/-  %.3f\n",i,j,conc,err); 
	}
      else
	{
	  printf("Oligo: i %5i j %5i conc %.3f\n",i,j,conc); 
	}
    }

  };
  class dist //container for distributions
  {
  public:
    string type;   //distribution type
    double *pars;  //parameters
    double *errs;  //errors
    int fit_flg=0; //fitting?
    int p;      //total parameters associated with distribution
    int dim;    //equal to 1 or 2 (no client or client)
    double Norm(double x,double b,double c){
      return exp(-pow(x-b,2.0)/(2*c*c));
      //return exp(-(x-b)*(x-b)/(2*c*c));
    }
    double SkewNorm(double x,double b,double c,double d){
      return exp(-pow(x-b,2.0)/(2*c*c))* 0.5*(1+gsl_sf_erf(d*(x-b)/(2*c)/sqrt(2)));
    }
    int findmax(double *array,int N,int j)
    {
      double  max=array[0+j*N]*1.0;
      int imax=0;
      for (int i=0;i<N;i++)
	if(array[i+j*N]*1.0>max)
	  {
	    max=array[i+j*N];
	    imax=i;
	  }
      return imax;
    }

    void SetDist()  //set distribution parameters
    {
      switch(dim){
      case 1:
	if(type=="Gaussian")
	  p=2;
	if(type=="SkewGaussian")
	  p=3;
	if(type=="alphaB")
	  p=2;
	break;
      case 2:
	if(type=="Gaussian")
	  p=4;
	if(type=="SkewGaussian")
	  p=5;
	if(type=="DoubleSkew")
	  p=6;
	if(type=="ClientSkew")
	  p=9;
	break;
      }
      pars=new double[p]; //set memory
      errs=new double[p]; //set memory
    }
    void CalcDist(double *input,int limHSP,int limSub) //calculate distribution and pour into input
    {
      switch(dim){
      case 1:
	if(type=="Gaussian")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      if(j==0)
		input[i+j*limHSP]+=(Norm(i*1.0,pars[0],pars[1]));
	if(type=="SkewGaussian")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      if(j==0)
		input[i+j*limHSP]+=(SkewNorm((i)*1.0,pars[0],pars[1],pars[2]));
	if(type=="alphaB")
	  {
	    double *tmp;
	    tmp=new double[limHSP*limSub];
	    for(int i=0;i<limHSP;++i)
	      for(int j=0;j<limSub;++j)
		if(j==0)
		  {
		    //i+1 is the number of sHSPs
		    if(i==0)
		      tmp[i+j*limHSP]=1;
		    //calculate distribution
		    else
		      {
			if((i+1)%2==0) //then  even
			  tmp[i+j*limHSP]=fabs(tmp[i-1+j*limHSP]*(fabs(pars[0])/(i+1)));
			if((i+1)%2==1) //then odd
			  tmp[i+j*limHSP]=fabs(tmp[i-1+j*limHSP]*(fabs(pars[0])/(i+1-1+fabs(pars[1]))));
		      }
		  }
	    //cout << pars[0] << " " << pars[1] << endl;

	    int imax=findmax(tmp,limHSP*limSub,0);
	    double max=tmp[imax];
	    for(int i=0;i<limHSP*limSub;i++)
	      input[i]+=tmp[i]/max;
	  }
	break;
      case 2:
	if(type=="Gaussian")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      if(j==0)
		input[i+j*limHSP]+=(Norm(i*1.0,pars[0],pars[1]))*(Norm(j*1.0,pars[2],pars[3]));
	if(type=="SkewGaussian")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      if(j==0)
		input[i+j*limHSP]+=(SkewNorm((i)*1.0,pars[0],pars[1],pars[2])*Norm(j*1.0,pars[3],pars[4]) );
	if(type=="DoubleSkew")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      input[i+j*limHSP]+=(SkewNorm(j*1.0,pars[3],pars[4],pars[5]) * SkewNorm(i*1.0,pars[0],pars[1],pars[2]) );
	if(type=="ClientSkew")
	  for(int i=0;i<limHSP;++i)
	    for(int j=0;j<limSub;++j)
	      if(j==0)
		{input[i+j*limHSP]+= fabs(pars[3]) * (SkewNorm((i)*1.0,pars[0],pars[1],pars[2])  );}
	      else
		{input[i+j*limHSP]+=(Norm(j*1.0,pars[7],pars[8]) * SkewNorm((i)*1.0,pars[4],pars[5],pars[6]) );}
	break;      
      }
    }
    void ShowPar(string lab,int i) //print parameters to screen
    {
      if(fit_flg)
	printf ("%s = %.5f +/- %.5f\n", lab.c_str(),pars[i], errs[i]);
      else
	printf ("%s = %.5f\n", lab.c_str(),pars[i]);
    }

    void ShowPars() //print parameters to screen
    {
      printf("Distribution: %s\n",type.c_str());
      switch(dim){
      case 1:
	if(type=="Gaussian")
	  {
	    ShowPar("Hx0 ",0);
	    ShowPar("Hsig",1);
	  }
	if(type=="SkewGaussian")
	  {
	    ShowPar("Hx0 ",0);
	    ShowPar("Hsig",1);
	    ShowPar("alpha",2);
	  }
	if(type=="alphaB")
	  {
	    ShowPar("rat ",0);
	    ShowPar("mult",1);
	  }
	break;
      case 2:
	if(type=="Gaussian")
	  {
	    ShowPar("Hx0 ",0);
	    ShowPar("Hsig",1);
	    ShowPar("Sx0 ",2);
	    ShowPar("Ssig",3);
	  }
	if(type=="SkewGaussian")
	  {
	    ShowPar("Hx0  ",0);
	    ShowPar("Hsig ",1);
	    ShowPar("alpha",2);
	    ShowPar("Sx0  ",3);
	    ShowPar("Ssig ",4);
	  }
	if(type=="DoubleSkew")
	  {
	    ShowPar("Hx0   ",0);
	    ShowPar("Hsig  ",1);
	    ShowPar("Halpha",2);
	    ShowPar("Sx0   ",3);
	    ShowPar("Ssig  ",4);
	    ShowPar("Salpha",5);
	  }
	if(type=="ClientSkew")
	  {
	    ShowPar("Hx0     ",0);
	    ShowPar("Hsig0   ",1);
	    ShowPar("Halpha0 ",2);
	    ShowPar("rel0    ",3);
	    ShowPar("Hx      ",4);
	    ShowPar("Hsig    ",5);
	    ShowPar("alpha   ",6);
	    ShowPar("Sx      ",7);
	    ShowPar("Ssig    ",8);
	  }
	break;      
      }
    }


  };

  //to add a new optional fitting parameter:
  //add it here. Then search for all intances of the 'last' parameter
  //such as shspMas_flg. Then add a new entry for the new parameter whenever
  //you encounter it in the code.
  class fitPars  //fitting parameter class
  {
  public:
    int dist_flg;       //fit distributions?
    int adduction_flg;  //fit adduction?
    int MRes_flg;       //fit MRes?
    int ResFudge_flg;   //fit ResFudge?
    int Zfudge_flg;     //fit Zfudge?
    int Zwidth_flg;     //fit Zwidth?
    int shspMass_flg;   //fit shspMass?
    int distparam;      //number of distribution parameters
    int specparam;      //number of spectral parametres
    int p;  //total parameters
    
    void GetSpecParams() //calculate how many spectral parameters to worry about
    {
      specparam=adduction_flg+Zfudge_flg+MRes_flg+ResFudge_flg+Zwidth_flg+shspMass_flg; //add on spectral params
      p=specparam;
      if(dist_flg) //if optimising distribution, then add parameters
	p+=distparam;
      cout << "Spectral parameters to fit: " << specparam << endl;
      cout <<"Total parameters to fit:  " << p << endl;
    }
    void Reset() //turn all fitting parameters off
    {
      adduction_flg=0;
      MRes_flg=0;
      ResFudge_flg=0;
      Zfudge_flg=0;
      Zwidth_flg=0;
      shspMass_flg=0;
      dist_flg=0;
    }
    //set fitting parameters from vector pars
    void SetPars(vector<string> pars,vector<dist> &distList,vector<free> &freeList)
    {
      Reset(); //reset flags for fitting parameters
      
      string dist="dist";
      string free="free";
      
      for(int k=0;k<pars.size();++k)
	{
	  //cout << pars[k] << " " << pars[k].find(dist) << endl;
	  if(pars[k]=="adduction")
	    adduction_flg=1;
	  else if(pars[k]=="Zfudge")
	    Zfudge_flg=1;
	  else if(pars[k]=="Zwidth")
	    Zwidth_flg=1;
	  else if(pars[k]=="MRes")
	    MRes_flg=1;
	  else if(pars[k]=="ResFudge")
	    ResFudge_flg=1;
	  else if(pars[k]=="shspMass")
	    shspMass_flg=1;
	  else if(pars[k].find(dist)==0)
	    {
	      //char * v=pars[k][pars[k].size()-1];
	      //char v=pars[k][0];
	      //int ind=atoi(&v);
	      for(int j=0;j<distList.size();++j)
		distList[j].fit_flg=1;
	      dist_flg=1;
	    }
	  else if(pars[k].find(free)==0)
	    {
	      //char v=pars[k][0];
	      //char * v=pars[k][pars[k].size()-1];
	      //int ind=atoi(&v);
	      //freeList[ind].fit_flg=1;
	      for(int j=0;j<freeList.size();++j)
		freeList[j].fit_flg=1;
	      dist_flg=1;
	    }
	  else
	    {
	      cout << "unrecognised fitting flag:" << endl;
	      cout << pars[k] <<endl;
	      exit(100);
	    }
	}
    }
    void ShowFlags() //print parameters to fit:
    {
      cout << "Parameters to fit: " << endl;
      if(dist_flg)
	cout << " dist";
      if(adduction_flg)
	cout << " adduction";
      if(Zfudge_flg)
	cout << " Zfudge";
      if(Zwidth_flg)
	cout << " Zwidth";
      if(MRes_flg)
	cout << " MRes";
      if(ResFudge_flg)
	cout << " ResFudge";
      if(shspMass_flg)
	cout << " shspMass";
      cout << endl;

    }
  };
  

  fitPars fitpar; //class to store fitting flags

  vector<free> freeList; //vector to store individual oligomeric species
  vector<dist> distList; //vector to store all required distributions

  void AddFree(int i,int j,double conc)  //add a free species to freeList
    {
      free entry;
      entry.i=i;
      entry.j=j;
      entry.conc=conc;
      freeList.push_back(entry);
    }
  void AddDist(string type) //add a distribution to distList
    {
      dist entry;
      entry.type=type;
      entry.dim=dim;
      entry.SetDist(); //create memory required
      distList.push_back(entry); //store
    }
  
  void make_input() //calculate current input distribution
  {

    init_array(input,limHSP*limSub);    //first, initialise input array    
    //for(int i=0;i<distList.size();++i)
    for(std::vector<dist>::iterator it = distList.begin(); it != distList.end(); ++it) 
      it->CalcDist(input,limHSP,limSub); //add on each distribution
    //for(int i=0;i<freeList.size();++i)
    for(std::vector<free>::iterator it = freeList.begin(); it != freeList.end(); ++it) 
      it->CalcDist(input,limHSP,limSub); //add on each distribution
  }
  
  //initialise array
  void init_array(double *array,int N)
  {
    memset(array, 0, sizeof(double)*N); 
    return;
  }
  
  /*
  vector<string> parList; //string that will contain entries from line input
  //(deppreciated way of reading in files)
  void AddLine(vector<string> linec) //edit line from master input file.
  {
    raw=linec[0];
    ident=linec[1];
    
    smooth=atoi(linec[2].c_str());

    shspMass=atof(linec[3].c_str());
    minMZ=atof(linec[4].c_str());
    maxMZ=atof(linec[5].c_str());
    //chap=atoi(linec[6].c_str());
    //xinit0=atof(linec[7].c_str());
    //xinit1=atof(linec[8].c_str());
    Zwidth=atof(linec[9].c_str());
    MRes=atof(linec[10].c_str());
    Zfudge=atof(linec[11].c_str());
    //max0=atof(linec[12].c_str());
    //min0=atof(linec[13].c_str());

    parList.push_back("raw");
    parList.push_back("ident");
    parList.push_back("smooth");
    parList.push_back("shspMass");
    parList.push_back("minMZ");
    parList.push_back("maxMZ");
    parList.push_back("chap");
    parList.push_back("xinit0");
    parList.push_back("xinit1");
    parList.push_back("Zwidth");
    parList.push_back("MRes");
    parList.push_back("Zfudge");
    parList.push_back("max0");
    parList.push_back("min0");

    for(int i=0;i<14;++i)
      {
	cout << parList[i] << " " << linec[i] << endl;
      }
  }
  */

  //read in raw data file, smooth and normalise
  int readfile()
  {
    vector<vector<string> > indat;
    indat=MakeFileVec(raw);
    lines=0;
    for(int i=0;i<indat.size();++i)
      	if(indat[i].size()==2)
	  lines+=1;
    
    cout << "Number of data lines in " << raw << ": " << lines << endl;

    input=new double[limHSP*limSub];       //input distribution
    input_comp=new double[limHSP*limSub];  //subset of the input distribution
    inputprop=new double[limSub];         //slice from the input distribution
    par=new double[limHSP*limSub+10];  //max parameters to be fitted
    err=new double[limHSP*limSub+10];  //max errors to be obtained
    candidateZ=new double[(maxZ+1-minZ)];  // Array for all possible charge states.

    datatmp=new double [lines*2];    //temp data: unsmoothed
    
    int ii=0;
    for(int i=0;i<indat.size();++i)
      	if(indat[i].size()==2)
	  {
	    datatmp[ii]=atof(indat[i][0].c_str());      //mz values
	    datatmp[ii+lines]=atof(indat[i][1].c_str()); //intensity values
	    ii++;
	  }
    cout << "Smoothing..." << endl;
    
    int NewLines=lines/smooth;
    data=new double[NewLines*3];
    for(int i=0;i<NewLines;i++)
      {
	double exval=0.0;
	double eyval=0.0;
	double cnt=0.0;
	for (int j=0;j<smooth;j++)
	  {
	    int ieff=j+i*smooth;
	    exval+=datatmp[ieff];
	    eyval+=datatmp[ieff+lines];
	    cnt+=1;
	  }
	data[i]=exval/cnt;
	data[i+NewLines]=eyval/cnt;
      }
    cout << "Final number of points in datafile: " << NewLines << endl;
    lines=NewLines;
    cout << "Normalising" << endl;
    normdata();

    cout << "writing.."<< endl;
    prin_spec(raw);  //print file with adjusted input spectrum

    rawmax=data[findmax(data,lines,1)+1*lines]; //get maximum of raw data
    
  }

  //normalise highest point on input spectrum to 100
  void normdata()
  {
    int imax=findmax(data,lines,1);
    double max=data[imax+lines*1];
    for(int i=0;i<lines;i++)
      data[i+lines*1]=data[i+lines*1]/max*100;
    return;
  }   


  int findmax(double *array,int N,int j)
  {
    double  max=array[0+j*N]*1.0;
    int imax=0;
    for (int i=0;i<N;i++)
      if(array[i+j*N]*1.0>max)
	{
	  max=array[i+j*N];
	  imax=i;
	}
    return imax;
  }
  
  
  int findmin(double *array,int N,int j)
  {
    double  min=array[0+j*N]*1.0;
    int imin=0;
    for (int i=0;i<N;i++)
      if(array[i+j*N]*1.0<min)
	{
	  min=array[i+j*N];
	  imin=i;
	}
    return imin;
  }

  /**************MS FUNCTIONS FOR CALCULATING SPECTRUM***************************/
  
  //# Controls centroid of charge state distributions.
  //#averageZ = x* tempM ** y
  double avgZ(double tempM,double Zfudge){
    return 0.0467 * pow(tempM,0.533)+Zfudge; //Values from Justin are 0.0467, 0.533;  
  }
  
  //# Detector response
  //#DetectionEfficient =  x*(1-exp(-y*(tempZ*9.1/tempM)**z))
  //only apply for QTof
  double DetectionEfficiency(double tempM,double tempZ){
      return 60*(1-exp(-1620 * pow((tempZ*9.1/tempM),1.75)));
  }
  
  //# Evaluate Gaussian function
  //# See http://en.wikipedia.org/wiki/Gaussian_function
  //# x = current m/z
  //# b = centroid m/z
  //# c = width = FWHM/(2*sqrt(2*ln(2)))
  double Norm(double x,double b,double c){
    return exp(-pow(x-b,2.0)/(2*c*c));
  }
  
  
  double SkewNorm(double x,double b,double c,double d){
    return exp(-pow(x-b,2.0)/(2*c*c))* 0.5*(1+gsl_sf_erf(d*(x-b)/(2*c)/sqrt(2)));
  }
  

  void eval_spec(double *inp)
  {
    reset_sim();   //reset calculated spectrum    

    for (int i=0;i<limHSP;i++)   //for all i to be considered...
      for (int j=0;j<limSub;j++)  //for all j to be considered...
	{
	  double mass=1.0*(((i)*shspMass*1.0)+((j)*subMass*1.0)); //sequence mass
	  mass=mass+fabs(adduction)*pow(mass,0.760);  //adducted mass
	  double avZ=avgZ(mass,Zfudge); //average charge state distribution
	  double conc=inp[i+limHSP*j];  //concentration from submitted array

	  double TotalWeight = 0;              // Intensity for all charge states.
	  for(int j=0;j<maxZ+1-minZ;j++)       // loop over all charge states:
	    {           
	      double CurrentWeight=Norm(minZ+j,avZ,Zwidth); //get gaussian height of charge state
	      candidateZ[j]=CurrentWeight;
	      TotalWeight += CurrentWeight;   	//Updates the total weight for all z of this component.
	    }
	  
	  for(int j=0;j<maxZ+1-minZ;j++){  //  loop over all charge states
	    // Normalizes intensities for charge state peaks based on the abundance from the input files
	    double CorrectWeight=conc*candidateZ[j]/TotalWeight;
	    if (CorrectWeight >= thresh)
	      {
		if(type=="QTof")
		  CorrectWeight=CorrectWeight*DetectionEfficiency(mass,candidateZ[j]);
		
		// Define bounds for current charge state
		double centroid = (mass )/(minZ+j); //mass over charge
		double cSTD     = centroid*((1/MRes*1.0) +1.0*((fabs(ResFudge)/100000000))*centroid);
		double lowMZ    = centroid - (cWid * cSTD);  //lowest MZ to conside
		double highMZ   = centroid + (cWid * cSTD);  //highest MZ to consider
		//cout << i << "low " << lowMZ << " " << highMZ << " " << centroid << " " << Complex_test[i+testmax*2] << endl;
		// Evaluate contribution of current charge state
		for (int j=0;j<lines;j++)//for each point in spectrum
		  if(data[j]>=lowMZ && data[j]<=highMZ) //if it is within range...
		    data[j+lines*2]+= CorrectWeight*Norm(data[j],centroid,cSTD); //# Evaluate Gaussian function for (current m/z, centroid m/z, MS resolution)
	      }
	  }
	}
    
    double simmax=data[findmax(data,lines,2)+2*lines];
    //cout << "datamax " << rawmax << " simmax " << simmax << endl;
    for (int j=0;j<lines;j++) //normalise spectrum so that max point is same in both
      data[j+lines*2]=data[j+lines*2]/simmax*rawmax;
  }
  
  
  void reset_sim() //reset spectrum simulation
  {
    for (int i=0;i<lines;i++)
      data[i+lines*2]=0;
    return;
  }
  

  /*******************END OF MS FUNCTIONS **********************/  

  void eval_spec_norm(double *data,int lines,double sum,double *inputprop,int limSub,int i)
  {
    double sup=0;
    for (int j=0;j<lines;j++)
      sup+=data[j+lines*(2)];  
    
    double stuff=0;
    for (int j=0;j<limSub;j++)
      stuff+=inputprop[j];
    printf("    Fraction with %i clients bound: %.2f \n",i,inputprop[i]/stuff);
    
    for (int j=0;j<lines;j++)  
      data[j+lines*(2)]=    data[j+lines*(2)]*sum/sup*(inputprop[i]/stuff);  
    return;
  }
  
  
  void eval_spec_norm_i(double sum,double poppy)
  {
    double sup=0;
    for (int j=0;j<lines;j++)
      sup+=data[j+lines*(2)];  
    
    double stuff=0;
    for (int j=0;j<limSub;j++)
      stuff+=inputprop[j];
    //printf("    Guy with %i clients bound: %.2f \n",i,poppy/stuff);
    
    for (int j=0;j<lines;j++)  
      data[j+lines*(2)]=    data[j+lines*(2)]*sum/sup*(poppy/stuff);  
    return;
  }
  
  void CalcChi2()
  {
    chi2=0;
    for (int j=0;j<lines;j++)
      chi2+=pow((data[j+lines*(1)]-data[j+lines*(2)]),2.0); ///100;
  }
  


  //# Write only spectrum to output file
  void prin_spec(string outfile)
  {        
    fp=fopen((outfile+".out").c_str(),"w");
    for (int i=0;i<lines;i++)
      fprintf(fp,"%f\t%f\n",data[i],data[i+lines*1]);
    fclose(fp); 
    
    /*  fp=fopen(outfile+'.gp'.c_str(),'w');
	fprintf(fp,'set term post eps enh solid color 20\n');
	fprintf(fp,'set output \''+outfile+'.eps\'\n');
	fprintf(fp,'set title \''+outfile+' \'\n');
	fprintf(fp,'unset key\n');
	fprintf(fp,'set xlabel \'m/z\'\n');
	fprintf(fp,'plot \''+outfile+'\' u 1:2 w li');
	fclose(fp);*/
    return;
  }
  
  //# Write spectrum to output file
  //append a spectrum raw and simulated to an output file
  void prin_spec3(string outfile)
  {        
    fp=fopen((outfile+identTag).c_str(),"a");
    for (int i=0;i<lines;i++)
      fprintf(fp,"%f\t%f\t%f\n",data[i],data[i+lines*1],data[i+lines*2]);
    fprintf(fp,"\n\n");
    fclose(fp); 
  }
  

  //# Write spectrum raw and fitted to output file
  //make gnuplot script to plot this
  void prin_spec2(string outfile)
  {        
    fp=fopen((outfile+identTag).c_str(),"w");
    for (int i=0;i<lines;i++)
      fprintf(fp,"%f\t%f\t%f\n",data[i],data[i+lines*1],data[i+lines*2]);
    fclose(fp); 
    
    fp=fopen("figs/spectrafit.gp","a");
    fprintf(fp,"reset\n");
    fprintf(fp,"set term post eps enh solid color 20\n");
    fprintf(fp,"set output \'figs/test.out%s.eps\'\n",identTag.c_str());
    fprintf(fp,"set title \'Fitting spectrum: %s\'\n",identTag.c_str());
    fprintf(fp,"set label \"MZ per point: %.2f\" at graph 0.02,graph 0.95\n",(maxMZ-minMZ)/(1.0*lines));
    fprintf(fp,"set label \"chi2/dof:     %.2f\" at graph 0.02,graph 0.9\n",chi2/lines);
    fprintf(fp,"set label \"ave Error:    %.2f%s\" at graph 0.02,graph 0.85\n",sqrt(chi2/lines),"%");
    
    fprintf(fp,"set xlabel \'m/z\'\n");
    fprintf(fp,"set xrange [*:*]\n");
    fprintf(fp,"set yrange [-20:110]\n");
    fprintf(fp,"plot \'out/test.out%s\' u 1:2 ti 'raw' w li lt 1,\\\n",identTag.c_str());
    fprintf(fp,"\'\' u 1:3 ti 'fitted' w li lt 2,\\\n");
    fprintf(fp,"\'\' u 1:(($2-$3)-10) ti 'difference' w li lt 3\n");
    fprintf(fp,"unset label\n");
    fclose(fp);
    return;
  }


  //initialise comparison file and produce gnuplot script for comparisons
  void prin_compfile(string infile,string infile2)
  {
    initfile((infile+identTag)); //initialise output file
    initfile((infile2+identTag));//initialise output file

    double stuff=0;
    for (int j=0;j<limSub;j++)
      stuff+=inputprop[j];
  
    //make gnuplot file
    fp=fopen("figs/spectrafit.gp","a");
    fprintf(fp,"reset\n");
    fprintf(fp,"set term post eps enh solid color 20\n");
    fprintf(fp,"set output \'figs/test.out.comp.%s.eps\'\n",identTag.c_str());
    fprintf(fp,"set title \'Fitted mass spectrum: %s\'\n",identTag.c_str());
    fprintf(fp,"set xlabel \'m/z\'\n");
    fprintf(fp,"set xrange [%f:%f]\n",minMZ,maxMZ);
    //fprintf(fp,"set yrange [-20:110]\n");
    
    fprintf(fp,"set view map\n");
    fprintf(fp,"set cblabel \'Oligomer\'\n");
    fprintf(fp,"set label \"Mode:      %s \" at graph 0.02,graph 0.95\n",mode.c_str());
    if(mode=="sim")
      {
	fprintf(fp,"set label \"Zfudge:    %.2f \" at graph 0.02,graph 0.9\n",Zfudge);
	fprintf(fp,"set label \"adduction: %.2e \" at graph 0.02,graph 0.85\n",fabs(adduction));
	fprintf(fp,"set label \"MRes:      %.2f \" at graph 0.02,graph 0.8\n",MRes);
	fprintf(fp,"set label \"Resfudge:  %.2e \" at graph 0.02,graph 0.75\n",fabs(ResFudge));
	fprintf(fp,"set label \"Zwidth:    %.2e \" at graph 0.02,graph 0.7\n",fabs(Zwidth));
      }
    else
      {
	fprintf(fp,"set label \"Zfudge:    %.2f (%i)\" at graph 0.02,graph 0.9\n",Zfudge,fitpar.Zfudge_flg);
	fprintf(fp,"set label \"adduction: %.2e (%i)\" at graph 0.02,graph 0.85\n",fabs(adduction),fitpar.adduction_flg);
	fprintf(fp,"set label \"MRes:      %.2f (%i)\" at graph 0.02,graph 0.8\n",MRes,fitpar.MRes_flg);
	fprintf(fp,"set label \"Resfudge:  %.2e (%i)\" at graph 0.02,graph 0.75\n",fabs(ResFudge),fitpar.ResFudge_flg);
	fprintf(fp,"set label \"Zwidth:    %.2e (%i)\" at graph 0.02,graph 0.7\n",fabs(Zwidth),fitpar.Zwidth_flg);
	fprintf(fp,"set label \"shspMass:  %.2e (%i)\" at graph 0.02,graph 0.65\n",fabs(shspMass),fitpar.shspMass_flg);
      }
	
    fprintf(fp,"splot \\\n");
    //  fprintf(fp,"\'out/test.out%s\' u 1:2 ti \'raw\' w li,\\\n",lab.c_str());
    fprintf(fp,"\'out/test.out%s\' u 1:3:(1) ti 'full' w li lt 2,\\\n",identTag.c_str());
    
    //int i=0;
    //for(i=0; i< limSub-1;i++)
    //  fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10):(%i) ti \'%i (%.2f)\' w li %i,\\\n ",lab.c_str(),i,i,i,inputprop[i]/stuff,i+3);
    //fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10):(%i) ti \'%i (%.2f)\' w li %i\n ",lab.c_str(),i,i,i,inputprop[i]/stuff,i+3);
    // fclose(fp);
    
    /*int i=0;
      for(i=0; i< limHSP-1;i++)
      fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3):(%i) noti w li lc palette,\\\n ",lab.c_str(),i,i+1);
      fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3):(%i) noti w li lc palette\n ",lab.c_str(),i,i+1);
      fclose(fp);*/
    
    int i=0;
    if(dim==2)
      {
	for(i=0; i< limHSP-1;i++)
	  fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette,\\\n ",identTag.c_str(),i,i,i+1);
	fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette\n ",identTag.c_str(),i,i,i+1);
      }
    else
      {
	for(i=0; i< limHSP-1;i++)
	  fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette,\\\n ",identTag.c_str(),i,i,i+1);
	fprintf(fp,"\'out/test.out.indiv.%s\' i %i u 1:($3-%i*5):(%i) noti w li lc palette\n ",identTag.c_str(),i,i,i+1);
      }
    fclose(fp);
    
    return;
  }



  

  //# Write input to output file for 3D gnuplot plotting
  /*  void prin_complexes(string outfile)
  {
    fp=fopen( (outfile+identTag).c_str(),"w"); //print input distribution in a square format
    for (int i=0;i<CNT;i++)//print oligomeric state, mass, charge state and weight of ion
      fprintf(fp,"%e\t%e\t%e\t%e\n",Complex_test[i+testmax*0]/shspMass,Complex_test[i+testmax*0],Complex_test[i+testmax*1],Complex_test[i+testmax*2]);
    fclose(fp); 
    }*/
  
  //# Write input to output file for 3D gnuplot plotting
  void prin_input(string outfile)
  {
    fp=fopen((outfile+identTag).c_str(),"w"); //print input distribution in a square format
    for (int i=0;i<limHSP;i++){
      for (int j=0;j<limSub;j++){
	fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0-0.5,input[i+j*limHSP]);
	fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0+0.5,input[i+j*limHSP]);}
      fprintf(fp,"\n");
      for (int j=0;j<limSub;j++){
	fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0-0.5,input[i+j*limHSP]);
	fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0+0.5,input[i+j*limHSP]);}
      fprintf(fp,"\n");
    }
    fclose(fp); 

    fp=fopen("figs/spectrafit.gp","a");
    fprintf(fp,"reset\n");
    fprintf(fp,"set term post eps enh solid color 20\n");
    fprintf(fp,"unset key\n");
    fprintf(fp,"set output \'figs/test.inp%s.eps\'\n",identTag.c_str());
    fprintf(fp,"set title \'Fitted distribution: %s\'\n",identTag.c_str());
    
    /*
    if(chap==0){
      fprintf(fp,"set label \"Gaussian fit\" at graph 0.02,graph 0.95\n");
      fprintf(fp,"set label \"centriod: %.2f\" at graph 0.02,graph 0.9\n",par[0]);
      fprintf(fp,"set label \"sigma:    %.2f\" at graph 0.02,graph 0.85\n",par[1]);
    }
    if(chap==2){
      fprintf(fp,"set label \"AlphaB Poisson fit\" at graph 0.02,graph 0.95\n");
      fprintf(fp,"set label \"MA: %.2f\" at graph 0.02,graph 0.9\n",par[0]);
      fprintf(fp,"set label \"N:  %.2f\" at graph 0.02,graph 0.85\n",par[1]);
    }
    if(chap==7){
      fprintf(fp,"set label \"Single Species fit\" at graph 0.02,graph 0.95\n");
      fprintf(fp,"set label \"Primary Species: %.0f\" at graph 0.02,graph 0.9\n",par[0]);
      for (int i=0;i<par[1];i++)
	fprintf(fp,"set label \"Spec %i: %.2f (%.2f)\" at graph 0.02,graph 0.85\n",i+1,par[2*i+2],par[2*i+3]);
	}*/
    
    
    if(dim==1)
      {
	fprintf(fp,"unset cbtics\n");
	//fprintf(fp,"set xrange[-0.5:%f]\n",II+0.5);
	//fprintf(fp,"set yrange[0:1.05]\n");
	fprintf(fp,"set xlabel \'number of sHSPs\'\n");
	fprintf(fp,"plot \'out/test.inp3d%s\' u 1:2:(1) ti 'distibution' w boxes\n",identTag.c_str());
      }
    else
      {
	fprintf(fp,"set pm3d\n");
	fprintf(fp,"unset key\n");
	fprintf(fp,"set xlabel '[Sub]y'\n");
	fprintf(fp,"set ylabel '[sHSP]x'\n");
	fprintf(fp,"set palette defined (0'white',1'blue')\n");
	fprintf(fp,"set pm3d map\n");
	fprintf(fp,"set size square\n");
	fprintf(fp,"splot 'out/test.inp%s' u 2:1:3\n",identTag.c_str());
      }
    //OutPutFile.close();  
    // os.system('gnuplot '+outfile+'.gp');*/
    
    fprintf(fp,"unset label\n");
    fclose(fp);
    
    return;
  }


  //# Write input to output file for 3D gnuplot plotting
  void prin_input2(string outfile)
  {
    fp=fopen((outfile+identTag).c_str(),"w");
    for (int j=0;j<limSub;j++){
      for (int i=0;i<limHSP;i++){
	fprintf(fp,"%i\t%e\n",i+1,input[i+j*limHSP]);}
      fprintf(fp,"\n\n");
    }
    fclose(fp); 
  }
  
  
  //calculate the sum of total intensity of the best fit spectrum
  double specsum()
  {
    double sum=0;
    for (int j=0;j<lines;j++)
      sum+=data[j+lines*(2)];  
    
    init_array(inputprop,limSub);
    for (int i=0;i<limHSP;i++)
      for (int j=0;j<limSub;j++)
	inputprop[j]+=input[i+j*limHSP];
    
    return sum;
  }


  //extract a specific client composition from the input array
  void input_extract(double *input_comp,int j)
  {
    init_array(input_comp,limHSP*limSub);  	  //make a new input matrix containing only the desired client bound state
    for (int i=0;i<limHSP;i++)
      input_comp[i+j*limHSP]=input[i+j*limHSP];
    return;
  }


  //extract a specific client composition from the input array
  void input_extract_i(double *input_comp,int i,int j)
  {
    init_array(input_comp,limHSP*limSub);  	  //make a new input matrix containing only the desired client bound state
    input_comp[i+j*limHSP]=input[i+j*limHSP];
    return;
  }



  //read distribution in from a text file.
  //store oligomers in freeList and populate input array
  void read_input()
  {
    cout << " Reading input distribution: " << input_file << endl;
    init_array(input,limHSP*limSub);    //first, initialise input array
    
    //string inputfile="raw/equil.txt";
    vector<vector<string> > raw;
    raw = MakeFileVec(input_file); //read file
    for(int ii=0;ii<raw.size();++ii)
      if(raw[ii].size()>0)
	{
	  //cout << "line: " << raw[ii][0] << " " << raw[ii][1] << " " << raw [ii][2] << endl;
	  int i= atoi(raw[ii][0].c_str());
	  int j= atoi(raw[ii][1].c_str());
	  double conc= atof(raw[ii][2].c_str());
	  if(i<limHSP && j<limSub)
	    {
	      //if(i==10 && j==20)
	      input[i+j*limHSP]=conc;
	      //input[i+j*limHSP]=1.;
	    }
	  else
	    {
	      //cout << "Excluding " << i  << " " << limHSP << " " << j << " " << limSub << " " << conc << endl;
	      //exit(100);
	    }
	}
    
    double Sum=0.0;
    for (int i=0;i<limHSP;i++)
      for (int j=0;j<limSub;j++)
	Sum+=input[i+j*limHSP];
    for (int i=0;i<limHSP;i++) //normalise
      for (int j=0;j<limSub;j++)
	input[i+j*limHSP]=input[i+j*limHSP]/Sum;

    /*Sum=0.0;
    for (int i=0;i<limHSP;i++)
      for (int j=0;j<limSub;j++)
	Sum+=input[i+j*limHSP];
	cout << "read sum: " << Sum << endl;*/

    for(int i=0;i<limHSP;++i) //populate freeList
      for(int j=0;j<limSub;++j)
	AddFree(i,j,input[i+j*limHSP]);
  }

  //count the number of distribution parameters we need to fit, if we're fitting.
  void GetDistParams() 
  {
    if(fitpar.dist_flg) //if fitting distribution
      {
	fitpar.distparam=0;
	//for(int i=0;i<distList.size();++i)
	for(std::vector<dist>::iterator it = distList.begin(); it != distList.end(); ++it) 
	  if(it->fit_flg)
	    fitpar.distparam+=it->p;
	for(std::vector<free>::iterator it = freeList.begin(); it != freeList.end(); ++it) 
	  //for(int i=0;i<freeList.size();++i)
	  if(it->fit_flg)
	    fitpar.distparam+=1;
      }
    cout << "Distribution parameter to fit: " << fitpar.distparam << endl;
  }
  
  //print parameter to the screen.
  void ShowSpecPar(string lab,double val,int flg,int &ii)
  {
    if(flg)
      { printf ("%s= %.5f +/- %.5f\n",lab.c_str(),val,err[ii]);ii++;}
    else
      printf ("%s= %.5f\n", lab.c_str(),val);
  }

  //print fitted parameters to the screen.
  void ShowPars()
  {
    int ii=0;
    if(fitpar.dist_flg)
      {
	for(int i=0;i<distList.size();++i)
	  {
	    if(distList[i].fit_flg) //sync error
	      for(int j=0;j<distList[i].p;++j)
		{distList[i].errs[j]=err[ii];++ii;}
	    distList[i].ShowPars();
	  }
	for(int i=0;i<freeList.size();++i)
	  {
	    if(freeList[i].fit_flg) //sync error
	      {freeList[i].err=err[ii];++ii;}
	    if(freeList.size()<10)
	      freeList[i].ShowPars();
	  }
      }
    ShowSpecPar("Adduction",adduction,fitpar.adduction_flg,ii);
    ShowSpecPar("Zfudge   ",Zwidth,fitpar.Zwidth_flg,ii);
    ShowSpecPar("MRes     ",MRes,fitpar.MRes_flg,ii);
    ShowSpecPar("ResFudge ",ResFudge,fitpar.ResFudge_flg,ii);
    ShowSpecPar("Zwidth   ",Zwidth,fitpar.Zwidth_flg,ii);
    ShowSpecPar("shspMass ",shspMass,fitpar.shspMass_flg,ii);

    //for(int i=0;i<fitpar.p;++i)
    //cout << " par " << par[i] << " " << err[i] << endl;
  }

  //transfer par-> spectrum parameters
  void UnPack()
  {
    int ii=0;
    for(std::vector<dist>::iterator it = distList.begin(); it != distList.end(); ++it) 
      //for(int i=0;i<distList.size();++i)
      if(it->fit_flg)
	for(int j=0;j<it->p;++j)
	  {it->pars[j]=fabs(par[ii]);++ii;}
    for(std::vector<free>::iterator it = freeList.begin(); it != freeList.end(); ++it) 
      //for(int i=0;i<freeList.size();++i)
      if(it->fit_flg)
	{it->conc=fabs(par[ii]);++ii;}

    if(fitpar.adduction_flg) //if we're minimising adduction, take the param
      {adduction=fabs(par[ii]);ii++;}
    if(fitpar.Zfudge_flg)    //if we're minimising Zfudge, take the param
      {Zfudge=par[ii];ii++;}
    if(fitpar.MRes_flg)    //if we're minimising MRes, take the param
      {MRes=fabs(par[ii]);ii++;}
    if(fitpar.ResFudge_flg)    //if we're minimising ResFudge, take the param
      {ResFudge=fabs(par[ii]);ii++;}
    if(fitpar.Zwidth_flg)    //if we're minimising residual 12mer, take the param
      {Zwidth=fabs(par[ii]);ii++;}
    if(fitpar.shspMass_flg)    //if we're minimising residual 12mer, take the param
      {shspMass=fabs(par[ii]);ii++;}

    //printf("zz %.10f \n",Zfudge); 
    //cout << "zz " << Zfudge << endl;
    //ShowPars();
  }

  //take values from spectrum and store in par
  void Pack()
  {//begin initiation loop...
    int ii=0;//parameter increment counter
    for(std::vector<dist>::iterator it = distList.begin(); it != distList.end(); ++it) 
      //for(int i=0;i<distList.size();++i)
	if(it->fit_flg) //if we're fitting this distribution
	  for(int j=0;j<it->p;++j)
	    {par[ii]=it->pars[j];++ii;}
    for(std::vector<free>::iterator it = freeList.begin(); it != freeList.end(); ++it) 
      //for(int i=0;i<freeList.size();++i)
      if(it->fit_flg)
	{par[ii]=it->conc;++ii;}
    //place in initial guesses for the mass spec parameters
    if(fitpar.adduction_flg==1)
      {par[ii]=adduction;ii++;}
    if(fitpar.Zfudge_flg==1)
      {par[ii]=Zfudge;ii++;}
    if(fitpar.MRes_flg==1)
      {par[ii]=MRes;ii++;}
    if(fitpar.ResFudge_flg==1)
      {par[ii]=ResFudge;ii++;}
    if(fitpar.Zwidth_flg==1)
      {par[ii]=Zwidth;ii++;}
    if(fitpar.shspMass_flg==1)
      {par[ii]=shspMass;ii++;}
  }//end initiation loop
  
  

    
  //make file to create pretty output
  void make_summary_init(char * infile)
  {
    fp=fopen(infile,"a");
    fprintf(fp,"arraygraph.py 3 4 0 0 0 0 \\\n");
    fclose(fp);
  }
  void make_summary_add(char * infile)
  {
    fp=fopen(infile,"a");
    fprintf(fp,"figs/test.out%s.eps figs/test.inp%s.eps figs/test.out.comp.%s.eps \\\n",identTag.c_str(),identTag.c_str(),identTag.c_str());
    fclose(fp);
  }
  void make_summary_end(char * infile)
  {
    fp=fopen(infile,"a");
    fprintf(fp,"\n");
    fprintf(fp,"mv summary.pdf pdf/%s.pdf\n",ident.c_str());
    fclose(fp);
  }

  //count population in input
  double SumInput()
  {
    double sum=0.0;
    for(int i=0;i<limHSP;++i)
      for(int j=0;j<limSub;++j)
	sum+=input[i+j*limHSP];
    return sum;
  }
  
  //master function to calculate distribution with current state
  void run_calc(int flg)
  {
    if(fitpar.dist_flg!=0) //only do this if fitting a distribution parameter
      make_input(); //Generate input matrix from parameters
    //complex_anal(input);  //analyse input species and calculate ion
    eval_spec(input);     //evaluate spectrum based on Complexes

    //cout << "input sum: " << SumInput() << " ion sum " << CNT << " chi2 " << chi2 << endl;
    if(flg==1) //make outputs if the flag is turned on
      {
	CalcChi2(); //numerically calculate chi2

	prin_input("out/test.inp");    //print input matrix
	prin_input2("out/test.inp3d"); //print input matrix
	prin_spec2("out/test.out");    //print output spectrum
	//prin_complexes("out/test.complex"); //print list of relevant complexes with charge state etc
	double sum=specsum(); //calculate the sum of total intensity of the best fit spectrum
	prin_compfile("out/test.out.comp.","out/test.out.indiv."); //initialise the gnuplot and output files
	if(dim==2)
	  {
	    for (int j=0;j<limSub;j++) //for each complex, calculate the individual contributions for output file
	      {
		input_extract(input_comp,j); //take only input distribution with subunit number j
		//complex_anal(input_comp);    //analyse ions
		eval_spec(input_comp);  //evaluate spectrum
		//eval_spec_norm(spec.datum,spec.lines,sum,inputprop,limSub,j); //optional: normalse?
		prin_spec3("out/test.out.comp.");           //append output spectrum
	      }
	  }
	else if(dim==1)
	  {
	    for (int i=0;i<limHSP;i++) //for each individual complex calculate its mass spectrum
	      {
		//  int i=20; //Hsp number
		//for (int j=0;j<limSub;j++) //for each complex, calculate the individual contributions for output file
		int j=0; //Substrate number
		input_extract_i(input_comp,i,j);
		//complex_anal(input_comp);  //analyse input species
		eval_spec(input_comp);              //evaluate spectrum based on Complexes
		eval_spec_norm_i(sum,input_comp[i+j*limHSP]);
		prin_spec3("out/test.out.indiv.");           //print output spectrum (append each new one to the bottom)
	      }
	  }
	make_summary_add("figs/make_summary.com");
      }
  }
};



