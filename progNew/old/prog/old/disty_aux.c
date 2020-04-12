/*################################################################*/
//# Functions for disty.c

struct mass {
  double thresh;
  double Zwidth; 
  double minMZ;
  double maxMZ;
  double MRes;
  double adduction;
  int minZ;
  int maxZ;
  int limHSP;
  int limSub;
  int testmax;
  double shspMass;
  double Zfudge;
  double ResFudge;
  double tw;
};


//counts in the number of lines in a file (assumes last line includes data)
int countlines(char *filey)
{
  FILE *fp;

  fp=fopen(filey,"r");
  int c=0;char ch='\0';
  while(ch!=EOF) {
    ch=fgetc(fp);
    if(ch=='\n')  c++;
  }
  fclose(fp);
  c++;
  std::cout << "Number of lines in " << filey << ": " << c << std::endl;
  return c;
}

int countfields(char *filey)
{
  ifstream raw1 (filey);
  int i=0;while(raw1.good()){double bg;raw1 >> bg; i++;}
  raw1.close();
  std::cout << "Number of fields in " << filey << ": " << i << std::endl;
  return i;
}


//counts in the number of lines in a file (assumes last line includes data)
int countlines_string(string filey)
{

  FILE *fp;
  fp=fopen(filey.c_str(),"r");
  int c=0;//char ch='\0';
  int ch;
  cout << fp << endl;
  //while(ch!=EOF) {
  {
    ch=fgetc(fp);
    if(ch=='\n')  c++;
  }
  fclose(fp);
  c++;
  std::cout << "Number of lines in " << filey << ": " << c << std::endl;
  return c;
}

int countfields_string(string filey)
{
  ifstream raw1 (filey.c_str());
  int i=0;while(raw1.good()){double bg;raw1 >> bg;i++;}
  raw1.close();
  std::cout << "Number of fields in " << filey << ": " << i << std::endl;
  return i;
}


//initialise array
void init_array(double *array,int N)
{
  for (int i=0;i<N;i++)
    array[i]=0.0;
  return;
}


//initialise a new file
void initfile(char * infile)
{
  FILE *fp;
  fp=fopen(infile,"w");
  fclose(fp); 
  return;
}



//helper function to read in files and output the corresponding array
void readfile(double *data,int lines,string infile)
{
  string raw_in=infile;
  //char *add = ".txt";
  //raw_in.append(add,4);

  {
    //read from the input file
    ifstream raw (infile.c_str());
    int i=0;
    while(raw.good())
      {
	raw >> data[i] >> data[i+lines];
	i++;
      }
    raw.close();
    cout << "i is equal to " << i << endl;
    
    for (int j=0;j<i;j++)
          printf("%f\t%f\n",data[i],data[i+lines]);
  }

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


//normalise highest point on input spectrum to 100
void normdata(double *data,int lines)
{
  
  int imax=findmax(data,lines,1);
  double max=data[imax+lines*1];
  for(int i=0;i<lines;i++)
    data[i+lines*1]=data[i+lines*1]/max*100;
  return;
}   




//# Controls centroid of charge state distributions.
//#averageZ = x* tempM ** y
double avgZ(double tempM,double Zfudge){
  return 0.0467 * pow(tempM,0.533)+Zfudge; //Values from Justin are 0.0467, 0.533;  
}



//# Detector response
//#DetectionEfficient =  x*(1-exp(-y*(tempZ*9.1/tempM)**z))
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


//# Generate 1D list of components from 2D input
void complex_anal(double *input,double *Complexes,int limHSP,int limSub,double adduction,double Zwidth,double shspMass,double subMass,double Zfudge){
  for (int i=0;i<limHSP;i++)
    for (int j=0;j<limSub;j++)
      {
	//double mass=1.0*adduction*(((i+1)*shspMass*1.0)+((j)*subMass*1.0));
	double mass=1.0*(((i+1)*shspMass*1.0)+((j)*subMass*1.0));
	mass=mass+fabs(adduction)*pow(mass,0.760);
	//mass=mass*(adduction+1.0);
	Complexes[(i+j*limHSP)+limHSP*limSub*0]=mass;             //mass
	Complexes[(i+j*limHSP)+limHSP*limSub*1]=input[i+limHSP*j];//abundance
	Complexes[(i+j*limHSP)+limHSP*limSub*2]=avgZ(mass,Zfudge);       //average charge state
	Complexes[(i+j*limHSP)+limHSP*limSub*3]=Zwidth;           //width of the charge state
      }
  return;
}



//# Finds all of the charge states for each component that might contribute to the mass spectrum.   
int complex_def(double *Complexes,double *Complex_test,int limHSP,int limSub,int testmax,int minZ,int maxZ,double thresh)
{
  int CNT=0;
  for (int i=0;i<limHSP*limSub;i++)//for each entry in Complexes array
    { 
      double TotalWeight = 0;              // Intensity for all charge states.
      double candidateZ[(maxZ+1-minZ)*2];  // Array for all possible charge states.
      
      // Assigns an intensity for all candidate charge states, based on a Gaussian distribution in z-space
      for(int j=0;j<maxZ+1-minZ;j++)       // loop over all charge states:
	{           
	  candidateZ[j+(maxZ+1-minZ)*0]=minZ+j; //current charge state
	  //weighting - x=charge state,centroid=adjusted mass,resolution?
	  double CurrentWeight = Norm(candidateZ[j+(maxZ+1-minZ)*0],Complexes[i+limHSP*limSub*2],Complexes[i+limHSP*limSub*3]);
	  candidateZ[j+(maxZ+1-minZ)*1]=CurrentWeight;
	  TotalWeight += CurrentWeight;   	//Updates the total weight for all z of this component.

	}
      
      //Adds entry for all charge states above the defined threshold to element for the complex.
      for(int j=0;j<maxZ+1-minZ;j++){  //  loop over all charge states
        // Normalizes intensities for charge state peaks based on the abundance from the input files
	double CorrectWeight = Complexes[i+limHSP*limSub*1]*candidateZ[j+(maxZ+1-minZ)*1]/TotalWeight;
	if (CorrectWeight >= thresh)
	  {
	    Complex_test[CNT+testmax*0]=Complexes[i];   //mass of candidate complex
	    Complex_test[CNT+testmax*1]=candidateZ[j];  //charge state of candidate complex
	    Complex_test[CNT+testmax*2]=CorrectWeight*DetectionEfficiency(Complexes[i],candidateZ[j]); //weight of candidate complex
	    //	    cout << CNT << " " << Complexes[i] << " " << candidateZ[j] << endl;
	    CNT++;
	  }
      }
    }
  //cout << "Number of complexes above detection threshold: " << CNT << endl;
  if(CNT > testmax) cout << "PROBLEM - need to increase size of complex trial array" << endl;

  //  for (int i=0; i<CNT;i++)
  // {cout << CNT << " " << Complex_test[i] << " " << Complex_test[i+testmax*1]<< " " << Complex_test[i+testmax*2]<<endl;}

  return CNT;
}
      
/*
//# Create array for spectrum, consisting of [m/z, intensity] elements.    
double init_spec(double *Spectrum,double MZmode,double minMZ,double maxMZ,double MZstep,double MZsample,double *datainp)
{
  double currentMZ = minMZ;
  if (MZmode == 2)
    while (currentMZ <= maxMZ)
      {
	Spectrum.append([currentMZ,0]);
	currentMZ += MZstep;
      }
  if (MZmode == 0)
    while (currentMZ <= maxMZ)
      {
	Spectrum.append([currentMZ,0]);
	currentMZ += currentMZ / (MRes * MZsample);
      }
  if (MZmode==3)
    for (i in range(len(datainp)))
      Spectrum.append([float(datainp[i][0]),0]);

  return Spectrum;
}
*/

void reset_sim(double *data,int lines)
{
  for (int i=0;i<lines;i++)
    data[i+lines*2]=0;
  return;
}




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


void eval_spec_norm_i(double *data,int lines,double sum,double *inputprop,double poppy,int limSub,int i)
{
  double sup=0;
  for (int j=0;j<lines;j++)
    sup+=data[j+lines*(2)];  

  double stuff=0;
  for (int j=0;j<limSub;j++)
    stuff+=inputprop[j];
  printf("    Guy with %i clients bound: %.2f \n",i,poppy/stuff);

  for (int j=0;j<lines;j++)  
    data[j+lines*(2)]=    data[j+lines*(2)]*sum/sup*(poppy/stuff);  
  return;
}



//# Evaluate spectrum for each trial component
double eval_spec(double *Complex_test,int cnt,int testmax,double *data,int lines,double MRes,double ResFudge)
{
  for (int i=0;i<cnt;i++)//for all trial complexes
    {
      // Define bounds for current charge state
      double centroid = (Complex_test[i+testmax*0] )/Complex_test[i+testmax*1]; //mass over charge
      double lowMZ = centroid - (5 * centroid*((1/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));     //lowest MZ is 
      double highMZ = centroid + (5 *centroid*((1/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));
    
      // Evaluate contribution of current charge state
      for (int j=0;j<lines;j++)//for each point in spectrum
	if(data[j]>=lowMZ && data[j]<=highMZ) //if it is within range...
	    //# Evaluate Gaussian function for (current m/z, centroid m/z, MS resolution)
	  data[j+lines*2]+= Complex_test[i+testmax*2]*Norm(data[j],centroid,centroid*((1/MRes*1.0)+1.0*((fabs(ResFudge)/100000000))*centroid));
    }

  double rawmax=data[findmax(data,lines,1)+1*lines];
  double simmax=data[findmax(data,lines,2)+2*lines];
  //cout << "datamax " << rawmax << " simmax " << simmax << endl;
  for (int j=0;j<lines;j++)
    data[j+lines*2]=data[j+lines*2]/simmax*rawmax;
  double chi2=0;
  for (int j=0;j<lines;j++)
    chi2+=pow((data[j+lines*(1)]-data[j+lines*(2)]),2.0)/100;
  return chi2;
}







//# Write spectrum to output file
void prin_spec(string infile,char* lab,int tag,double *array,int N,int col)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab,tag);
  
  fp=fopen(outfile.c_str(),"w");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\n",array[i],array[i+N*col]);
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
void prin_spec2(string infile,string lab,int tag,double *array,int N,double chi2)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  
  fp=fopen(outfile.c_str(),"w");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\t%f\n",array[i],array[i+N*1],array[i+N*2]);
  fclose(fp); 
  
  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"set output \'figs/test.out%s.eps\'\n",lab.c_str());
  fprintf(fp,"set title \'%s chi2/dof: %.2f ave error: %.2f%\'\n",lab.c_str(),chi2,sqrt(chi2*100));
  fprintf(fp,"set xlabel \'m/z\'\n");
  fprintf(fp,"set xrange [6000:14000]\n");
  fprintf(fp,"set yrange [-20:110]\n");
  fprintf(fp,"plot \'out/test.out%s\' u 1:2 ti 'raw' w li,\\\n",lab.c_str());
  fprintf(fp,"\'\' u 1:3 ti 'fitted' w li,\\\n",lab.c_str());
  fprintf(fp,"\'\' u 1:(($2-$3)-10) ti 'difference' w li\n",lab.c_str());
  fclose(fp);
  return;
}



//initialise comparison file and produce gnuplot script for comparisons
void prin_compfile(string infile,string infile2,string infile3,string lab,int tag,int limSub,double *inputprop,double Zfudge,double adduction,double MRes,double ResFudge)
{

  //initialise file
  FILE *fp;

  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  fp=fopen(outfile.c_str(),"w");  
  fclose(fp);


  string outfile2=infile2;
  //char *add = ".out";
  outfile2.append(lab.c_str(),tag);
  fp=fopen(outfile2.c_str(),"w");  
  fclose(fp);

  string outfile3=infile3;
  //char *add = ".out";
  outfile3.append(lab.c_str(),tag);
  fp=fopen(outfile3.c_str(),"w");  
  fclose(fp);



  double stuff=0;
  for (int j=0;j<limSub;j++)
    stuff+=inputprop[j];
  
  //make gnuplot file
  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"set output \'figs/test.out.comp.%s.eps\'\n",lab.c_str());
  fprintf(fp,"set title \'%s\'\n",lab.c_str());
  fprintf(fp,"set xlabel \'m/z\'\n");
  fprintf(fp,"set xrange [6000:14000]\n");
  fprintf(fp,"set yrange [-20:110]\n");
  fprintf(fp,"set label \"Zfudge:    %.2f\" at graph 0.02,graph 0.95\n",Zfudge);
  fprintf(fp,"set label \"adduction: %.2e\" at graph 0.02,graph 0.9\n",adduction);
  fprintf(fp,"set label \"MRes:      %.2f\" at graph 0.02,graph 0.85\n",MRes);
  fprintf(fp,"set label \"Resfudge:  %.2e\" at graph 0.02,graph 0.8\n",ResFudge);
  fprintf(fp,"plot \\\n");
  //  fprintf(fp,"\'out/test.out%s\' u 1:2 ti \'raw\' w li,\\\n",lab.c_str());
  fprintf(fp,"\'out/test.out%s\' u 1:3 ti 'fitted' w li 2,\\\n",lab.c_str());
  int i=0;
  for(i=0; i< limSub-1;i++)
    fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10) ti \'%i (%.2f)\' w li %i,\\\n ",lab.c_str(),i,i,inputprop[i]/stuff,i+3);
  fprintf(fp,"\'out/test.out.comp.%s\' i %i u 1:($3-10) ti \'%i (%.2f)\' w li %i\n ",lab.c_str(),i,i,inputprop[i]/stuff,i+3);
  fclose(fp);

  return;
}



//# Write spectrum to output file
void prin_spec3(string infile,string lab,int tag,double *array,int N,int col)
{        
  FILE *fp;
  string outfile=infile;
  //char *add = ".out";
  outfile.append(lab.c_str(),tag);
  
  fp=fopen(outfile.c_str(),"a");
  for (int i=0;i<N;i++)
    fprintf(fp,"%f\t%f\t%f\n",array[i],array[i+N*1],array[i+N*2]);
  fprintf(fp,"\n\n");
  fclose(fp); 
  
  return;
}






//# Write input to output file for 3D gnuplot plotting
void prin_input(string outfile,string label,int tag,double *input,int II,int JJ)
{
  FILE *fp;


  outfile.append(label.c_str(),tag);

  //  fp=fopen(outfile.c_str(),"w");
  //for (int i=0;i<II;i++){
  //  for (int j=0;j<JJ;j++){
  //    fprintf(fp,"%i\t%i\t%e\n",i+1,j,input[i+j*II]);}
  //  fprintf(fp,"\n");
  // }
  //fclose(fp); 


  fp=fopen(outfile.c_str(),"w"); //print input distribution in a square format
  for (int i=0;i<II;i++){
    for (int j=0;j<JJ;j++){
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0-0.5,input[i+j*II]);
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1-0.5,j*1.0+0.5,input[i+j*II]);}
    fprintf(fp,"\n");
    for (int j=0;j<JJ;j++){
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0-0.5,input[i+j*II]);
      fprintf(fp,"%e\t%e\t%e\n",i*1.0+1+0.5,j*1.0+0.5,input[i+j*II]);}
    fprintf(fp,"\n");
  }
  //fclose(fp); 


  fp=fopen("figs/spectrafit.gp","a");
  fprintf(fp,"reset\n");
  fprintf(fp,"set term post eps enh solid color 20\n");
  fprintf(fp,"set pm3d map\n");
  fprintf(fp,"unset key\n");
  fprintf(fp,"set palette defined (0\'white\',0.4\'white\',1\'blue\')\n");
  fprintf(fp,"set output \'figs/test.inp%s.eps\'\n",label.c_str());
  fprintf(fp,"set title \'%s\'\n",label.c_str());
  fprintf(fp,"unset cbtics\n");
  fprintf(fp,"set yrange[-0.5:%f]\n",JJ+0.5);
  fprintf(fp,"set xlabel \'number of sHSPs\'\n");
  fprintf(fp,"set ylabel \'number of clients\'\n");
  fprintf(fp,"splot \'out/test.inp%s\' u 1:2:3 ti 'raw'\n",label.c_str());
  fclose(fp);


  /*OutPutFile=open(outfile+'.gp','w');
  fprintf('set term post eps enh solid color 20\n');
  fprintf('set output \''+outfile+'.eps\'\n');
  fprintf('set title \''+outfile+' \'\n');
  fprintf('set pm3d\n');
  fprintf('unset key\n');
  fprintf('set xlabel \'[Sub]y\'\n');
  fprintf('set ylabel \'[sHSP]x\'\n');
  fprintf('set palette defined (0\'white\',1\'blue\')\n');
  fprintf('set pm3d map\n');
  fprintf('set size square\n');
  fprintf('splot \''+outfile+'\' u 2:1:3');
  OutPutFile.close();  
  os.system('gnuplot '+outfile+'.gp');*/
  return;
}


//# Write input to output file for 3D gnuplot plotting
void prin_input2(string outfile,string label,int tag,double *input,int II,int JJ)
{
  FILE *fp;
  outfile.append(label.c_str(),tag);
  fp=fopen(outfile.c_str(),"w");
  for (int j=0;j<JJ;j++){
    for (int i=0;i<II;i++){
      fprintf(fp,"%i\t%e\n",i+1,input[i+j*II]);}
    fprintf(fp,"\n\n");
  }
  fclose(fp); 

  return;
}



//calculate the sum of total intensity of the best fit spectrum
double specsum(double *data,int lines,double *input,double *inputprop,int limHSP,int limSub)
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
void input_extract(double *input_comp,double *input,int SubSpec,int limHSP)
{
  for (int i=0;i<limHSP;i++)
    input_comp[i+SubSpec*limHSP]=input[i+SubSpec*limHSP];
  return;
}


//extract a specific client composition from the input array
void input_extract_i(double *input_comp,double *input,int limHSP,int HspSpec,int SubSpec)
{
  input_comp[HspSpec+SubSpec*limHSP]=input[HspSpec+SubSpec*limHSP];
  return;
}



void addtw(double *input,int limHSP,int limSub,double tw)
{
  //find sum over all input array
  double sum=0;
  for (int i=0; i<limHSP;i++)
    for (int j=0; j<limSub;j++)
      sum+=input[i+j*limHSP];
  input[11]=sum*fabs(tw);
  return;
}


		   

//# let distribution of each Hsp x0,sigma, Sub x0,sigma
void make_input(double *input,int limHSP,int limSub,int chap,double* par,double tw)
{
  for (int i=0;i<limHSP;i++)
    {
      for (int j=0;j<limSub;j++)
	{
	  
	  if(chap==0 || chap==8)//2d Gaussian (chap8 is if we don't want to minimise par0 and 2 in grid search)
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      double Sx0=par[2];
	      double Ssig=par[3];
	      input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig))*(Norm((i+1)*1.0,Hx0,Hsig));
	    }
	  if(chap==1 || chap==9)//2d Skewed Gaussian (chap9 is if we don't want to minimise par0 and 2 in grid search)
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      double Sx0=par[2];
	      double Ssig=par[3];
	      double alpha=par[4];
	      input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig))*(SkewNorm((i+1)*1.0,Hx0,Hsig,alpha));
	    }
	  if(chap==2)//2d Skewed twisted Gaussian
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      double Sx0=par[2];
	      double Ssig=par[3];
	      double alpha=par[4];
	      double skew=par[5];
	      input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig) * SkewNorm((i+1)*1.0,Hx0,Hsig,alpha) * exp(skew*(j*1.0-Sx0)*(i*1.0+1.0-Hx0)/(Ssig*Hsig)) );
	    }
	  if(chap==3)//Skewed in both dimensions
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      double Sx0=par[2];
	      double Ssig=par[3];
	      double alpha=par[4];
	      double alpha2=par[5];
	      input[i+j*limHSP]=(SkewNorm(j*1.0,Sx0,Ssig,alpha2) * SkewNorm((i+1)*1.0,Hx0,Hsig,alpha) );
	    }
	  if(chap==4)//have a 2d skewed Gaussian but have client=0 free
	    {
	      if(j==0)
		{
		  double Hx01=par[5];
		  double Hsig1=par[6];
		  double alpha1=par[7];
		  double rel = par[8];
		  input[i+j*limHSP]= fabs(rel) * (SkewNorm((i+1)*1.0,Hx01,Hsig1,alpha1)  );
		}
	      else
		{
		  double Hx0=par[0];
		  double Hsig=par[1];
		  double Sx0=par[2];
		  double Ssig=par[3];
		  double alpha=par[4];
		  input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig) * SkewNorm((i+1)*1.0,Hx0,Hsig,alpha) );
		}
	    }
	  
	  
	  
	  if(chap==5)//Have each client number slice have its own gaussian
	    {
	      double Hx0 =par[0+3*j];
	      double Hsig=par[1+3*j];
	      double rel =par[2+3*j];
	      if(rel>1E-4)
		input[i+j*limHSP]=rel*(Norm((i+1)*1.0,Hx0,Hsig));
	      else
		input[i+j*limHSP]=0;
	    }
	  
	  
	  
	  if(chap==6)//Have each client number slice have its skewed Gaussian
	    {
	      double Hx0  =par[0+4*j];
	      double Hsig =par[1+4*j];
	      double alpha=par[2+4*j];
	      double rel  =par[3+4*j];

	      if(rel>1E-4)
		input[i+j*limHSP]=rel*(SkewNorm((i+1)*1.0,Hx0,Hsig,alpha));
	      else
		input[i+j*limHSP]=0.0;
	      
	      
	    }
	  

	  if(chap==7)//2d Skewed twisted Gaussian
	    {
	      double Hx0=par[0];
	      double Hsig=par[1];
	      double Sx0=par[2];
	      double Ssig=par[3];
	      double alpha=par[4];
	      input[i+j*limHSP]=(Norm(j*1.0,Sx0,Ssig) * Norm((i+1)*1.0,Hx0,Hsig) * exp((j*1.0-Sx0)*(i*1.0+1-Hx0)/(Ssig*Hsig)) );
	    }

	  if( chap>=10 )
	    {
	      input[i+j*limHSP]=par[i+j*limHSP];
	    }
	}
    }

  if(tw>0) addtw(input,limHSP,limSub,tw); //if tw>0, add 12mer to the count


  return;
}



void make_summary_init(char * infile)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"arraygraph.py 3 4 0 0 0 0 \\\n");
  fclose(fp);

  return;
}


void make_summary_add(char * infile,string lab)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"figs/test.out%s.eps figs/test.inp%s.eps figs/test.out.comp.%s.eps \\\n",lab.c_str(),lab.c_str(),lab.c_str());
  fclose(fp);

  return;
}


void make_summary_end(char * infile,string lab)
{
  FILE *fp;

  fp=fopen(infile,"a");
  fprintf(fp,"\n");
  fprintf(fp,"mv summary.pdf pdf/%s.pdf\n",lab.c_str());
  fclose(fp);

  return;
}





  
double run_calc(int flg,int chap,string lab,int tag,int limHSP,int limSub,double *par,double adduction,double Zwidth,double shspMass,double subMass,int testmax,int minZ,int maxZ,double thresh,double* data,int lines,double MRes,double Zfudge,double ResFudge,double tw)
{
  double input[limHSP*limSub];init_array(input,limHSP*limSub);
  double Complexes[limHSP*limSub*5];init_array(Complexes,limHSP*limSub*5);
  double Complex_test[testmax*3];init_array(Complex_test,testmax*3);
  


  reset_sim(data,lines);   //reset calculated spectrum
  make_input(input,limHSP,limSub,chap,par,tw);                                           //Generate input matrix
  complex_anal(input,Complexes,limHSP,limSub,adduction,Zwidth,shspMass,subMass,Zfudge);  //analyse input species
  int cnt=complex_def(Complexes,Complex_test,limHSP,limSub,testmax,minZ,maxZ,thresh);    //find all charge states for input distribution
  double chi2=eval_spec(Complex_test,cnt,testmax,data,lines,MRes,ResFudge);              //evaluate spectrum based on Complexes

  if(flg==1) //make outputs if the flag is turned on
  {
	  prin_input("out/test.inp",lab,tag,input,limHSP,limSub);    //print input matrix
	  prin_input2("out/test.inp3d",lab,tag,input,limHSP,limSub); //print input matrix
	  prin_spec2("out/test.out",lab,tag,data,lines,chi2/lines);  //print output spectrum


	  double inputprop[limSub];double sum=specsum(data,lines,input,inputprop,limHSP,limSub); //calculate the sum of total intensity of the best fit spectrum
	  prin_compfile("out/test.out.comp.","out/test.out.indiv1.","out/test.out.indiv2.",lab,tag,limSub,inputprop,Zfudge,adduction,MRes,ResFudge); //initialise the gnuplot and output files

	  for (int i=0;i<limSub;i++) //for each complex, calculate the individual contributions for output file
	    {
	      init_array(Complexes,limHSP*limSub*5);
	      init_array(Complex_test,testmax*3);
	      reset_sim(data,lines);   //reset calculated spectrum
	      double input_comp[limHSP*limSub];init_array(input_comp,limHSP*limSub);  	  //make a new input matrix containing only the desired client bound state
	      input_extract(input_comp,input,limHSP,i);
	      complex_anal(input_comp,Complexes,limHSP,limSub,adduction,Zwidth,shspMass,subMass,Zfudge);  //analyse input species
	      int cnt_comp=complex_def(Complexes,Complex_test,limHSP,limSub,testmax,minZ,maxZ,thresh);    //find all charge states for input distribution
	      double chi2_comp=eval_spec(Complex_test,cnt_comp,testmax,data,lines,MRes,ResFudge);              //evaluate spectrum based on Complexes
	      eval_spec_norm(data,lines,sum,inputprop,limSub,i);
	      prin_spec3("out/test.out.comp.",lab,tag,data,lines,3);           //print output spectrum
	    }
	  

	  for (int i=0;i<limHSP;i++) //for each individual complex calculate its mass spectrum
	  {

	    //  int i=20; //Hsp number
	    //	    for (int j=0;j<limSub;j++) //for each complex, calculate the individual contributions for output file
	    for (int ij=0;ij<2;ij++)

	      {
		int j=0;
		if(ij==0)
		  j=1; //Substrate number
		if(ij==1)
		  j=2;//Substrate number

		init_array(Complexes,limHSP*limSub*5);
		init_array(Complex_test,testmax*3);
		reset_sim(data,lines);   //reset calculated spectrum
		double input_comp[limHSP*limSub];init_array(input_comp,limHSP*limSub);  	  //make a new input matrix containing only the desired client bound state
		input_extract_i(input_comp,input,limHSP,i,j);

		if(j==1)
		  prin_input("out/test.inp.indiv1.",lab,tag,input_comp,limHSP,limSub);    //print input matrix
		if(j==2)
		  prin_input("out/test.inp.indiv2.",lab,tag,input_comp,limHSP,limSub);    //print input matrix


		complex_anal(input_comp,Complexes,limHSP,limSub,adduction,Zwidth,shspMass,subMass,Zfudge);  //analyse input species
		int cnt_comp=complex_def(Complexes,Complex_test,limHSP,limSub,testmax,minZ,maxZ,thresh);    //find all charge states for input distribution
		double chi2_comp=eval_spec(Complex_test,cnt_comp,testmax,data,lines,MRes,ResFudge);              //evaluate spectrum based on Complexes

		//  input_comp[HspSpec+SubSpec*limHSP]=input[HspSpec+SubSpec*limHSP];
		eval_spec_norm_i(data,lines,sum,inputprop,input_comp[i+j*limHSP],limSub,j);

		if(j==1)
		  prin_spec3("out/test.out.indiv1.",lab,tag,data,lines,3);           //print output spectrum (append each new one to the bottom)
		if(j==2)
		  prin_spec3("out/test.out.indiv2.",lab,tag,data,lines,3);           //print output spectrum (append each new one to the bottom)
	      }
	  }

	  make_summary_add("figs/make_summary.com",lab);
  }
  return chi2/lines;
}


