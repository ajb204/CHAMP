
#include <gsl/gsl_rng.h>      //used by jiggler
#include <gsl/gsl_randist.h>  //used by jiggler
#include <gsl/gsl_multifit_nlinear.h> //used by fitty
#include <gsl/gsl_blas.h>     //used by fitty

//move gsl vector intoa spec array
void UnPack(const gsl_vector *x,mass *spec)
{
   for(int i=0;i<(*spec).fitpar.p;++i)
     {
        (*spec).par[i]=gsl_vector_get(x,i);
     }
   (*spec).UnPack();

}

// CALCULATE CHI^2 MATRIX
int expb_f (const gsl_vector * x, void *data,gsl_vector * f)
{
  mass *spec=(mass* )data; //recast void pointer
  UnPack(x,spec);     //move parameters into class
  (*spec).run_calc(0); //calculate model with current settings
  for (size_t i = 0; i < (*spec).lines; i++)   //record and store chi values
    //chi(i) = (Ydata(i)-Ysim(i))/sigma
    gsl_vector_set (f, i,  ((*spec).data[i+(*spec).lines*1]-(*spec).data[i+(*spec).lines*2]) );
  return GSL_SUCCESS;
}


void callback(const size_t iter, void * data,const gsl_multifit_nlinear_workspace *w)
{
  mass *spec=(mass* )data; //recast void pointer
  gsl_vector * x = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);

  printf ("iter: %3zu x = ",iter);
//if(p<10)
//    for(int i=0;i<p;i++)
//     printf("% 15.8f ",gsl_vector_get(s->x,i));
//  printf(" |f(x)| = %g\n", gsl_blas_dnrm2 (s->f));
  /* print out current location */
  //cout << " " << (*spec).fitpar.distparam << endl; // +(*spec).fitpar.specparam << endl;
  //if((*spec).fitpar.distparam+(*spec).fitpar.specparam<10)
  //for(int i=0;i<(*spec).fitpar.distparam+(*spec).fitpar.specparam;i++)
  //    printf("% 15.8f ",gsl_vector_get(x,i));
  //printf("% 15.8f ",gsl_vector_get(x,0));
   printf(" |f(x)| = %g\n", gsl_blas_dnrm2 (f));
  //        printf("%f %f\n",
//        gsl_vector_get(x, 0),
//        gsl_vector_get(x, 1));
}


//MAIN FUNCTION. SETS UP AND CALLS MINIMSER
void fitty(int flg,mass &spec,double maxIter)
{
  spec.GetDistParams(); //set number of distribution parameters
  spec.fitpar.GetSpecParams(); //set number of spec pars, and total.
  spec.Pack(); //fill up parameters into initial 'par' array

  if(flg)
    spec.ShowPars();

  size_t p = spec.fitpar.p; //total parameters in model.

  //if(spec.fitpar.distparam>0 && spec.chap>=10 ) //for completely free distribution parameters
  //  p=p-1;  //subtract 1

  const size_t n = spec.lines; //number of datapoints

  gsl_matrix *covar = gsl_matrix_alloc (p, p);   //initialise minimiser
  gsl_matrix *J;

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
  gsl_multifit_nlinear_workspace *w  = gsl_multifit_nlinear_alloc(T,&params,n,p);
  gsl_multifit_nlinear_fdf f;

  const double xtol = 1e-8; //convergence condition1
  const double gtol = 1e-8; //convergence condition2
  const double ftol = 0;    //convergence condition3

  int status,info; //integers to follow progress

  gsl_vector_view x = gsl_vector_view_array (spec.par, p);  //link initial parameters
  f.f = &expb_f;
  f.df= NULL; //do this numerically
  f.n = n;
  f.p = p;
  f.params = &spec;
  if(flg==1)
    cout << "Fitting with "<< p << " total parameters " << endl;

  gsl_multifit_nlinear_init(&x.vector,&f,w);

  if(flg) //callbacks
    status=gsl_multifit_nlinear_driver(maxIter,xtol,gtol,ftol,callback,NULL,&info,w);
  else //no callbacks
    status=gsl_multifit_nlinear_driver(maxIter,xtol,gtol,ftol,NULL,NULL,&info,w);

  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar ( J,0.0, covar);

  const gsl_vector * x0 = gsl_multifit_nlinear_position(w);
  gsl_vector *f0 = gsl_multifit_nlinear_residual(w);
  UnPack(x0,&spec); //load the fitting parameters back into the model.
  spec.run_calc(flg); //run the model to get final state
  gsl_vector_view err=gsl_matrix_diagonal(covar); //unpack covarience matrix for errors
  if(flg==1){
    printf ("internal chi2/dof %f\n",spec.chi2/(spec.lines-spec.fitpar.p));
    printf ("status = %s\n", gsl_strerror (status));}

  double chi = gsl_blas_dnrm2(f0);  //chi
  double dof = n - p;               //dof
  double c = GSL_MAX_DBL(1, chi / sqrt(dof));  //reduced chi^2

  if(flg==1)
    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

  for(int i=0;i<spec.fitpar.p;++i) //get errors
    { //save parameter and errors
      spec.err[i]=c*sqrt(gsl_matrix_get(covar,i,i));
      spec.par[i]=gsl_vector_get(x0,i);
    }
  if(flg)
    spec.ShowPars();

  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);

}

/*************************/
/* NOTE gridsearch functions need updating */

    /*
void gridsearch1D(int i,int j,int trial,int lines,double *data,double *x_init,struct mass& mass_init,struct control& params,string *raw,string *ident,double Max0,double Min0,double Max1,double Min1,char *GridType)
  {

    int maxsearch=(Max0-Min0)/2.0+1; //stepsize  for centroid


    //initiate a grid search output file and mode
    FILE *fp;
    string labby="out/chi2init.out";
    labby.append(ident[i].c_str(),4);

    int jj=0;
    if(j==0){//for Gaussian
      jj=0;
      char *lab="_0";labby.append(lab,2);}
    if(j==1){//for Skewed Gaussian
      jj=1;
      char *lab="_1";labby.append(lab,2);}
    fp=fopen(labby.c_str(),"w");
    fclose(fp);

    if(params.spec.chap==0)
      params.spec.chap=8;
    if(params.spec.chap==2)
      params.spec.chap=3;
    //run the gridsearch to find optimum centre of substrate and client distributions
    double testym[maxsearch*maxsearch*3];
    int cnt=0;
    for (int k1=0;k1<maxsearch;k1++)
      {
	cout << "    Gridsearch progress: " << k1+1 << " of " << maxsearch << endl;
	for (int k2=0;k2<maxsearch;k2++)
	  {
	    testym[cnt+maxsearch*maxsearch*0]=Min0+(k1)/(maxsearch*1.0-1.0)*(Max0-Min0);  //grid search from 10 to 50 on 'rat'
	    if(GridType=="log")
	      testym[cnt+maxsearch*maxsearch*1]=Min1*pow(10,k2/(maxsearch*1.0-1.0)*log10(Max1/Min1));
	    else
	      testym[cnt+maxsearch*maxsearch*1]=Min1+(k2)/(maxsearch*1.0-1.0)*(Max1-Min1);  //grid search from 10 to 50 on 'rat'


	    //x_init[0] = 25.0;    //Hx0    substrate centre
	    //x_init[1] = 2.0;     //Hsig   substrate width
	    //x_init[2] = 2.0;     //Sxo    client centre
	    //x_init[3] = 1.0;     //Ssig   client width
	    //x_init[4] = 0.01;    //alpha  skew factor

	    x_init[0]=testym[cnt+maxsearch*maxsearch*0];
	    x_init[1]=testym[cnt+maxsearch*maxsearch*1];
	    x_init[2]=2.0;
	    x_init[3]=0.5;
	    x_init[4]=0.01;

	    struct mass massy;
	    massy=mass_init;
	    double chi2=runchap("fitty",trial,lines,data,"Null",params,x_init,massy,0);  // skewed gaussian with client=0 seperate + all params
	    //	      if(x_init[2]>4.0)//push the chi^2 up if finding a solution with too many clients
	    //	chi2=chi2+10;
	    //if(x_init[0]<12.0)//push the chi^2 up if finding a solution with too few HSP
	    //	chi2=chi2+10;
	    if(chi2!=chi2)
	      chi2=1E6;
	    testym[cnt+maxsearch*maxsearch*2]=chi2;

	    fp=fopen(labby.c_str(),"a");
	    fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",testym[cnt+maxsearch*maxsearch*0],testym[cnt+maxsearch*maxsearch*1],chi2,x_init[0],x_init[1],x_init[2],x_init[3]);
	    fclose(fp);
	    cnt++;

	  }
	fp=fopen(labby.c_str(),"a");
	fprintf(fp,"\n");
	fclose(fp);
      }

    int imin=findmin(testym,maxsearch*maxsearch,2);//find the lowest chi2 element in the test array


    //Initialise the final minimisation
    x_init[0]=testym[imin+maxsearch*maxsearch*0];   //update initialise matrix with minimsed value
    x_init[1]=testym[imin+maxsearch*maxsearch*1];   //update initialise matrix with minimised value
    x_init[2]=2.0;
    x_init[3]=0.5;
    x_init[4]=0.01;


    //run the final fit with these parameters
    struct mass massy;
    massy=mass_init;
    double chi2=runchap("fitty",trial,lines,data,"Null",params,x_init,massy,0);  // skewed gaussian with client=0 seperate + all params


    //update the struct with the new fitted values ready for main minimisation. Minimiser should start very close to the bottom of the well
    mass_init.MRes      =massy.MRes;
    mass_init.adduction =massy.adduction;
    mass_init.Zfudge    =massy.Zfudge;
    mass_init.ResFudge  =massy.ResFudge;

    if(params.chap==8)
      params.chap=0;
    if(params.chap==3)
      params.chap=2;

    return;
  }




  //Gridsearch over oligomers and zfudge
  void gridsearchSingle(int i,int j,int trial,int lines,double *data,double *x_init,struct mass& mass_init,struct control& params,string *raw,string *ident,double Max0,double Min0)
  {

    int maxsearch=(Max0-Min0)+1; //stepsize  for centroid

    //initiate a grid search output file and mode
    FILE *fp;
    string labby="out/chi2init.out";
    labby.append(ident[i].c_str(),4);

    int jj=0;
    if(j==0){//for Gaussian
      jj=0;
      char *lab="_0";labby.append(lab,2);}
    if(j==1){//for Skewed Gaussian
      jj=1;
      char *lab="_1";labby.append(lab,2);}
    fp=fopen(labby.c_str(),"w");
    fclose(fp);
    cout << "    Gridsearching... " << endl;
    //run the gridsearch to find optimum centre of substrate and client distributions
    int Zsearch=5;
    int AddSearch=2;
    double testym[Zsearch*maxsearch*AddSearch*4];
    int cnt=0;
    for (int k1=0;k1<maxsearch;k1++)
      {
	for  (int k2=0;k2<Zsearch;k2++)
	  {
	    for  (int k3=0;k3<AddSearch;k3++)
	      {

		testym[cnt+maxsearch*Zsearch*AddSearch*0]=Min0+(k1)/(maxsearch*1.0-1.0)*(Max0-Min0);  //grid search from 10 to 50 on 'rat'
		testym[cnt+maxsearch*Zsearch*AddSearch*1]=-2.+k2*1.;  //zfudge value
		if(k3==0)
		  testym[cnt+maxsearch*Zsearch*AddSearch*2]=1E-3;  //adduction value
		if(k3==1)
		  testym[cnt+maxsearch*Zsearch*AddSearch*2]=0.1;  //adduction value


		//	  testym[cnt+maxsearch*Zsearch*0]=24.0;  //grid search from 10 to 50 on 'rat'
		struct mass massy;
		massy=mass_init;
		x_init[0]      =testym[cnt+maxsearch*Zsearch*AddSearch*0];
		massy.Zfudge   =testym[cnt+maxsearch*Zsearch*AddSearch*1];
		massy.adduction=testym[cnt+maxsearch*Zsearch*AddSearch*2];

		double chi2=runchap("fitty",trial,lines,data,"Null",params,x_init,massy,0);  // skewed gaussian with client=0 seperate + all params
		if(chi2!=chi2)
		  chi2=1E6;
		testym[cnt+maxsearch*Zsearch*AddSearch*3]=chi2;
		fp=fopen(labby.c_str(),"a");
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",testym[cnt+maxsearch*Zsearch*AddSearch*0],testym[cnt+maxsearch*Zsearch*AddSearch*1],testym[cnt+maxsearch*Zsearch*AddSearch*2],chi2,x_init[0],x_init[1],x_init[2],x_init[3]);
		fclose(fp);
		cnt++;
	      }
	  }
	fp=fopen(labby.c_str(),"a");
	fprintf(fp,"\n");
	fclose(fp);
      }
    int imin=findmin(testym,maxsearch*Zsearch*AddSearch,3);//find the lowest chi2 element in the test array
    struct mass massy;//run the final fit with these parameters
    massy=mass_init;
    x_init[0]      =testym[imin+maxsearch*Zsearch*AddSearch*0];   //update initialise matrix with minimsed value
    massy.Zfudge   =testym[imin+maxsearch*Zsearch*AddSearch*1];
    massy.adduction=testym[imin+maxsearch*Zsearch*AddSearch*2];


    cout << " Best fitting oligomer: " << x_init[0] << endl;
    cout << " Optimum Zfudge       : " << massy.Zfudge << endl;
    cout << " Optimum adduction    : " << massy.adduction << endl;
    double chi2=runchap("fitty",trial,lines,data,"Null",params,x_init,massy,0);  // skewed gaussian with client=0 seperate + all params
    //update the struct with the new fitted values ready for main minimisation. Minimiser should start very close to the bottom of the well
    mass_init.MRes      =massy.MRes;
    mass_init.adduction =massy.adduction;
    mass_init.Zfudge    =massy.Zfudge;
    mass_init.ResFudge  =massy.ResFudge;

    return;
  }

    */


//pack up current pars and save jiggle condition
void PackJiggle(double *parsCurr,mass &spec)
{
  spec.Pack();
  for(int i=0;i<spec.fitpar.p;++i)
    parsCurr[i]=spec.par[i];
}

//unpack pars and pour them into the spec class
void UnpackJiggle(double *parsCurr,mass &spec)
{
  for(int i=0;i<spec.fitpar.p;++i)
    spec.par[i]=parsCurr[i];
  spec.UnPack();
}

//add noise to parsNew from parsCurr
void AddNoise(double *parsNew,double *parsCurr,double sigma,gsl_rng *r,int parNo)
{
  for(int i=0;i<parNo;++i)
    {
      double rando=(gsl_ran_gaussian(r,sigma)+1);
      parsNew[i]=parsCurr[i]*rando;
      cout << "   rando: " << rando << " oldpar: " << parsCurr[i] << " newpar: " << parsNew[i] << endl;
    }
}

//master function to control the jiggler
void Jiggler(int flg,mass &spec,double maxIter,int jiggles,double jiggleSigma)
{
  fitty(0,spec,maxIter); //run a fitty to setup the jiggler
  spec.ShowPars();  //show the initial parameters

  //initialise random number generator
  const gsl_rng_type *T;
  gsl_rng *r;
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);

  spec.run_calc(0); //get starting chi2
  spec.CalcChi2(); //numerically calculate chi2
  double chi2=spec.chi2/spec.lines;
  cout << "Unleashing the jiggler!" << endl;
  cout << "Starting chi2: " << chi2 << endl;

  int p=spec.fitpar.p;

  double parsCurr[p];
  double parsNew[p];
  PackJiggle(parsCurr,spec);

  for(int i=0;i<p;++i) //show starting pars
      cout << "   ->   Startpar " << parsCurr[i] << endl;

  int go=0;   //control the while loop
  int cnt=0; //control the while loop
  while(go==0) //keep going until go is broken
    {
      AddNoise(parsNew,parsCurr,jiggleSigma,r,p); //add noise to parameteres
      UnpackJiggle(parsNew,spec);  //unpack the new parameters into the spectrum.
      fitty(0,spec,maxIter); //optimise pars
      spec.CalcChi2(); //numerically calculate chi2
      double chi2new=spec.chi2/spec.lines;
      cout << "oldchi2 "<< chi2 << " newchi2 " << chi2new << " counts " << cnt << endl;

      if(chi2new<chi2 && isnan(chi2new)==0) //keep new one if it's any good
	{
	  cnt=0;//reset counter
	  PackJiggle(parsCurr,spec); //store parameters
	  //if( fabs(chi2new-chi2)>0.01)
	  cout << "  yay! Have lowerered chi2 from " << chi2 << " to " << chi2new << endl;
	  cout << "  current zfudge: " << spec.Zfudge << endl;
	  cout << endl;
	  spec.ShowPars();
	  cout << endl;
	  chi2=chi2new;
	}
      cnt++; //increment counter
      if(cnt==jiggles) //if it's big enough, break loop.
	go=1;
    }
  gsl_rng_free(r);

  UnpackJiggle(parsCurr,spec); //unpack best parameters
  spec.run_calc(1);//sim and save
  cout << "final jiggler chi2: " << spec.chi2/spec.lines << endl;
  //cout << "final Zfudge " << spec.Zfudge << endl;
}



//analysis class to run the calculation.
class anal
{
 public:

  vector<mass> files;   //will store as many spectra as needed for the calculation
  int jiggles=20;         //default jiggler convergence condition
  double jiggleSigma=0.5; //default jiggler noise condition
  double maxIter=50;      //default steepest descent interaction number

  //class to store the specific optimisation protocol
  class protocol
  {
  public:
    string type;
    vector<string> pars;
    void AddStep(vector<string> liney)
    {
      type=liney[0];
      cout << "adding protocol line: " << endl;
      cout << type << endl;
      cout << liney.size() << endl;
      for(int i=0;i<liney.size();++i)
	cout << liney[i] << " " ;
      cout << endl;
      for(int i=0;i<liney.size()-2;++i)
	pars.push_back(liney[i+2]);

    }
  };

  vector< vector<protocol> > protocols; //store the required calculation

  void ParseInpFile(string inputfile) //parse input file
  {
    cout << endl << "Reading inputfile: " << inputfile << endl;
    vector<vector<string> > infile;
    infile=MakeFileVec(inputfile.c_str());

    char tost='#';

    int prot=0;
    for(int i=0;i<infile.size();++i)
	if(infile[i].size()>0)
	    if(infile[i][0]=="PROTOCOL")
	      prot++;
    cout << "Protocols to follow: " << prot << endl;

    for(int i=0;i<prot;++i) //setup 'files' vector with required number of runs.
      {
	mass spec;
	files.push_back(spec);
	vector<protocol> pros;
	protocols.push_back(pros);
      }

    int prot_flg=0;  //have we reach the protocol divergence statements yet?
    //monster loop to read entries in input file and load into relevant place
    for(int i=0;i<infile.size();++i) //now parse the file and add parameters
      {
	if(infile[i].size()>0 && infile[i][0][0]!=tost)
	  {
	    if(infile[i][0]=="raw")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].raw=infile[i][2];
	      else
		files[prot_flg-1].raw=infile[i][2];

	    else if(infile[i][0]=="ident")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].ident=infile[i][2];
	      else
		files[prot_flg-1].ident=infile[i][2];

	    else if(infile[i][0]=="type")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].type=infile[i][2];
	      else
		files[prot_flg-1].type=infile[i][2];

	    else if(infile[i][0]=="gnuplot")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].gnuplot=infile[i][2];
	      else
		files[prot_flg-1].gnuplot=infile[i][2];



	    else if(infile[i][0]=="thresh")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		    files[j].thresh=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].thresh=atof(infile[i][2].c_str());


	    else if(infile[i][0]=="limHSP")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].limHSP=atoi(infile[i][2].c_str());
	      else
		files[prot_flg-1].limHSP=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="limSub")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].limSub=atoi(infile[i][2].c_str());
	      else
		files[prot_flg-1].limSub=atoi(infile[i][2].c_str());

	    else if(infile[i][0]=="dim")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].dim=atoi(infile[i][2].c_str());
	      else
		files[prot_flg-1].dim=atoi(infile[i][2].c_str());

	    else if(infile[i][0]=="shspMass")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].shspMass=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].shspMass=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="subMass")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].subMass=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].subMass=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="minMZ")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].minMZ=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].minMZ=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="maxMZ")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].maxMZ=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].maxMZ=atof(infile[i][2].c_str());


	    else if(infile[i][0]=="smooth")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].smooth=atoi(infile[i][2].c_str());
	      else
		files[prot_flg-1].smooth=atoi(infile[i][2].c_str());

	    else if(infile[i][0]=="Zfudge")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].Zfudge=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].Zfudge=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="adduction")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].adduction=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].adduction=atof(infile[i][2].c_str());


	    else if(infile[i][0]=="Zwidth")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].Zwidth=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].Zwidth=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="MRes")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].MRes=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].MRes=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="ResFudge")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].ResFudge=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].ResFudge=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="cWid")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].cWid=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].cWid=atof(infile[i][2].c_str());

	    else if(infile[i][0]=="zWid")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  files[j].zWid=atof(infile[i][2].c_str());
	      else
		files[prot_flg-1].zWid=atof(infile[i][2].c_str());


	    else if(infile[i][0]=="ReadFree")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  {
		    files[j].input_file=infile[i][2];
		    files[j].input_read=1;
		  }
	      else
		{
		  files[prot_flg-1].input_file=infile[i][2];
		  files[prot_flg-1].input_read=1;
		}

	    else if(infile[i][0]=="AddDist")
	      if(prot_flg==0)
		for(int j=0;j<prot;++j)
		  {
		    files[j].AddDist(infile[i][2]);
        for(int k=0; k<files[j].distList[files[j].distList.size()-1].p;++k){
		      files[j].distList[files[j].distList.size()-1].pars[k]=atof(infile[i][k+3].c_str());
        }
		  }
	      else
		{
		  files[prot_flg-1].AddDist(infile[i][2]);
      for(int k=0; k<files[prot_flg-1].distList[files[prot_flg-1].distList.size()-1].p;++k){
        files[prot_flg-1].distList[files[prot_flg-1].distList.size()-1].pars[k]=atof(infile[i][k+3].c_str());
      }


		}


	    else if(infile[i][0]=="sim" || infile[i][0]=="fitty" || infile[i][0]=="jiggle")
	      {

    //cout << prot_flg << infile[i][0] << endl;
		protocol pro;
		pro.AddStep(infile[i]);
		if(prot_flg==0)
		  for(int j=0;j<prot;++j)
		    {
		      protocols[j].push_back(pro);
		    }
		else
		  protocols[prot_flg-1].push_back(pro);

	      }

	    else if(infile[i][0]=="maxIter")
	      {
		maxIter=atof(infile[i][2].c_str());
	      }

	    else if(infile[i][0]=="jiggles")
	      {
		jiggles=atoi(infile[i][2].c_str());
	      }

	    else if(infile[i][0]=="jiggleSigma")
	      {
		jiggleSigma=atof(infile[i][2].c_str());
	      }

	    else if(infile[i][0]=="PROTOCOL")
	      {
		prot_flg++;  //increment protocol
	      }


	    else
	      {
		cout << "Unparsed command: " << infile[i][0] << endl;
		cout << "line " << i+1 << endl;
		exit(100);
	      }


	  }
      }
    cout << endl;
    cout << "Calculation plan: " << endl;
    cout << endl << "Files to analyse: " << files.size() << endl;
    for(int i=0;i<protocols.size();++i)
      {
	cout << files[i].ident << " " << files[i].raw << " " << files[i].dim << endl;
	for(int j=0;j<protocols[i].size();++j)
	  {
	  cout << "step " << j+1 << " : " << protocols[i][j].type << endl;
	  cout << "     flags: ";
	  for(int k=0;k<protocols[i][j].pars.size();++k)
	    cout << protocols[i][j].pars[k] << " ";
	  cout << endl;
	  }

      }
    cout << endl;
  }

  //execute a calculation plan
  void RunProtocol(int i)
  {
    for(int j=0;j<protocols[i].size();++j)
      {
	cout << "Running protocol step: " << protocols[i][j].type << " on file " << files[i].raw << endl;
	if(protocols[i][j].type=="sim")
	  {
	    runchap("sim",i,1);   //run minimisation (j specifies skew or not)
	  }
	else if(protocols[i][j].type=="fitty")
	  {
	    files[i].fitpar.SetPars(protocols[i][j].pars,files[i].distList,files[i].freeList); //reset flags
	    files[i].fitpar.ShowFlags();
	    runchap("fitty",i,1); //run minimisation (j specifies skew or not)
	  }

	else if(protocols[i][j].type=="jiggle")
	  {
	    files[i].fitpar.SetPars(protocols[i][j].pars,files[i].distList,files[i].freeList); //reset flags
	    files[i].fitpar.ShowFlags();
	    runchap("jiggle",i,1); //run minimisation (j specifies skew or not)
	  }


      }
  }


  //function execute specific functions sim fitty or jiggler
  void runchap(string mode,int i,int verb)
  {
    cout << "Running mode: " << mode << endl;

    //Do this only if the outfile isn't set to 'Null'
    if(verb)
      SetTrial(i);
    files[i].mode=mode;

    if(mode=="fitty")
      fitty(1,files[i],maxIter);
    else if(mode=="sim")
      files[i].run_calc(1);
    else if(mode=="jiggle")
      Jiggler(1,files[i],maxIter,jiggles,jiggleSigma);
    else
      {
	cout << "Mode not recognised. Aborting." << endl;
	cout  << mode << endl;
	exit(100);
      }
    files[i].trial++; //increment trial number
  }

  //convert trial number into a tag, _a, _b etc.
  void SetTrial(int i)
  {
    string lab;
    //append tag depending on run number
    switch(files[i].trial){
    case 0:
      lab="_a";
      break;
    case 1:
      lab="_b";
      break;
    case 2:
      lab="_c";
      break;
    case 3:
      lab="_d";
      break;
    case 4:
      lab="_e";
      break;
    case 5:
      lab="_f";
      break;
    case 6:
      lab="_g";
      break;
    }
    files[i].identTag=files[i].ident+lab; //create the identity tag from ident and label
  }


};
