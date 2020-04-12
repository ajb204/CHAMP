

struct data {
  int chap;
  int pp;
  double * x_init;
  double * datum;
  int lines;
  int limHSP;
  int limSub;
  double Zwidth;
  double shspMass;
  double subMass;
  int testmax;
  int minZ;
  int maxZ;
  double thresh;
  double adduction;
  int adduction_flg;
  double Zfudge;
  int Zfudge_flg;
  double MRes;
  int MRes_flg;
  double ResFudge;
  int ResFudge_flg;
};




int
expb_f (const gsl_vector * x, void *data, 
	gsl_vector * f)
{
  int chap=((struct data *)data)->chap;
  int pp=((struct data *)data)->pp;

  int lines=((struct data *)data)->lines;
  double* x_init=((struct data *)data)->x_init;
  double* datum=((struct data *)data)->datum;

  int limHSP=((struct data *)data)->limHSP;
  int limSub=((struct data *)data)->limSub;
  double Zwidth=((struct data *)data)->Zwidth;
  double shspMass=((struct data *)data)->shspMass;
  double subMass=((struct data *)data)->subMass;
  int testmax=((struct data *)data)->testmax;
  int minZ=((struct data *)data)->minZ;
  int maxZ=((struct data *)data)->maxZ;
  double thresh=((struct data *)data)->thresh;
  double adduction=((struct data *)data)->adduction;
  int adduction_flg=((struct data *)data)->adduction_flg;
  double Zfudge=((struct data *)data)->Zfudge;
  int Zfudge_flg=((struct data *)data)->Zfudge_flg;
  double MRes=((struct data *)data)->MRes;
  int MRes_flg=((struct data *)data)->MRes_flg;
  double ResFudge=((struct data *)data)->ResFudge;
  int ResFudge_flg=((struct data *)data)->ResFudge_flg;




  double par[limHSP*limSub];    

  {
    int i=0;
    if(pp>0)//if we're minimsing the distribution, transfer params
      for(i=0;i<pp;i++)
	par[i]=gsl_vector_get (x, i);
    else //otherwise take the default params
      for(int j=0;j<10;j++)
	par[j]=x_init[j];
    if(adduction_flg==1) //if we're minimising adduction, take the param
      {adduction=gsl_vector_get (x, i);i++;}
    if(Zfudge_flg==1)    //if we're minimising Zfudge, take the param
      {Zfudge=gsl_vector_get (x, i);i++;}
    if(MRes_flg==1)    //if we're minimising MRes, take the param
      {MRes=gsl_vector_get (x, i);i++;}
    if(ResFudge_flg==1)    //if we're minimising ResFudge, take the param
      {ResFudge=gsl_vector_get (x, i);i++;}

  }



  double chi2=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge,ResFudge);
  //printf("chi2 calc called from f(x) %f\n",chi2);
  //printf("Hx0 %f\tHsig %f\tSx0 %f\tSsig %f alpha %f\n",par[0],par[1],par[2],par[3],par[4]);
  //cout << "f(x) MRes " << MRes << " Zfudge " << Zfudge << " adduction " << adduction << " ResFudge " << ResFudge << " chi2 " << chi2 << endl;

  for (size_t i = 0; i < lines; i++)
    //chi^2(i) = (Ydata(i)-Ysim(i))**2/sigma**2
    gsl_vector_set (f, i,  (datum[i+lines*1]-datum[i+lines*2])/10 );
    
  
  return GSL_SUCCESS;
}



int
expb_df (const gsl_vector * x, void *data, 
	 gsl_matrix * J)
{
  //cout << "CALCULATING DERIVATIES" << endl;
  int chap=((struct data *)data)->chap;  
  int pp=((struct data *)data)->pp;
  int lines=((struct data *)data)->lines;
  double* datum=((struct data *)data)->datum;
  double* x_init=((struct data *)data)->x_init;

  int limHSP=((struct data *)data)->limHSP;
  int limSub=((struct data *)data)->limSub;
  double Zwidth=((struct data *)data)->Zwidth;
  double shspMass=((struct data *)data)->shspMass;
  double subMass=((struct data *)data)->subMass;
  int testmax=((struct data *)data)->testmax;
  int minZ=((struct data *)data)->minZ;
  int maxZ=((struct data *)data)->maxZ;
  double thresh=((struct data *)data)->thresh;
  double adduction=((struct data *)data)->adduction;
  int adduction_flg=((struct data *)data)->adduction_flg;
  double Zfudge=((struct data *)data)->Zfudge;
  int Zfudge_flg=((struct data *)data)->Zfudge_flg;
  double MRes=((struct data *)data)->MRes;
  int MRes_flg=((struct data *)data)->MRes_flg;  
  double ResFudge=((struct data *)data)->ResFudge;
  int ResFudge_flg=((struct data *)data)->ResFudge_flg;  

  double STEP=1E-6;
  

  double par[limHSP*limSub];    
  {
    int i=0;
    if(pp>0)//if we're minimsing the distribution, transfer params
      for(i=0;i<pp;i++)
	par[i]=gsl_vector_get (x, i);
    else //otherwise take the default params
      for(int j=0;j<10;j++)
	par[j]=x_init[j];
    if(adduction_flg==1) //if we're minimising adduction, take the param
      {adduction=gsl_vector_get (x, i);i++;}
    if(Zfudge_flg==1)    //if we're minimising Zfudge, take the param
      {Zfudge=gsl_vector_get (x, i);i++;}
    if(MRes_flg==1)    //if we're minimising Zfudge, take the param
      {MRes=gsl_vector_get (x, i);i++;}
    if(ResFudge_flg==1)    //if we're minimising ResFudge, take the param
      {ResFudge=gsl_vector_get (x, i);i++;}
  }

  
  double chi2_0=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge,ResFudge);
  double datum0[lines*3];for (int i=0;i<lines*3;i++) datum0[i]=datum[i];
  //cout << "df(x) MRes " << MRes << " Zfudge " << Zfudge << " adduction " << adduction << " ResFudge " << ResFudge <<  " chi2 " << chi2_0 << endl;
  //printf("chi2 calc called from df(x) %f\n",chi2_0);
  //printf("Hx0 %f\tHsig %f\tSx0 %f\tSsig %f alpha %f\n",par[0],par[1],par[2],par[3],par[4]);

  //if we're optimising the distribution:
  {
    int i=0;
    if(pp>0)
      for (i=0;i<pp;i++)
	{
	  for(int k=0;k<pp;k++)//reset distribution parameters
	    par[k]=gsl_vector_get (x, k);
	  par[i]=par[i]+STEP;
	  double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge,ResFudge);
	  //printf("Hx0 %f\tHsig %f\tSx0 %f\tSsig %f alpha %f %f %f\n",par[0],par[1],par[2],par[3],par[4],chi2_1,chi2_0);
	  for (int j = 0;j < lines; j++)
	    //gradient will be Y(p)-Y(p+dp)
	    gsl_matrix_set (J, j, i, -1.0*( datum[j+lines*2]-datum0[j+lines*2] )/(STEP*10));
	}
  
    if(adduction_flg==1)
      {
	for(int k=0;k<pp;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction+STEP,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge,ResFudge);
	for (int j = 0;j < lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( datum[j+lines*2]-datum0[j+lines*2] )/(STEP*10));
	i++;
      }
    
    if(Zfudge_flg==1)
      {
	for(int k=0;k<pp;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge+STEP,ResFudge);
	for (int j = 0;j < lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( datum[j+lines*2]-datum0[j+lines*2] )/(STEP*10));
	i++;
      }
    if(MRes_flg==1)
      {
	for(int k=0;k<pp;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes+STEP,Zfudge,ResFudge);
	for (int j = 0;j < lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( datum[j+lines*2]-datum0[j+lines*2] )/(STEP*10));
	i++;
      }
    if(ResFudge_flg==1)
      {
	for(int k=0;k<pp;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes,Zfudge,ResFudge+STEP);
	for (int j = 0;j < lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( datum[j+lines*2]-datum0[j+lines*2] )/(STEP*10));
	i++;
      }
  }
  
  return GSL_SUCCESS;
}


int
expb_fdf (const gsl_vector * x, void *data,
	  gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, data, f);
  expb_df (x,data, J);
  
  return GSL_SUCCESS;
}


void
print_state (size_t iter, gsl_multifit_fdfsolver * s,int p)
{
    printf ("iter: %3u x = ",iter);
  if(p<10)
     for(int i=0;i<p;i++)
      printf("% 15.8f ",gsl_vector_get(s->x,i));
  printf(" |f(x)| = %g\n", gsl_blas_dnrm2 (s->f));

}



void
fitty (char* outfile,struct data& spec,int pp)
{

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = spec.lines;
  const size_t p = (pp+spec.adduction_flg+spec.Zfudge_flg+spec.MRes_flg+spec.ResFudge_flg);
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  double x_init2[p];


  if(spec.chap==5 || spec.chap==6)
    {
      int i=0;

      double input[spec.limHSP*spec.limSub];init_array(input,spec.limHSP*spec.limSub);
      make_input(input,spec.limHSP,spec.limSub,spec.chap-5,spec.x_init);   //calculate the 2d skewed gaussian distribution
      //prin_input("testy1.inp","",0,input,spec.limHSP,spec.limSub);    //print input matrix
      
      double Hx0 =spec.x_init[0];  //extract distribution parameters
      double Hsig=spec.x_init[1]; //extract distribution parameters
      double Sx0 =spec.x_init[2];  //extract distribution parameters
      double Ssig=spec.x_init[3]; //extract distribution parameters
      double alpha=0;
      if(spec.chap==6)
	alpha=spec.x_init[4]; //extract distribution parameters

      //double chi2=run_calc(0,0,"Null",0,spec.limHSP,spec.limSub,spec.x_init,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes,spec.Zfudge,spec.ResFudge);
      

      int ii=int(Hx0);//take nearest integer to Sx0      
      //set initial parameters for fitting
      for(int j=0;j<spec.limSub;j++)
	{

	  if(spec.chap==5){
	    x_init2[i]=Hx0;i++;
	    x_init2[i]=Hsig;i++;
	    x_init2[i]=input[ii+spec.limHSP*j];i++;
	    printf ("%i Hx0_0      = %.5f \n",j, x_init2[0+3*j]);
	    printf ("%i Hsig_0     = %.5f \n",j, x_init2[1+3*j]);
	    printf ("%i rel        = %.5f \n",j, x_init2[2+3*j]);}
	  else{
	    x_init2[i]=Hx0;i++;
	    x_init2[i]=Hsig;i++;
	    x_init2[i]=alpha;i++;
	    x_init2[i]=input[ii+spec.limHSP*j];i++;
	    printf ("%i Hx0_0     = %.5f \n", j,x_init2[0+4*j]);
	    printf ("%i Hsig_0    = %.5f \n", j,x_init2[1+4*j]);
	    printf ("%i alpha     = %.5f \n", j,x_init2[2+4*j]);
	    printf ("%i rel       = %.5f \n", j,x_init2[3+4*j]);}
	}
      
      if(spec.adduction_flg==1)
	{
	  cout << "Fitting adduction" << spec.adduction <<endl;
	  x_init2[i]=spec.adduction;i++;
	}
      if(spec.Zfudge_flg==1)
	{
	  cout << "Fitting Zfudge" << spec.Zfudge <<endl;
	  x_init2[i]=spec.Zfudge;i++;
	}
      if(spec.MRes_flg==1)
	{
	  cout << "Fitting MRes "  << spec.MRes << endl;
	  x_init2[i]=spec.MRes;i++;
	}
      if(spec.ResFudge_flg==1)
	{
	  cout << "Fitting ResFudge "  << spec.ResFudge <<endl;
	  x_init2[i]=spec.ResFudge;i++;
	}


      //init_array(input,spec.limHSP*spec.limSub);
      //make_input(input,spec.limHSP,spec.limSub,5,x_init2);   //calculate the 2d skewed gaussian distributio
      //prin_input("testy2.inp","",0,input,spec.limHSP,spec.limSub);    //print input matrix
      //double chi22=run_calc(0,5,"Null",0,spec.limHSP,spec.limSub,x_init2,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes,spec.Zfudge,spec.ResFudge);
      //cout << chi22 << endl;
      
      
    }
  
  
  
  
  if(spec.chap<10 && spec.chap!=6 && spec.chap!=5)
  {//setup initial parameter matrix
    int i=0;
    if(pp>0){cout << "Fitting " << pp << " distribution parameters" << endl;
      for(i=0;i<pp;i++){
	x_init2[i]=spec.x_init[i];}}
    
    if(spec.adduction_flg==1)
      {
	cout << "Fitting adduction" << spec.adduction <<endl;
	x_init2[i]=spec.adduction;i++;
      }
    if(spec.Zfudge_flg==1)
      {
	cout << "Fitting Zfudge" << spec.Zfudge <<endl;
	x_init2[i]=spec.Zfudge;i++;
      }
    if(spec.MRes_flg==1)
      {
	cout << "Fitting MRes " << spec.MRes << endl;
	x_init2[i]=spec.MRes;i++;
      }
    if(spec.ResFudge_flg==1)
      {
	cout << "Fitting MRes " << spec.ResFudge << endl;
	x_init2[i]=spec.ResFudge;i++;
      }
  }

  if(spec.chap>=10)
    //completely free distribution parameters
    {
      int i=0;
      cout << "Running uber fit with " << pp << " parameters" << endl;
      double input[spec.limHSP*spec.limSub];init_array(input,spec.limHSP*spec.limSub);
      make_input(input,spec.limHSP,spec.limSub,spec.chap-10,spec.x_init);
      for (int j=0;j<5;j++)
	cout << spec.x_init[j] << endl;
      for(i=0;i<pp;i++)
	x_init2[i]=input[i];

      if(spec.adduction_flg==1)
	{
	  cout << "Fitting adduction" << spec.adduction <<endl;
	  x_init2[i]=spec.adduction;i++;
	}
      if(spec.Zfudge_flg==1)
	{
	  cout << "Fitting Zfudge" << spec.Zfudge <<endl;
	  x_init2[i]=spec.Zfudge;i++;
	}
      if(spec.MRes_flg==1)
	{
	  cout << "Fitting MRes "  << spec.MRes << endl;
	  x_init2[i]=spec.MRes;i++;
	}
      if(spec.ResFudge_flg==1)
	{
	  cout << "Fitting ResFudge "  << spec.ResFudge <<endl;
	  x_init2[i]=spec.ResFudge;i++;
	}
      
    }
  
  

  
  gsl_vector_view x = gsl_vector_view_array (x_init2, p);      

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &spec;

  cout << "Fitting with "<< p << " total parameters " << endl;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
  
  print_state (iter, s,p);
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
      //printf ("status = %s\n", gsl_strerror (status));
      print_state (iter, s,p);
      if (status)
	break;
      status = gsl_multifit_test_delta (s->dx, s->x,
					1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 50);
  
  gsl_multifit_covar (s->J, 0.0, covar);
  
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
  
  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
    
    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
    

    if(spec.chap<10 && spec.chap!=6 && spec.chap!=5)
      {
	printf("\nCovariance matrix\n");
	for (int i=0;i<pp;i++)
	  {
	    for(int j=0;j<pp;j++)
	      printf("%f\t",c*gsl_matrix_get(covar,i,j));
	    printf("\n");
	  }
      }
    
    

    double par[spec.limHSP*spec.limSub];
    {
      int ii=0;
      
      if(spec.chap<10 && pp>0 && spec.chap!=5 && spec.chap!=6)
	{
	  printf ("Hx0      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	  printf ("Hsig     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	  printf ("Sx0      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	  printf ("Ssig     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	  if(spec.chap==1){
	    printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  if(spec.chap==2){
	    printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("skew      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  if(spec.chap==3){
	    printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("alpha2    = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  if(spec.chap==4){
	    printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("Hx0_0     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("Hsig_0    = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("alpha_0   = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("rel       = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  
	  for (int j=0;j<pp;j++)
	    par[j]=spec.x_init[j];
	  
	}
      
      
      if(spec.chap==6)
	{
	  
	  for(int j=0;j<spec.limSub;j++){
	    printf ("%i Hx0_0     = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("%i Hsig_0    = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("%i alpha_0   = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("%i rel       = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  
	  for (int j=0;j<pp;j++)
	    par[j]=spec.x_init[j];
	}


      if(spec.chap==5)
	{
	  
	  for(int j=0;j<spec.limSub;j++){
	    printf ("%i Hx0_0     = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("%i Hsig_0    = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
	    printf ("%i alpha_0   = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
	  
	  for (int j=0;j<pp;j++)
	    par[j]=spec.x_init[j];
	}
      
      
      
      
      if(spec.chap<10 && pp==0)//do this if not fitting distribution at all
	for (int j=0;j<pp;j++)
	  par[j]=spec.x_init[j];
      
      
      if(spec.chap>=10 && pp>0)//do this if doing the 'all param' dist fit
	for (ii=0;ii<pp;ii++)
	  par[ii]=FIT(ii);
      
      
      
      
      if(spec.adduction_flg==1){
	printf ("Adduction= %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.adduction=FIT(ii);ii++;  }
      if(spec.Zfudge_flg==1){
	printf ("Zfudge   = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.Zfudge=FIT(ii);ii++;}
      if(spec.MRes_flg==1){
	printf ("Mres     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.MRes=FIT(ii);ii++;}
      if(spec.ResFudge_flg==1){
	printf ("ResFudge = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.ResFudge=FIT(ii);ii++;}
      
    }
    
        
    // printf("chap = %i\t%f\t%f\n",spec.chap,par[0],par[4]);
    double chi2=run_calc(1,spec.chap,outfile,6,spec.limHSP,spec.limSub,par,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes,spec.Zfudge,spec.ResFudge);
    printf ("internal chi2/dof %f\n",chi2);
  }
  
  printf ("status = %s\n", gsl_strerror (status));
  
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  return;
}
  



  void runchap(string infile,char* outfile,int chap,int distparam,int adduction_flg,int Zfudge_flg,int MRes_flg,int ResFudge_flg,double* x_init,struct mass& massy,double subMass)
{
  if(distparam!=0)
    {
      if(chap==0)
	distparam=4;
      if(chap==1)
	distparam=5;
      if(chap==2)
	distparam=6;
      if(chap==3)
	distparam=6;
      if(chap==4)
	distparam=9;
      if(chap==5)
	distparam=massy.limSub*3;
      if(chap==6)
	distparam=massy.limSub*4;
      if(chap==7)
	distparam=5;
      if(chap>=10)
	distparam=massy.limSub*massy.limHSP;
    }
  
  int lines=countlines_string(infile);

  //lines=lines+1000;
  double data[(lines+1)*3];init_array(data,3*(lines+1));    
  readfile(data,lines,infile);            //load in input data
  prin_spec(infile,".out",4,data,lines,1);//print out input spectrum
  struct data spec={chap,distparam,x_init,data,lines,massy.limHSP,massy.limSub,massy.Zwidth,massy.shspMass,subMass,massy.testmax,massy.minZ,massy.maxZ,massy.thresh,massy.adduction,adduction_flg,massy.Zfudge,Zfudge_flg,massy.MRes,MRes_flg,massy.ResFudge,ResFudge_flg};


  fitty(outfile,spec,distparam);

  //cout << spec.adduction << " " << massy.adduction << endl;
  //cout << spec.Zfudge << " " << massy.Zfudge << endl;
  //cout << spec.MRes << " " << massy.MRes << endl;

  massy.adduction=spec.adduction;
  massy.Zfudge=spec.Zfudge;
  massy.MRes=spec.MRes;
  massy.ResFudge=spec.ResFudge;


  //cout << spec.adduction << " " << massy.adduction << endl;
  //cout << spec.Zfudge << " " << massy.Zfudge << endl;
  //cout << spec.MRes << " " << massy.MRes << endl;

  return;
}
