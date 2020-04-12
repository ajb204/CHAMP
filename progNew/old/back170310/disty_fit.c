struct data {
  int chap;
  int p;
  double * datum;
  int lines;
  int limHSP;
  int limSub;
  double adduction;
  double Zwidth;
  double shspMass;
  double subMass;
  int testmax;
  int minZ;
  int maxZ;
  double thresh;
  double MRes;
};


int
expb_f (const gsl_vector * x, void *data, 
	gsl_vector * f)
{
  int chap=((struct data *)data)->chap;
  int lines=((struct data *)data)->lines;
  double* datum=((struct data *)data)->datum;

  int limHSP=((struct data *)data)->limHSP;
  int limSub=((struct data *)data)->limSub;
  double adduction=((struct data *)data)->adduction;
  double Zwidth=((struct data *)data)->Zwidth;
  double shspMass=((struct data *)data)->shspMass;
  double subMass=((struct data *)data)->subMass;
  int testmax=((struct data *)data)->testmax;
  int minZ=((struct data *)data)->minZ;
  int maxZ=((struct data *)data)->maxZ;
  double thresh=((struct data *)data)->thresh;
  double MRes=((struct data *)data)->MRes;

  double Hx0 = gsl_vector_get (x, 0);
  double Hsig= gsl_vector_get (x, 1);
  double Sx0 = (gsl_vector_get (x, 2));
  double Ssig= gsl_vector_get (x, 3);
  double alpha=0;
  if(chap==1)
    alpha= gsl_vector_get (x, 4);
  //calculate simulated spectrum with the current parameters
  double par[10];

  par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;
  if(chap==1)
    par[4]=alpha;
  
  double chi2=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes);

  //printf("chi2 calc called from f(x) %f\n",chi2);
  //printf("Hx0 %f\tHsig %f\tSx0 %f\tSsig %f alpha %f\n",Hx0,Hsig,Sx0,Ssig,alpha);

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
  int p=((struct data *)data)->p;
  int lines=((struct data *)data)->lines;
  double* datum=((struct data *)data)->datum;

  int limHSP=((struct data *)data)->limHSP;
  int limSub=((struct data *)data)->limSub;
  double adduction=((struct data *)data)->adduction;
  double Zwidth=((struct data *)data)->Zwidth;
  double shspMass=((struct data *)data)->shspMass;
  double subMass=((struct data *)data)->subMass;
  int testmax=((struct data *)data)->testmax;
  int minZ=((struct data *)data)->minZ;
  int maxZ=((struct data *)data)->maxZ;
  double thresh=((struct data *)data)->thresh;
  double MRes=((struct data *)data)->MRes;


  double STEP=1E-6;
  
  double Hx0 = gsl_vector_get (x, 0);
  double Hsig= gsl_vector_get (x, 1);
  double Sx0 = (gsl_vector_get (x, 2));
  double Ssig= gsl_vector_get (x, 3);
  double alpha=0;
  if(chap==1)
    alpha= gsl_vector_get (x, 4);

  double par[10];
  par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;
  if(chap==1)
    par[4]=alpha;
  double chi2_0=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes);
  double datum0[lines*3];for (int i=0;i<lines*3;i++) datum0[i]=datum[i];

  // printf("chi2 calc called from df(x) %f\n",chi2_0);
  //printf("Hx0 %f\tHsig %f\tSx0 %f\tSsig %f alpha %f\n",Hx0,Hsig,Sx0,Ssig,alpha);

  for(int k=0;k<p;k++)
    {
  //do kth parameter
      par[0]=Hx0;par[1]=Hsig;par[2]=Sx0;par[3]=Ssig;par[4]=alpha;
      par[k]=par[k]+STEP;

      double chi2_1=run_calc(0,chap,"Null",0,limHSP,limSub,par,adduction,Zwidth,shspMass,subMass,testmax,minZ,maxZ,thresh,datum,lines,MRes);
      for (int i = 0; i < lines; i++)
        //gradient will be Y(p)-Y(p+dp)
        gsl_matrix_set (J, i, k, -1.0*( datum[i+lines*2]-datum0[i+lines*2] )/(STEP*10));
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
  for(int i=0;i<p;i++)
    printf("% 15.8f ",gsl_vector_get(s->x,i));
  printf(" |f(x)| = %g\n", gsl_blas_dnrm2 (s->f));

}



void
fitty (char* outfile,struct data spec,double *x_init,int pp)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = spec.lines;
  const size_t p = pp;
  spec.p=p;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;

  double x_init2[p];
  for (int i=0;i<p;i++)
      x_init2[i]=x_init[i];
  
  gsl_vector_view x = gsl_vector_view_array (x_init2, p);
  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &spec;
  
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
    

    printf("\nCovariance matrix\n");
    for (int i=0;i<p;i++)
      {
	for(int j=0;j<p;j++)
	  printf("%f\t",c*gsl_matrix_get(covar,i,j));
	printf("\n");
      }
    

    printf ("Hx0      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    printf ("Hsig     = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    printf ("Sx0      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
    printf ("Ssig     = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
    if(spec.chap==1)
      printf ("alpha     = %.5f +/- %.5f\n", FIT(4), c*ERR(4));
    
    double par[10];
    par[0]=FIT(0);par[1]=FIT(1);par[2]=FIT(2);par[3]=FIT(3);
    if(spec.chap==1)
      par[4]=FIT(4);
    printf("chap = %i\t%f\t%f\n",spec.chap,par[0],par[4]);
    double chi2=run_calc(1,spec.chap,outfile,6,spec.limHSP,spec.limSub,par,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes);
    printf ("internal chi2/dof %f\n",chi2);
  }
  
  printf ("status = %s\n", gsl_strerror (status));
  
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  return;
}


