
//STRUCT FOR HOLDING ALL THE RELEVANT INFO


/* CALCULATE CHI^2 MATRIX */
int
expb_f (const gsl_vector * x, void *data, 
		gsl_vector * f)
{

  struct data *spec=(struct data* )data;


  
  double par[(*spec).massy.limHSP*(*spec).massy.limSub];
  
  {
    int i=0;  //start incrementor for paramter number

    if( (*spec).params.chap==7){//unique to the single species fitting...
      par[0]=(*spec).x_init[0]; //get the principle oligomer...
      par[1]=(*spec).x_init[1]; //get the number of extra species...
      for(int j=0;j<par[1];j++)
	{
	  par[2*j+2]=(*spec).x_init[2*j+2]; //add the size of the next species...
	  //in par array, odd numbers that are not 1 are optimised
	  par[2*j+3]=gsl_vector_get (x, j); //get the concentration to optimise.
	  i++; //increment the distribution parameters
	}
    }
    else
      if((*spec).params.distparam>0)//if we're minimising the distribution, transfer params
	for(i=0;i<(*spec).params.distparam;i++)
	  {par[i]=gsl_vector_get (x, i);
	  }
      else //otherwise just take the default params
	for(int j=0;j<10;j++)
	  par[j]=(*spec).x_init[j];
    
    if((*spec).params.adduction_flg==1) //if we're minimising adduction, take the param
      {spec->massy.adduction=gsl_vector_get (x, i);i++;}
    if((*spec).params.Zfudge_flg==1)    //if we're minimising Zfudge, take the param
      {spec->massy.Zfudge=gsl_vector_get (x, i);i++;}
    if((*spec).params.MRes_flg==1)    //if we're minimising MRes, take the param
      {spec->massy.MRes=gsl_vector_get (x, i);i++;}
    if((*spec).params.ResFudge_flg==1)    //if we're minimising ResFudge, take the param
      {spec->massy.ResFudge=gsl_vector_get (x, i);i++;}
    if((*spec).params.tw_flg==1)    //if we're minimising residual 12mer, take the param
      {spec->massy.tw=gsl_vector_get (x, i);i++;}
    
  }

  double chi2=run_calc(0,"Null",0,par,*spec);
  
  for (size_t i = 0; i < (*spec).lines; i++)   //record and store chi values
    //chi(i) = (Ydata(i)-Ysim(i))/sigma
    gsl_vector_set (f, i,  ((*spec).datum[i+(*spec).lines*1]-(*spec).datum[i+(*spec).lines*2])/10 );
  
  return GSL_SUCCESS;
}



/* CALCULATE JACOBIAN MATRIX */
int
expb_df (const gsl_vector * x, void *data, 
		gsl_matrix * J)
{
  //cout << "CALCULATING DERIVATIES" << endl;
  struct data *spec=(struct data* )data;
  double STEP=1E-6;   //define step size for gradient differences
  
  double par[(*spec).massy.limHSP*(*spec).massy.limSub];
  {
    int i=0;  //start incrementor for paramter number

    if( (*spec).params.chap==7){//unique to the single species fitting...
      par[0]=(*spec).x_init[0]; //get the principle oligomer...
      par[1]=(*spec).x_init[1]; //get the number of extra species...
      for(int j=0;j<par[1];j++)
	{
	  par[2*j+2]=(*spec).x_init[2*j+2]; //add the size of the next species...
	  //in par array, odd numbers that are not 1 are optimised
	  par[2*j+3]=gsl_vector_get (x, j); //get the concentration to optimise.
	  i++;
	}
      //if( (*spec).params.distparam!=i){
      //cout << "PROBLEM WTIH PARAM COUNTING!" << endl;
      //exit(100);}
    }
    else
      if((*spec).params.distparam>0)//if we're minimising the distribution, transfer params
	for(i=0;i<(*spec).params.distparam;i++)
	  {par[i]=gsl_vector_get (x, i);
	  }
      else //otherwise just take the default params
	for(int j=0;j<10;j++)
	  par[j]=(*spec).x_init[j];


    

    if((*spec).params.adduction_flg==1) //if we're minimising adduction, take the param
      {spec->massy.adduction=gsl_vector_get (x, i);i++;}
    if((*spec).params.Zfudge_flg==1)    //if we're minimising Zfudge, take the param
      {spec->massy.Zfudge=gsl_vector_get (x, i);i++;}
    if((*spec).params.MRes_flg==1)    //if we're minimising MRes, take the param
      {spec->massy.MRes=gsl_vector_get (x, i);i++;}
    if((*spec).params.ResFudge_flg==1)    //if we're minimising ResFudge, take the param
      {spec->massy.ResFudge=gsl_vector_get (x, i);i++;}
    if((*spec).params.tw_flg==1)    //if we're minimising residual 12mer, take the param
      {spec->massy.tw=gsl_vector_get (x, i);i++;}
  }

  double chi2=run_calc(0,"Null",0,par,*spec);
  double datum0[(*spec).lines*3];for (int i=0;i<(*spec).lines*3;i++) datum0[i]=(*spec).datum[i];  //store the STEP=0 output for safe keeping...
  
  //if we're optimising the distribution:
  {
    int i=0;
    if((*spec).params.distparam>0)
      for (i=0;i<(*spec).params.distparam;i++)//loop over each distribution parameter (pp)
	{
	  if((*spec).params.chap!=7){
	    for(int k=0;k<(*spec).params.distparam;k++)//reset distribution parameters
	      {
		par[k]=gsl_vector_get (x, k);
	      }
	    par[i]=par[i]+STEP;
	  }
	  else
	    {//only increment the concentration value
	      par[0]=(*spec).x_init[0]; //get the principle oligomer...
	      par[1]=(*spec).x_init[1]; //get the number of extra species...
	      for(int j=0;j<par[1];j++){
		  par[2*j+2]=(*spec).x_init[2*j+2]; //add the size of the next species...
		  par[2*j+3]=gsl_vector_get (x, j); //get the concentration to optimise.
	      }
	      par[2*i+3]=par[2*i+3]+STEP*100; //increment distribution parameter
	    }
	  
	  
	  double chi2_1=run_calc(0,"Null",0,par,*spec);
	  
	  for (int j = 0;j < (*spec).lines; j++)//gradient will be Y(p)-Y(p+dp)
	    gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	  
	}
    
    
    if((*spec).params.chap==8 || (*spec).params.chap==9 || (*spec).params.chap==3)  //set gradient of two parameters to zero for gridsearch
      {
	for (int j = 0;j < (*spec).lines; j++)
	  gsl_matrix_set (J, j, 0, 0.0 );
	for (int j = 0;j < (*spec).lines; j++)
	  gsl_matrix_set (J, j, 1, 0.0 );
      }
    
    
    //now check the gradients for each of the additional model parameters required by the minimisation
    if((*spec).params.adduction_flg==1)
      {
	for(int k=0;k<(*spec).params.distparam;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	spec->massy.adduction=(*spec).massy.adduction+STEP;
	double chi2_1=run_calc(0,"Null",0,par,*spec);
	
	for (int j = 0;j < (*spec).lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	i++;
	spec->massy.adduction=(*spec).massy.adduction-STEP;
      }
    
    if((*spec).params.Zfudge_flg==1)
      {
	for(int k=0;k<(*spec).params.distparam;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	spec->massy.Zfudge=(*spec).massy.Zfudge+STEP;
	double chi2_1=run_calc(0,"Null",0,par,*spec);
	for (int j = 0;j < (*spec).lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	i++;
	spec->massy.Zfudge=(*spec).massy.Zfudge-STEP;
      }
    if((*spec).params.MRes_flg==1)
      {
	for(int k=0;k<(*spec).params.distparam;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	spec->massy.MRes=(*spec).massy.MRes+STEP;
	double chi2_1=run_calc(0,"Null",0,par,*spec);
	for (int j = 0;j < (*spec).lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	i++;
	spec->massy.MRes=(*spec).massy.MRes-STEP;
      }
    if((*spec).params.ResFudge_flg==1)
      {
	for(int k=0;k<(*spec).params.distparam;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	spec->massy.ResFudge=(*spec).massy.ResFudge+STEP;
	double chi2_1=run_calc(0,"Null",0,par,*spec);
	for (int j = 0;j < (*spec).lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	i++;
	spec->massy.ResFudge=(*spec).massy.ResFudge-STEP;
      }
    if((*spec).params.tw_flg==1)
      {
	for(int k=0;k<(*spec).params.distparam;k++)
	  par[k]=gsl_vector_get (x, k);//reset distribution parameters
	spec->massy.tw=(*spec).massy.tw+STEP;
	double chi2_1=run_calc(0,"Null",0,par,*spec);
	for (int j = 0;j < (*spec).lines; j++)
	  //gradient will be Y(p)-Y(p+dp)
	  gsl_matrix_set (J, j, i, -1.0*( (*spec).datum[j+(*spec).lines*2]-datum0[j+(*spec).lines*2] )/(STEP*10));
	i++;
	spec->massy.tw=(*spec).massy.tw-STEP;
      }
  }
  
  return GSL_SUCCESS;
}


//CALLS BOTH FUNCTION AND DERIVATIVE CALCULATOR
int
expb_fdf (const gsl_vector * x, void *data,
		gsl_vector * f, gsl_matrix * J)
{
	expb_f (x, data, f);
	expb_df (x,data, J);

	return GSL_SUCCESS;
}


//PRINTS CURRENT STATE OF ITERATOR
void
print_state (size_t iter, gsl_multifit_fdfsolver * s,int p)
{
	printf ("iter: %3u x = ",iter);
	if(p<10)
		for(int i=0;i<p;i++)
			printf("% 15.8f ",gsl_vector_get(s->x,i));
	printf(" |f(x)| = %g\n", gsl_blas_dnrm2 (s->f));

}


//MAIN FUNCTION. SETS UP AND CALLS MINIMSER
double fitty (string outfile,struct data& spec)
{
  int pp=spec.params.distparam;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = spec.lines;
  size_t p = (pp+spec.params.adduction_flg+spec.params.Zfudge_flg+spec.params.MRes_flg+spec.params.ResFudge_flg+spec.params.tw_flg); //add on spectral params
  if(pp>0 && spec.params.chap>=10 && spec.params.tw_flg==1) //for completely free distribution parameters
    p=p-1;


  
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  double x_init2[p];

	{//begin initiation loop...

		int i=0;//parameter increment counter

		if(spec.params.chap==5 || spec.params.chap==6) //each slice has its own gaussian (5) or skewed gaussian (6)
		{

			if(spec.runno==0)//if we're on the first run...
			{
				for(int j=0;j<spec.massy.limSub;j++)
				{

					x_init2[i]=spec.x_init[i];i++;
					x_init2[i]=spec.x_init[i];i++;
					x_init2[i]=spec.x_init[i];i++;
					x_init2[i]=spec.x_init[i];i++;
					printf ("%i Hx0      = %.5f \n",j, x_init2[i-4]);
					printf ("%i Hsig     = %.5f \n",j, x_init2[i-3]);
					printf ("%i alpha    = %.5f \n",j, x_init2[i-2]);
					printf ("%i rel      = %.5f \n",j, x_init2[i-1]);
				}
			}
			else  //otherwise...
			{

				double input[spec.massy.limHSP*spec.massy.limSub];init_array(input,spec.massy.limHSP*spec.massy.limSub);
				make_input(input,spec.massy.limHSP,spec.massy.limSub,spec.params.chap-5,spec.x_init,spec.massy.tw);   //calculate the 2d skewed gaussian distribution
				//	prin_input("testy1.inp","",0,input,spec.limHSP,spec.limSub);    //print input matrix

				double Hx0 =spec.x_init[0]; //extract distribution parameters, HSP centre
				double Hsig=spec.x_init[1]; //extract distribution parameters, HSP width
				double Sx0 =spec.x_init[2]; //extract distribution parameters, client centre
				double Ssig=spec.x_init[3]; //extract distribution parameters, client width
				double alpha=0;
				if(spec.params.chap==6)
					alpha=spec.x_init[4]; //extract distribution parameters

				//double chi2=run_calc(0,0,"Null",0,spec.limHSP,spec.limSub,spec.x_init,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes,spec.Zfudge,spec.ResFudge);

				int ii=int(Hx0);//take nearest integer to Sx0
				cout << " ii is " << ii << endl;
				//set initial parameters for fitting
				for(int j=0;j<spec.massy.limSub;j++)
				{

					if(spec.params.chap==5){
						x_init2[i]=Hx0;i++;
						x_init2[i]=Hsig;i++;
						x_init2[i]=input[ii+spec.massy.limHSP*j];i++;
						printf ("%i Hx0_0      = %.5f \n",j, x_init2[0+3*j]);
						printf ("%i Hsig_0     = %.5f \n",j, x_init2[1+3*j]);
						printf ("%i rel        = %.5f \n",j, x_init2[2+3*j]);}
					else{
						x_init2[i]=Hx0;i++;
						x_init2[i]=Hsig;i++;
						x_init2[i]=alpha;i++;
						x_init2[i]=input[ii+spec.massy.limHSP*j];i++;
						printf ("%i Hx0_0     = %.5f \n", j,x_init2[0+4*j]);
						printf ("%i Hsig_0    = %.5f \n", j,x_init2[1+4*j]);
						printf ("%i alpha     = %.5f \n", j,x_init2[2+4*j]);
						printf ("%i rel       = %.5f \n", j,x_init2[3+4*j]);}
				}
			}

			//init_array(input,spec.limHSP*spec.limSub);
			//make_input(input,spec.limHSP,spec.limSub,5,x_init2,spec.tw);   //calculate the 2d skewed gaussian distributio
			//prin_input("testy2.inp","",0,input,spec.limHSP,spec.limSub);    //print input matrix
			//double chi22=run_calc(0,5,"Null",0,spec.limHSP,spec.limSub,x_init2,spec.adduction,spec.Zwidth,spec.shspMass,spec.subMass,spec.testmax,spec.minZ,spec.maxZ,spec.thresh,spec.datum,spec.lines,spec.MRes,spec.Zfudge,spec.ResFudge);
			//cout << chi22 << endl;

		}


		if(spec.params.chap==4 ) //if 2D skewed Gaussian
		{
			if(spec.params.adduction_flg==0) //go here if coming straight from skewed gaussian (chap=1)
			{

				if(pp>0){
				  
				  if(spec.verb==1)
				    cout << "Fitting " << pp << " distribution parameters" << endl;

				for(i=0;i<5;i++){  //add Hx0, Hsig, Sxo, Ssig and alpha
					x_init2[i]=spec.x_init[i];}}
				x_init2[i]=spec.x_init[0];i++;  //add Hxo_0
				x_init2[i]=spec.x_init[1];i++;  //add Hsig_0
				x_init2[i]=spec.x_init[2];i++;  //add alpha_0

				double input[spec.massy.limHSP*spec.massy.limSub];init_array(input,spec.massy.limHSP*spec.massy.limSub);  //Initialise input matrix
				make_input(input,spec.massy.limHSP,spec.massy.limSub,1,x_init2,spec.massy.tw);                               //Generate input matrix from skewed gaussian
				int imax=findmax(input,spec.massy.limHSP,0);                                            //find maximum value of 0 client
				x_init2[i]=input[imax];i++;                                                       //set this for starting value for rel



			}
			else  //go down this route if already fitted to this distribution once.
			{
				if(pp>0){
				  if(spec.verb==1)
				    cout << "Fitting " << pp << " distribution parameters" << endl;
				for(i=0;i<pp;i++){  //add Hx0, Hsig, Sxo, Ssig and alpha
					x_init2[i]=spec.x_init[i];}}



			}

		}

		if(spec.params.chap<10 && spec.params.chap!=6 && spec.params.chap!=5 && spec.params.chap!=4)  //if chap is less than 10 and not 4,5 or 6 (exceptions above)
		{//setup initial parameter matrix
			if(pp>0){
			  if(spec.verb==1)
			    cout << "Fitting " << pp << " distribution parameters" << endl;
			  if(spec.params.chap!=7)
			    for(i=0;i<pp;i++){
			      x_init2[i]=spec.x_init[i];}
			  else
			    for(i=0;i<pp;i++){
			      x_init2[i]=spec.x_init[2*i+3];
			      cout << "Initialising oligo" << spec.x_init[2*i+2] << " with concentration " << x_init2[i] << endl;
			    }
			}
			
		}


		if(spec.params.chap>=10) //for completely free distribution parameters
		{
		  if(spec.verb==1)
		    {
		    cout << "Running uber fit with " << pp << " distribution parameters" << endl;
		    cout << "Running uber fit with " << p << " total parameters" << endl;
		    }

		  double input[spec.massy.limHSP*spec.massy.limSub];init_array(input,spec.massy.limHSP*spec.massy.limSub);
		  //	  make_input(input,spec.massy.limHSP,spec.massy.limSub,spec.params.chap-5,spec.x_init,spec.massy.tw); //re-calculate the distribution for the last distribution

		  //HERE IS THE CHANGE! READING IN INPUT HERE!
		  //		  make_input(input,spec.massy.limHSP,spec.massy.limSub,spec.params.chap-10,spec.x_init,spec.massy.tw); //re-calculate the distribution for the last distribution
		  read_input(input,spec.massy.limHSP,spec.massy.limSub); //re-calculate the distribution for the last distribution
		  prin_input("testy2.inp","",0,input,spec.massy.limHSP,spec.massy.limSub,x_init2,spec.params.chap);    //print input matrix

			//	for (int j=0;j<5;j++)
			//  cout << spec.x_init[j] << endl;
			for(i=0;i<pp;i++)
				x_init2[i]=input[i];

			if(spec.params.tw_flg==1)
			{
				spec.params.tw_flg=0;  //turn of the 12mer flag at this point, if it has been used

			}
		}


		//the extra model flags are stored at the end of the initation array

		//place in initial guesses for the mass spec parameters
		if(spec.params.adduction_flg==1)
		{
		  if(spec.verb==1)
		    cout << "Fitting adduction" << spec.massy.adduction <<endl;
			x_init2[i]=spec.massy.adduction;i++;
		}
		if(spec.params.Zfudge_flg==1)
		{
		  if(spec.verb==1)
		    cout << "Fitting Zfudge" << spec.massy.Zfudge <<endl;
			x_init2[i]=spec.massy.Zfudge;i++;
		}
		if(spec.params.MRes_flg==1)
		{
		  if(spec.verb==1)
		    cout << "Fitting MRes "  << spec.massy.MRes << endl;
			x_init2[i]=spec.massy.MRes;i++;
		}
		if(spec.params.ResFudge_flg==1)
		{
		  if(spec.verb==1)
		    cout << "Fitting ResFudge "  << spec.massy.ResFudge <<endl;
			x_init2[i]=spec.massy.ResFudge;i++;
		}
		if(spec.params.tw_flg==1)
		  {
		    if(spec.verb==1)
		      cout << "Fitting residual 12mer "  << spec.massy.tw <<endl;

			if(spec.massy.tw>0)
				x_init2[i]=spec.massy.tw;
			else
				x_init2[i]=0.1;
			i++;
		}

	}//end initiation loop


	gsl_vector_view x = gsl_vector_view_array (x_init2, p);

	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &spec;
	if(spec.verb==1)
	  cout << "Fitting with "<< p << " total parameters " << endl;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	if(spec.verb==1)
	  print_state (iter, s,p);
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		//printf ("status = %s\n", gsl_strerror (status));
		
		if(spec.verb==1)
		  print_state (iter, s,p);
		if (status)
			break;
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-8, 1e-8);
	}
	while (status == GSL_CONTINUE && iter < 50);

	//gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))


	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));

	if(spec.verb==1)
	  printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);


	if(spec.params.chap<10 && spec.params.chap!=6 && spec.params.chap!=5)
	{
	  if(spec.verb==1){
		printf("\nCovariance matrix\n");
		for (int i=0;i<pp;i++)
		{
			for(int j=0;j<pp;j++)
				printf("%f\t",c*gsl_matrix_get(covar,i,j));
			printf("\n");
		}}
	}



	double par[spec.massy.limHSP*spec.massy.limSub];

	{//a loop to store ii.
		int ii=0;

		if(spec.params.chap<10 && pp>0 && spec.params.chap!=5 && spec.params.chap!=6 && spec.params.chap!=7 && spec.params.chap!=2)
		{
		  if(spec.verb==1)printf ("Hx0      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  if(spec.verb==1)printf ("Hsig     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  if(spec.verb==1)printf ("Sx0      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  if(spec.verb==1)printf ("Ssig     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;

			if(spec.params.chap==1 || spec.params.chap==9){
			  if(spec.verb==1)printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
			if(spec.params.chap==3){
			  if(spec.verb==1)printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("alpha2    = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
			if(spec.params.chap==4){
			  if(spec.verb==1)printf ("alpha     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("Hx0_0     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("Hsig_0    = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("alpha_0   = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("rel       = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}

			for (int j=0;j<pp;j++)
				par[j]=spec.x_init[j];

		}


		if(spec.params.chap==7)//Modified Poisson for single distribution
		{
		  //printf ("Rat      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  //	printf ("Mult     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  par[0]=spec.x_init[0];
		  par[1]=spec.x_init[1];
		  for (int j=0;j<pp;j++){
		    par[2*j+2]=spec.x_init[2*j+2];
		    par[2*j+3]=FIT(ii);ii++;

		  }
		}


		if(spec.params.chap==2)//Poisson for single distribution
		{
		  if(spec.verb==1)printf ("Rat      = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  if(spec.verb==1)printf ("Mult     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
		  for (int j=0;j<pp;j++)
		    par[j]=spec.x_init[j];

		}




		if(spec.params.chap==5)
		{
			for(int j=0;j<spec.massy.limSub;j++){
			  if(spec.verb==1)printf ("%i Hx0_0     = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("%i Hsig_0    = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("%i alpha_0   = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
			for (int j=0;j<pp;j++)
				par[j]=spec.x_init[j];
		}



		if(spec.params.chap==6)
		{
			for(int j=0;j<spec.massy.limSub;j++){
			  if(spec.verb==1)printf ("%i Hx0_0     = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("%i Hsig_0    = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("%i alpha_0   = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;
			  if(spec.verb==1)printf ("%i rel       = %.5f +/- %.5f\n", j,FIT(ii), c*ERR(ii));spec.x_init[ii]=FIT(ii);ii++;}
			for (int j=0;j<pp;j++)
				par[j]=spec.x_init[j];
		}






		//    if(spec.params.chap<10 && pp==0)//do this if not fitting distribution at all
		//	for (int j=0;j<pp;j++)
		//	  par[j]=spec.x_init[j];


		if(spec.params.chap>=10 && pp>0)//do this if doing the 'all param' dist fit
			for (ii=0;ii<pp;ii++)
				par[ii]=FIT(ii);



		//add on the 'extra' parameters
		if(spec.params.adduction_flg==1){ if(spec.verb==1)printf ("Adduction= %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));   spec.massy.adduction=FIT(ii);ii++;}
		if(spec.params.Zfudge_flg==1){    if(spec.verb==1)printf ("Zfudge   = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));   spec.massy.Zfudge=FIT(ii);   ii++;}
		if(spec.params.MRes_flg==1){  	  if(spec.verb==1)printf ("Mres     = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));   spec.massy.MRes=FIT(ii);     ii++;}
		if(spec.params.ResFudge_flg==1){  if(spec.verb==1)printf ("ResFudge = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));   spec.massy.ResFudge=FIT(ii); ii++;}
		if(spec.params.tw_flg==1){	  if(spec.verb==1)printf ("resid 12 = %.5f +/- %.5f\n", FIT(ii), c*ERR(ii));   spec.massy.tw=FIT(ii);       ii++;}

	}//end of loop for ii


	// printf("chap = %i\t%f\t%f\n",spec.params.chap,par[0],par[4]);

	double chi2=0;
	if(pp==0)
	  for (int i =0; i<10; i++)
	    par[i]=spec.x_init[i];

	if(outfile=="Null")
	  chi2=run_calc(0,outfile,6,par,spec);
	else
	  chi2=run_calc(1,outfile,6,par,spec);

	if(spec.verb==1){
	  printf ("internal chi2/dof %f\n",chi2);
	  printf ("status = %s\n", gsl_strerror (status));}


	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	return chi2;
}



double runchap(string mode,int runno,int lines,double *data,string outfile,struct control& params,double* x_init,struct mass& massy,int verb)
{


//Do this only if the outfile isn't set to 'Null'
  if(outfile!="Null")
  {

    //append tag depending on run number
    if(runno==0){
      char *lab="_a";
      outfile.append(lab,2);}
    if(runno==1){
      char *lab="_b";
      outfile.append(lab,2);}
    if(runno==2){
      char *lab="_c";
      outfile.append(lab,2);}
    if(runno==3){
      char *lab="_d";
      outfile.append(lab,2);}
    if(runno==4){
      char *lab="_e";
      outfile.append(lab,2);}
    if(runno==5){
      char *lab="_f";
      outfile.append(lab,2);}
    if(runno==6){
      char *lab="_g";
      outfile.append(lab,2);}
  }

  //get number of distribution parameters to worry about
  if(params.distparam!=0)
  {
    if(params.chap==0) //1d Gaussian
      params.distparam=4;
    if(params.chap==1)
      params.distparam=5;
    if(params.chap==2) //1d alphaB model (2param)
      params.distparam=2;
    if(params.chap==3)
      params.distparam=6;
    if(params.chap==4)
      params.distparam=9;
    if(params.chap==5)
      params.distparam=massy.limSub*3;
    if(params.chap==6)
      params.distparam=massy.limSub*4;
    if(params.chap==7) //list of specific oligomers
      params.distparam=x_init[1];
    if(params.chap==8)
      params.distparam=4;
    if(params.chap==9)
      params.distparam=5;
    if(params.chap>=10)
      params.distparam=massy.limSub*massy.limHSP;
  }

  
  if(params.chap==7)   //Poisson modified for single dist
    params.distparam=x_init[1];

  if(runno==1 && params.chap==7 && x_init[1]!=0.0)
    //    if(params.distparam>0 && params.chap==7){
    {
      cout << "Number of distribution parameters to fit: " << params.distparam << endl;
      cout << "Principle oligo:     " << x_init[0] << endl;
      cout << "Extra oligos:        " << x_init[1] << endl;
      cout << "Oligo type:          " << x_init[2] << endl;
      cout << "Oligo conc:          " << x_init[3] << endl;
    }


  struct data spec={runno,lines,x_init,data,massy,params,verb};
  
  double chi2;
  if(mode=="fitty")
    chi2=fitty(outfile,spec);
  if(mode=="sim")
    chi2=run_calc(1,outfile,6,data,spec);


  //update parameters in massy for subsequent runs
  massy.adduction=spec.massy.adduction;
  massy.Zfudge=spec.massy.Zfudge;
  massy.MRes=spec.massy.MRes;
  massy.ResFudge=spec.massy.ResFudge;
  massy.tw=spec.massy.tw;

  return chi2;
}


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
  
  if(params.chap==0)
    params.chap=8;
  if(params.chap==2)
    params.chap=3;
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

