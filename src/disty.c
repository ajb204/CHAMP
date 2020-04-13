/*###############################################################*/
//#
//# CHAMP version 7.
//#
//  Chem & Biol. (2012) 19: 599-607  (1D MS citation)
//  Chem. Biol. (2010) 17: 1008-17
//  Proc. Natl. Acad. Sci. U.S.A. (2010), 107: 2007-12
//  J. Mol. Biol. (2011) 413: 297-309
//  J. Mol. Biol. (2011) 413: 310-20
//  Structure (2011) 19: 1855-63   (2D MS citation)
//  J. Am. Chem. Soc. (2012) 134(37) 15343-50
//  Science 359, 6378 930 2018
//
//  Based on code from Matt Bush to calculate MS spectra. Translated
//  into C and wired up to an optimiser. Optimisation of MS spectra is
//  challenging: accomplished here using a combination of steepest
//  descent (fitty) and stimulated annealing (jiggler).
//
//  A.Baldwin (c) University of Oxford 2020.


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
  anal inst;
  inst.ParseInpFile(inputfile);
  for(int i=0;i<inst.files.size();++i)
    {
      inst.files[i].make_summary_init();  //write first line of summary file
      inst.files[i].readfile();  //read datafile
      if(inst.files[i].input_read)
	inst.files[i].read_input();
      inst.files[i].make_input();
      inst.RunProtocol(i);       //execute protocol
      inst.files[i].make_summary_end();
    }
  return 0;
}
