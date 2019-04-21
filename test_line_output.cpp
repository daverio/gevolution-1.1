
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"

// f(R) headers
#include "fR_tools.hpp"
#include "background_fR.hpp"

#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "tools.hpp"
#include "output.hpp"
#include "hibernation.hpp"
#include "relaxation.hpp"

using namespace std;
using namespace LATfield2;


int main(int argc, char **argv)
{

    //-------- Initilization of the parallel object ---------
    int n,m,s;
    string filename;

    for (int i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'n':
				n = atoi(argv[++i]);
				break;
			case 'm':
				m =  atoi(argv[++i]);
				break;
      case 's':
  			s =  atoi(argv[++i]);
  			break;
      case 'f':
  			filename = argv[++i]; //output file name
  			break;
		}
	}

    parallel.initialize(n,m);

    Lattice lat(3,s,2);
    Field<Real> rho(lat);

    rho.updateHalo();

    Site x(lat);

    for(x.first();x.test();x.next())
    {
        rho(x)=x.coord(2);
    }

    output_lines(rho, 1, filename);

    //--------------------------------------------------------
}
