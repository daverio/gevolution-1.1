//////////////////////////
// mg_tools.hpp
//////////////////////////
//
// code components related to f(R)
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (CEICO - Czech Academy of Sciences, Prague)
//
//////////////////////////


/////////////////////////////////////////////////
// Check (scalar) field:
// Max = max(fabs(field))
// |hom| = sum(field)/number_of_points
// |avg| average of fabs(field), not fabs(average)!
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
double check_field(Field<FieldType> & field, string field_name, long n3, string message = "") // TODO: correct lattice size
{
  Site x(field.lattice());
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	int prec = 16;
	double max = 0.,
				 hom = 0.,
				 sum = 0.,
				 temp;

  for(x.first(); x.test(); x.next())
  {
    temp = field(x);
    hom += temp;
    sum += fabs(temp);
    if(fabs(temp) >= max)
    {
      max = fabs(temp);
    }
  }
  parallel.max(max);
  parallel.sum(sum);
  parallel.sum(hom);
  sum /= n3;
  hom /= n3;

	COUT << message << scientific << setprecision(prec)
			 << setw(20) << field_name
       << "  Max = " << max;
	if(hom < 0)
	{
		COUT << "  hom = " << hom << endl;
	}
	else
	{
		COUT << "  hom =  " << hom << endl;
	}

	std::cout.copyfmt(oldState);

	return max;
}


/////////////////////////////////////////////////
// Checks vector field, each component separately
// Quantities are the same as check_field
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void check_vector_field(Field<FieldType> & field, string field_name, long n3, string message = "") // TODO: correct lattice size
{
  Site x(field.lattice());
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
	int prec = 16;
  double max[3] = {0.,0.,0.},
				 hom[3] = {0.,0.,0.},
				 temp;

	COUT << message << scientific << setprecision(prec);

  for(int i=0; i<3; i++)
  {
    for(x.first(); x.test(); x.next())
    {
      temp = field(x,i);
      hom[i] += temp;
      if(fabs(temp) > max[i])
      {
        max[i] = fabs(temp);
      }
    }
  parallel.max(max[i]);
	hom[i] /= n3;

  COUT << setw(17) << field_name << "[" << i << "]"
			 << "  Max = " << setw(prec + 3) << max[i];

	if(hom[i] < 0)
	{
		COUT << "  hom = " << hom[i] << endl;
	}
	else
	{
		COUT << "  hom =  " << hom[i] << endl;
	}
  }

	std::cout.copyfmt(oldState);

	return;
}

void check_particles(metadata sim,
										 Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm
										 )
{
	int i, num[8] = {0, 262143, 1, 26, 30000, 180024, 220023, 262142};

	for(i=0; i<2; i++)
	{
		pcls_cdm->coutPart(num[i]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

/////////////////////////////////////////////////
// Flips (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void flip_fields(Field<FieldType> & field1, Field<FieldType> & field2)
{
	Site x(field1.lattice());
	FieldType temp;

	for(x.first(); x.test(); x.next())
	{
		temp = field1(x);
		field1(x) = field2(x);
		field2(x) = temp;
	}
	return;
}

// Writes field1 --> field2, field2 --> field3, field3 --> field1
template <class FieldType>
void flip_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & field3)
{
	Site x(field1.lattice());
	FieldType temp;

	for(x.first(); x.test(); x.next())
	{
		temp = field3(x);
		field3(x) = field2(x);
		field2(x) = field1(x);
		field1(x) = temp;
	}
	return;
}


/////////////////////////////////////////////////
// Sets whole field to zero
/////////////////////////////////////////////////
template <class FieldType>
void zero_field(Field<FieldType> & field)
{
	Site x(field.lattice());
	for(x.first(); x.test(); x.next())
	{
		field(x) = 0.;
	}
}

/////////////////////////////////////////////////
// Copies (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void copy_field(Field<FieldType> & source, Field<FieldType> & destination, double coeff = 1.)
{
	Site x(source.lattice());

	if(coeff == 1.)
	{
		for(x.first(); x.test(); x.next())
		{
			destination(x) = source(x);
		}
	}
	else
	{
		for(x.first(); x.test(); x.next())
		{
			destination(x) = coeff * source(x);
		}
	}
	return;
}

/////////////////////////////////////////////////
// Sums (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void add_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) + field2(x);
	}
	return;
}

template <class FieldType>
void subtract_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) - field2(x);
	}
	return;
}

template <class FieldType>
void add_fields(Field<FieldType> & field1, double coeff1, Field<FieldType> & field2, double coeff2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = coeff1 * field1(x) + coeff2 * field2(x);
	}
	return;
}

template <class FieldType>
void scatter_field(Field<FieldType> & field1, double c1, Field<FieldType> & field2, double c2, Field<FieldType> & result, double k)
{
	Site x(field1.lattice());
	double r;

	for(x.first(); x.test(); x.next())
	{
		r = (double) rand()/RAND_MAX;
		r *= k;
		result(x) = c1 * (1. - r) * field1(x) + c2 * r * field2(x);
	}
	return;
}


// Builds laplacian from field
template <class FieldType>
void build_laplacian(Field<FieldType> & field, Field<FieldType> & laplace, double dx)
{
	double dx2 = dx*dx;
	Site x(field.lattice());

	for(x.first(); x.test(); x.next())
	{
		laplace(x) = field(x+0) + field(x-0) + field(x+1) + field(x-1) + field(x+2) + field(x-2) - 6. * field(x);
		laplace(x) /= dx2;
	}
	return;
}
