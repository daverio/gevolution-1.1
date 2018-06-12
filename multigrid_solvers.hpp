void relax_modifiedPoisson(MultiGrid & mg_engine,
                           MultiField<Real> * source,
                           MultiField<Real> * phi,
                           Real dx2,
                           const Real modif,
                           int level,
                           Real dt = 0)
{
  if(mg_engine.isPartLayer(level))
  {
    if(dt==0)dt = dx2/16.0;
    phi[level].updateHalo();
    //double c = 0.25 - modif*dx2*0.125;
    SiteRedBlack3d x(source[level].lattice());
    for(x.first();x.test();x.next())
    {
      if(std::isnan(phi[level](x) ))
      {
        cout<<"phi[level](x) is nan! level: "<< level
            <<" modif*phi[level](x)*dt :" << modif*phi[level](x)*dt
            <<" source[level](x)*dt :" << source[level](x)*dt
            <<endl;
        parallel.abortForce();
      }
      if(std::isinf(phi[level](x) ))
      {
        cout<<"phi[level](x) is inf! level: "<< level
            <<" modif*phi[level](x)*dt :" << modif*phi[level](x)*dt
            <<" source[level](x)*dt :" << source[level](x)*dt
            <<endl;
        parallel.abortForce();
      }

      phi[level](x) = phi[level](x) +
                      (phi[level](x+0)+phi[level](x-0)
                      +phi[level](x+1)+phi[level](x-1)
                      +phi[level](x+2)+phi[level](x-2)
                      - 6.0*phi[level](x) )*dt/dx2
                      - modif*phi[level](x)*dt - source[level](x)*dt;

    }
  }
}

void get_residual(MultiGrid & mg_engine,
                  MultiField<Real> * source,
                  MultiField<Real> * phi,
                  MultiField<Real> * residual,
                  Real overdx2,
                  const Real modif,
                  int level)
{

  if(mg_engine.isPartLayer(level))
  {
    double maxres = 0;
    double res,modt,src,lap,ph;

    phi[level].updateHalo();
    Site x(source[level].lattice());
    Site xmaxlap(source[level].lattice());

    for(x.first();x.test();x.next())
    {
      residual[level](x) = source[level](x)
                           -(phi[level](x+0)+phi[level](x-0)
                           +phi[level](x+1)+phi[level](x-1)
                           +phi[level](x+2)+phi[level](x-2) - phi[level](x)*6.0) * overdx2
                           + modif * phi[level](x);

      if(abs(residual[level](x))>maxres)
      {
        maxres=abs(residual[level](x));
        lap = (phi[level](x+0)+phi[level](x-0)
              +phi[level](x+1)+phi[level](x-1)
              +phi[level](x+2)+phi[level](x-2)
              - phi[level](x)*6.0)  * overdx2;
        src = source[level](x);
        modt = modif * phi[level](x);
        res = src - lap + modt;
        ph = phi[level](x);
        xmaxlap = x;
      }
    }
    parallel.layer(mg_engine.player(level)).max(maxres);
    //if(parallel.layer(mg_engine.player(level)).isRoot())
    //    cout<<"get_residual maxres: "<<maxres
    //        << xmaxlap
    //        <<" lap: "<<lap
    //        <<" modt: "<<modt
    //        <<" src: "<<src
    //        <<" res: "<<res
    //        <<" ph: "<<ph
    //        <<endl;
  }
}

void fill0(MultiGrid & mg_engine,
           MultiField<Real> * field,
           int level)
{
  if(mg_engine.isPartLayer(level))
  {
    Site x(field[level].lattice());
    for(x.first();x.test();x.next())field[level](x) = 0.0;
  }
}

void solveCL(MultiGrid & mg_engine,
             MultiField<Real> * source,
             MultiField<Real> * phi,
             Real dx2,
             const Real modif,
             double wp)
{

  int level = mg_engine.nl()-1;
  if(mg_engine.isPartLayer(level))
  {
    double res,modt,src,lap,ph;

    double overdx2 = 1.0/dx2;
    Site x(phi[level].lattice());
    double error = 1;
    double temp;
    int count = 1;

    //first guess:
    if(modif!=0)
    {
      for(x.first();x.test();x.next())
      {
        phi[level](x) = -source[level](x)/modif;
        //cout<<x<<" "<< phi[level](x)<<endl;
      }
    }
    else
    {
      for(x.first();x.test();x.next())
      {
        phi[level](x) = 0;//source[level](x);
        //cout<<x<<" "<< phi[level](x)<<endl;
      }
    }


    if(phi[level].lattice().size(0)>2)
    {
      while(error>wp)
      {
        relax_modifiedPoisson(mg_engine,source,phi,dx2,modif,level);
        error = 0.0;
        for(x.first();x.test();x.next())
        {
          temp = abs(source[level](x)
                 -(phi[level](x+0)+phi[level](x-0)
                 +phi[level](x+1)+phi[level](x-1)
                 +phi[level](x+2)+phi[level](x-2) - phi[level](x)*6.0) * overdx2
                 + modif * phi[level](x));

          //if(phi[level](x)!=0)temp /= abs(phi[level](x));

          if(temp>error)error=temp;

          lap = (phi[level](x+0)+phi[level](x-0)
                +phi[level](x+1)+phi[level](x-1)
                +phi[level](x+2)+phi[level](x-2)
                - phi[level](x)*6.0)  * overdx2;
          src = source[level](x);
          modt = modif * phi[level](x);
          res = src - lap + modt;
          ph = phi[level](x);

        }
        parallel.layer(mg_engine.player(level)).max(error);

        //if(count%1000000==0 && parallel.layer(mg_engine.player(level)).isRoot())
        //    cout<<"SolveCL count: "<<count<< " error: "<<error
        //    <<" lap: "<<lap
        //    <<" modt: "<<modt
        //    <<" src: "<<src
        //    <<" res: "<<res
        //    <<" ph: "<<ph
        //    <<endl;
        //count++;
      }
      //cout<<"SolveCL count: "<<count<< endl;
      //    <<" lap: "<<lap
      //    <<" modt: "<<modt
      //    <<" src: "<<src
      //    <<" res: "<<res
      //    <<" ph: "<<ph
      //    <<endl;
    }


  }
}
void copy_field(MultiGrid & mg_engine,
                MultiField<Real> * source,
                MultiField<Real> *  output,
                int level)
{
  if(mg_engine.isPartLayer(level))
  {
    Site x(source[level].lattice());
    for(x.first();x.test();x.next())
    {
      output[level](x) = source[level](x);
    }
  }
}

void add_residual(MultiGrid & mg_engine,
                  MultiField<Real> * sol,
                  MultiField<Real> * error,
                  int level)
{
  if(mg_engine.isPartLayer(level))
  {
    Site x(sol[level].lattice());
    for(x.first();x.test();x.next())
    {
      sol[level](x) += error[level](x);
    }
  }
}

void solveModifiedPoisson_linearMGV(MultiGrid & mg_engine,
                                    MultiField<Real> * source,
                                    MultiField<Real> * phi,
                                    MultiField<Real> * residual,
                                    MultiField<Real> * rhs,
                                    double dx,
                                    const Real modif = 0.,
                                    const int preRelNum = 1,
                                    const int postRelNum = 1,
                                    const int numVcycle = 1)
{
  double dx2[mg_engine.nl()];
  dx2[0] = dx * dx;
  for(int i=1;i<mg_engine.nl();i++)
  {
    dx2[i]=dx2[i-1]*4.0;
    //COUT<<"level: "<<i<<" ; dx: "<<dx2[i]<<endl;
  }

  //copy_field(mg_engine,source,rhs,0);

  for(int nV = 0;nV< numVcycle;nV++)
  {
    for(int l = 0;l<mg_engine.nl()-1;l++)
    {
      for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,source,phi,dx2[l],modif,l);
      get_residual(mg_engine,source,phi,residual,1.0/dx2[l],modif,l);
      mg_engine.restrict(residual,source,l);
      fill0(mg_engine,phi,l+1);
    }
    solveCL(mg_engine,source,phi,dx2[mg_engine.nl()-1],modif,1.0e-16);
    for(int l = mg_engine.nl()-1 ;l>0;l--)
    {
      mg_engine.prolong(phi,residual,l);
      add_residual(mg_engine,phi,residual,l-1);
      for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,source,phi,dx2[l-1],modif,l-1);
    }
  }
  COUT << "--------- multigrid end ---------"<<endl;
}
void solveModifiedPoisson_linearGammaCycle(MultiGrid & mg_engine,
                                           MultiField<Real> * source,
                                           MultiField<Real> * phi,
                                           MultiField<Real> * residual,
                                           double * dx2,
                                           const Real modif = 0.,
                                           const int preRelNum = 1,
                                           const int postRelNum = 1,
                                           const int gamma = 1,
                                           const int lvl = 0)
{
  //COUT<< "Relaxing level: "<< lvl
  //    << " restricting to level: "<< lvl+1 << endl;

  for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,source,phi,dx2[lvl],modif,lvl);
  get_residual(mg_engine,source,phi,residual,1.0/dx2[lvl],modif,lvl);
  mg_engine.restrict(residual,source,lvl);
  fill0(mg_engine,phi,lvl+1);

  if(lvl+1 == mg_engine.nl()-1)
  {
    //COUT<< "solving level: "<< lvl+1 <<endl;
    solveCL(mg_engine,source,phi,dx2[mg_engine.nl()-1],modif,1.0e-15);

  }
  else
  {
    for(int i = 0;i<gamma;i++)
      solveModifiedPoisson_linearGammaCycle(mg_engine,
                                            source,
                                            phi,
                                            residual,
                                            dx2,
                                            modif,
                                            preRelNum,
                                            postRelNum,
                                            gamma,
                                            lvl+1);
  }
  //COUT<< "Prolonging from level: "<< lvl+1
  //    << " Relaxing level: "<< lvl << endl;
  mg_engine.prolong(phi,residual,lvl+1);
  add_residual(mg_engine,phi,residual,lvl);
  for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,source,phi,dx2[lvl],modif,lvl);
}

void solveModifiedPoisson_linearMGW(MultiGrid & mg_engine,
                                    MultiField<Real> * source,
                                    MultiField<Real> * phi,
                                    MultiField<Real> * residual,
                                    MultiField<Real> * rhs,
                                    double dx,
                                    const Real modif = 0.,
                                    const int preRelNum = 1,
                                    const int postRelNum = 1,
                                    const int numCycle = 1,
                                    const int gamma = 2)
{
  double dx2[mg_engine.nl()];
  dx2[0] = dx * dx;
  for(int i=1;i<mg_engine.nl();i++)
  {
    dx2[i]=dx2[i-1]*4.0;
    //COUT<<"level: "<<i<<" ; dx: "<<dx2[i]<<endl;
  }
  //COUT << "-------- multigrid start --------"<<endl;
  for(int nV = 0;nV< numCycle;nV++)
  {
    solveModifiedPoisson_linearGammaCycle(mg_engine,
                                          source,
                                          phi,
                                          residual,
                                          dx2,
                                          modif,
                                          preRelNum,
                                          postRelNum,
                                          gamma,
                                          0);
  }
  //COUT << "--------- multigrid end ---------"<<endl;
}



void solveModifiedPoisson_linearFMS(MultiGrid & mg_engine,
                                    MultiField<Real> * source,
                                    MultiField<Real> * phi,
                                    MultiField<Real> * residual,
                                    MultiField<Real> * rhs,
                                    Real dx,
                                    const Real modif = 0.,
                                    const int preRelNum = 1,
                                    const int postRelNum = 1,
                                    const int numVcycle = 1)
{

  int j,l,nc;

  double dx2[mg_engine.nl()];
  dx2[0] = dx * dx;
  for(int i=1;i<mg_engine.nl();i++)
  {
    dx2[i]=dx2[i-1]*4.0;
    //COUT<<"level: "<<i<<" ; dx: "<<dx2[i]<<endl;
  }

  for(int i = 0;i<mg_engine.nl()-1;i++)mg_engine.restrict(source,i);

  solveCL(mg_engine,source,phi,dx2[mg_engine.nl()-1],modif,1.0e-12);

  for(j=mg_engine.nl()-2 ; j>=0 ;j--)
  {
    mg_engine.prolong(phi,j+1); //prolong phi from level j+1 to level j, we will start V cycles from level j
    //COUT<< " starting V cycle from level: "<<j<<endl;
    copy_field(mg_engine,source,rhs,j);
    for(nc = 0; nc<numVcycle; nc++)
    {
      for(l=j ; l < mg_engine.nl()-1 ; l++)
      {
        for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,rhs,phi,dx2[l],modif,l);
        get_residual(mg_engine,rhs,phi,residual,1.0/dx2[l],modif,l);
        mg_engine.restrict(residual,rhs,l);
        fill0(mg_engine,phi,l+1);
      }
      solveCL(mg_engine,rhs,phi,dx2[mg_engine.nl()-1],modif,1.0e-11);
      for(l=mg_engine.nl()-1 ; l > j  ; l--)
      {
        mg_engine.prolong(phi,residual,l);
        add_residual(mg_engine,phi,residual,l-1);
        for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,rhs,phi,dx2[l-1],modif,l-1);
      }
    }
  }
  get_residual(mg_engine,source,phi,residual,1.0/dx2[0],modif,0);
  COUT << "--------- multigrid end ---------"<<endl;





}


void preparechiSource(Field<Real> &source,
                      Field<Real> &phi,
                      Field<Real> &phiprime,
                      Field<Real> &T0i,
                      double const H,
                      double const a2,
                      double const dx,
                      double const fourpiG)
{

  Real temp;
  double dx2 = dx*dx;
  Site x(phi.lattice());
  for(x.first();x.test();x.next())
  {



    temp = 0;
    for(int i = 0;i<3;i++) temp += T0i(x,i)-T0i(x-i,i);

    source(x) = fourpiG * temp / a2 / dx;

    source(x) += ( phiprime(x+0) + phiprime(x-0)
                   + phiprime(x+1) + phiprime(x-1)
                   + phiprime(x+2) + phiprime(x-2) - 6.0*phiprime(x) )/dx2;

    source(x) /= H;

    source(x) += ( phi(x+0) + phi(x-0)
                     + phi(x+1) + phi(x-1)
                     + phi(x+2) + phi(x-2) - 6.0*phi(x) )/ dx2;

    //source(x) *= dx2;

  }
}


// Check scalar fields
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
       << "  Max =  " << max;
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
