void relax_modifiedPoisson(MultiGrid & mg_engine,
                           MultiField<Real> * source,
                           MultiField<Real> * field,
                           Real dx2,
                           const Real modif,
                           int level,
                           Real dt = 0)
{
  if(mg_engine.isPartLayer(level))
  {
    if(dt==0)dt = dx2/16.0;
    field[level].updateHalo();
    //double c = 0.25 - modif*dx2*0.125;
    SiteRedBlack3d x(source[level].lattice());
    for(x.first(); x.test(); x.next())
    {
      if(std::isnan(field[level](x) ))
      {
        cout<<"field[level](x) is nan! level: "<< level
            <<" modif*field[level](x)*dt :" << modif*field[level](x)*dt
            <<" source[level](x)*dt :" << source[level](x)*dt
            <<endl;
        parallel.abortForce();
      }
      if(std::isinf(field[level](x) ))
      {
        cout<<"field[level](x) is inf! level: "<< level
            <<" modif*field[level](x)*dt :" << modif*field[level](x)*dt
            <<" source[level](x)*dt :" << source[level](x)*dt
            <<endl;
        parallel.abortForce();
      }

      field[level](x) = field[level](x) +
                      (field[level](x+0)+field[level](x-0)
                      +field[level](x+1)+field[level](x-1)
                      +field[level](x+2)+field[level](x-2)
                      - 6.0*field[level](x) )*dt/dx2
                      - modif*field[level](x)*dt - source[level](x)*dt;

    }
  }
}

void get_residual(MultiGrid & mg_engine,
                  MultiField<Real> * source,
                  MultiField<Real> * field,
                  MultiField<Real> * residual,
                  Real overdx2,
                  const Real modif,
                  int level)
{

  if(mg_engine.isPartLayer(level))
  {
    double maxres = 0;
    double res,modt,src,lap,ph;

    field[level].updateHalo();
    Site x(source[level].lattice());
    Site xmaxlap(source[level].lattice());

    for(x.first();x.test();x.next())
    {
      residual[level](x) = source[level](x)
                           -(field[level](x+0)+field[level](x-0)
                           +field[level](x+1)+field[level](x-1)
                           +field[level](x+2)+field[level](x-2) - field[level](x)*6.0) * overdx2
                           + modif * field[level](x);

      if(abs(residual[level](x))>maxres)
      {
        maxres=abs(residual[level](x));
        lap = (field[level](x+0)+field[level](x-0)
              +field[level](x+1)+field[level](x-1)
              +field[level](x+2)+field[level](x-2)
              - field[level](x)*6.0)  * overdx2;
        src = source[level](x);
        modt = modif * field[level](x);
        res = src - lap + modt;
        ph = field[level](x);
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
             MultiField<Real> * field,
             Real dx2,
             const Real modif,
             double wp)
{

  int level = mg_engine.nl()-1;
  if(mg_engine.isPartLayer(level))
  {
    double res,modt,src,lap,ph;

    double overdx2 = 1.0/dx2;
    Site x(field[level].lattice());
    double error = 1;
    double temp;
    int count = 1;

    //first guess:
    if(modif!=0)
    {
      for(x.first();x.test();x.next())
      {
        field[level](x) = -source[level](x)/modif;
        //cout<<x<<" "<< field[level](x)<<endl;
      }
    }
    else
    {
      for(x.first();x.test();x.next())
      {
        field[level](x) = 0;//source[level](x);
        //cout<<x<<" "<< field[level](x)<<endl;
      }
    }


    if(field[level].lattice().size(0)>2)
    {
      while(error>wp)
      {
        relax_modifiedPoisson(mg_engine,source,field,dx2,modif,level);
        error = 0.0;
        for(x.first();x.test();x.next())
        {
          temp = abs(source[level](x)
                 -(field[level](x+0)+field[level](x-0)
                 +field[level](x+1)+field[level](x-1)
                 +field[level](x+2)+field[level](x-2) - field[level](x)*6.0) * overdx2
                 + modif * field[level](x));

          //if(field[level](x)!=0)temp /= abs(field[level](x));

          if(temp>error)error=temp;

          lap = (field[level](x+0)+field[level](x-0)
                +field[level](x+1)+field[level](x-1)
                +field[level](x+2)+field[level](x-2)
                - field[level](x)*6.0)  * overdx2;
          src = source[level](x);
          modt = modif * field[level](x);
          res = src - lap + modt;
          ph = field[level](x);

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
                                    MultiField<Real> * field,
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
      for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,source,field,dx2[l],modif,l);
      get_residual(mg_engine,source,field,residual,1.0/dx2[l],modif,l);
      mg_engine.restrict(residual,source,l);
      fill0(mg_engine,field,l+1);
    }
    solveCL(mg_engine,source,field,dx2[mg_engine.nl()-1],modif,1.0e-16);
    for(int l = mg_engine.nl()-1 ;l>0;l--)
    {
      mg_engine.prolong(field,residual,l);
      add_residual(mg_engine,field,residual,l-1);
      for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,source,field,dx2[l-1],modif,l-1);
    }
  }
  COUT << "--------- multigrid end ---------"<<endl;
}


void solveModifiedPoisson_linearGammaCycle(MultiGrid & mg_engine,
                                           MultiField<Real> * source,
                                           MultiField<Real> * field,
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

  for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,source,field,dx2[lvl],modif,lvl);
  get_residual(mg_engine,source,field,residual,1.0/dx2[lvl],modif,lvl);
  mg_engine.restrict(residual,source,lvl);
  fill0(mg_engine,field,lvl+1);

  if(lvl+1 == mg_engine.nl()-1)
  {
    //COUT<< "solving level: "<< lvl+1 <<endl;
    solveCL(mg_engine,source,field,dx2[mg_engine.nl()-1],modif,1.0e-15);

  }
  else
  {
    for(int i = 0;i<gamma;i++)
      solveModifiedPoisson_linearGammaCycle(mg_engine,
                                            source,
                                            field,
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
  mg_engine.prolong(field,residual,lvl+1);
  add_residual(mg_engine,field,residual,lvl);
  for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,source,field,dx2[lvl],modif,lvl);
}

void solveModifiedPoisson_linearMGW(MultiGrid & mg_engine,
                                    MultiField<Real> * source,
                                    MultiField<Real> * field,
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
                                          field,
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
                                    MultiField<Real> * field,
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

  solveCL(mg_engine,source,field,dx2[mg_engine.nl()-1],modif,1.0e-12);

  for(j=mg_engine.nl()-2 ; j>=0 ;j--)
  {
    mg_engine.prolong(field,j+1); //prolong field from level j+1 to level j, we will start V cycles from level j
    //COUT<< " starting V cycle from level: "<<j<<endl;
    copy_field(mg_engine,source,rhs,j);
    for(nc = 0; nc<numVcycle; nc++)
    {
      for(l=j ; l < mg_engine.nl()-1 ; l++)
      {
        for(int i=0;i<preRelNum;i++)relax_modifiedPoisson(mg_engine,rhs,field,dx2[l],modif,l);
        get_residual(mg_engine,rhs,field,residual,1.0/dx2[l],modif,l);
        mg_engine.restrict(residual,rhs,l);
        fill0(mg_engine,field,l+1);
      }
      solveCL(mg_engine,rhs,field,dx2[mg_engine.nl()-1],modif,1.0e-11);
      for(l=mg_engine.nl()-1 ; l > j  ; l--)
      {
        mg_engine.prolong(field,residual,l);
        add_residual(mg_engine,field,residual,l-1);
        for(int i=0;i<postRelNum;i++)relax_modifiedPoisson(mg_engine,rhs,field,dx2[l-1],modif,l-1);
      }
    }
  }
  get_residual(mg_engine,source,field,residual,1.0/dx2[0],modif,0);
  COUT << "--------- multigrid end ---------"<<endl;
}


void preparechiSource(Field<Real> &source,
                      Field<Real> &phi,
                      Field<Real> &phiprime,
                      Field<Real> &a4_T0i,
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
    for(int i=0; i<3; i++) temp += a4_T0i(x,i)-a4_T0i(x-i,i);

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

void prepare_chi_extra_source(Field<Real> &source,
                              Field<Real> &phi,
                              Field<Real> &laplace_phi,
                              Field<Real> &Bi,
                              double const dx)
{
  Real temp1, temp2;
  double dx2 = dx*dx;
  int i, j;
  Site x(phi.lattice());
  for(x.first(); x.test(); x.next())
  {
    temp1 = 0.;
    temp2 = 0.;

    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
        if(j-i)
        {
          temp1 += (Bi(x+i+j, i) + Bi(x+j, i)) * (phi(x+i+j+j) - phi(x+i)) // computed in x+i+j
          - (Bi(x+i-j, i) + Bi(x-j, i)) * (phi(x+i) - phi(x+i-j-j)) // computed in x+i-j
          - (Bi(x-i+j, i) + Bi(x-i-i+j, i)) * (phi(x-i+j+j) - phi(x-i)) // computed in x-i+j
          + (Bi(x-i-j, i) + Bi(x-i-j-i, i)) * (phi(x-i) - phi(x-i-j-j)) // computed in x-i-j
          ;
        }
      }
      temp2 += (
               (Bi(x+i, i) + Bi(x, i)) * (phi(x+i+i) - phi(x))
            + (Bi(x-i, i) + Bi(x-i-i, i)) * (phi(x) - phi(x-i-i))
            -2. * ((Bi(x, i) + Bi(x-i, i)) * (phi(x+i) - phi(x-i)))
               ) / dx2;
      temp2 += (Bi(x, i) + Bi(x-i, i)) * (laplace_phi(x+i) - laplace_phi(x-i)) / 4.;
    }
    temp1 /= 4. * dx2;

    temp1 += temp2;
    temp1 /= 2. * dx;
    source(x) += temp1;
  }
}
