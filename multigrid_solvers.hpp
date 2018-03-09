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
    for(x.first();x.test();x.next())
    {
      phi[level](x) = -source[level](x)/modif;
      //cout<<x<<" "<< phi[level](x)<<endl;
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
