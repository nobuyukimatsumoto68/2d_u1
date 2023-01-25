class FieldTrsf {
public:
  const Lattice& lat;
  const Kernel& ker;
  const double eps;
  const int ITER_MAX = 1000;

  FieldTrsf
  (
   const Lattice& lat_,
   const Kernel& ker_,
   const double eps_
   )
    : lat(lat_)
    , ker(ker_)
    , eps(eps_)
  {
    assert(&lat==&ker.lat);
    assert(eps<0.5);
  };

  std::string info() const& {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    ss << "--- field transformation info ---" << std::endl;
    ss << "eps = " << eps << std::endl;
    ss << "--------------------";

    return ss.str();
  }


  ScalarField operator()
  (
   const ScalarField& phi,
   const bool is_even,
   const uint mu
   ) const& {
    ScalarField res(phi);
    const ScalarField dker = ker.grad(phi,is_even,mu);

    res += dker * eps;
    res.proj_u1();

    return res;
  }

  ScalarField inv
  (
   const ScalarField& phi,
   const bool is_even,
   const uint mu
   ) const& {
    ScalarField X(ker.grad(phi, is_even, mu));
    ScalarField res(phi);

    for(int i=0; i<ITER_MAX; ++i){
      res = phi;
      res -= X * eps;

      ScalarField Xold(X);
      X = ker.grad(res, is_even, mu);
      Xold -= X;

      if( std::sqrt( Xold.squaredNorm()/Xold.size )<1.0e-12 ) break;
    }
    res = phi;
    res -= X * eps;
    res.proj_u1();

    return res;
  }


  void star
  (
   ScalarField& dS_Fstar,
   const ScalarField& dS,
   const ScalarField& phi,
   const bool is_even,
   const uint mu
   ) const& {
    dS_Fstar = dS;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Idx gi=0; gi<lat.vol; ++gi) { // y
      const Coord y(lat,gi);
      for(uint nu=0; nu<2; ++nu){
        double sum = 0.0;
        for(int dx0=-1; dx0<=1; ++dx0){
          for(int dx1=-1; dx1<=1; ++dx1){
            Coord x(lat, gi);
            x.shift(0,dx0);
            x.shift(1,dx1);
            if(x.is_even()==is_even) {
              sum += ker.hess(phi,y,nu,x,mu)*dS(x,mu);
            }
          }
        }

        sum *= eps;
        dS_Fstar(gi,nu) += sum;
      }
    }
  }


  double star_det
  (
   const ScalarField& phi,
   const bool is_even,
   const uint mu
   ) const& {
    double res = 1.0;

    for(Idx gi=0; gi<lat.vol; ++gi) { // y
      const Coord x(lat,gi);
      if(x.is_even()==is_even) res *= 1.0+eps*ker.hess(phi,x,mu,x,mu);
    }
    return res;
  }


  double dlndet
  (
   const ScalarField& phi,
   const bool is_even,
   const uint mu,
   const Coord& y, const uint nu
   ) const& {
    double res = 0.0;

    for(int dx0=-1; dx0<=1; ++dx0){
      for(int dx1=-1; dx1<=1; ++dx1){
        Coord x(lat, y());
        x.shift(0,dx0);
        x.shift(1,dx1);
        if(x.is_even()!=is_even) continue;

        double tmp = eps * ker.dd_d(phi,x,mu,y,nu);
        tmp /= 1.0 + eps * ker.hess(phi,x,mu,x,mu);

        res += tmp;
      }
    }

    return res;
  }

  void force
  (
   ScalarField& dSp,
   const ScalarField& dS,
   const ScalarField& phi,
   const bool is_even,
   const uint mu
   ) const& {
    star(dSp, dS, phi, is_even,mu);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Idx gi=0; gi<lat.vol; ++gi) {
      const Coord y(lat,gi);
      for(uint nu=0; nu<2; ++nu){
        const double mdlndet = -dlndet(phi,is_even,mu,y,nu);
        dSp(y,nu) += mdlndet;
      }
    }
  }


};
