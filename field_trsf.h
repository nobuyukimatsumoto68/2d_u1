class FieldTrsf {
public:
  const Lattice& lat;
  const Kernel& Stilde;
  const double eps;
  const int ITER_MAX = 1000;

  FieldTrsf
  (
   const Lattice& lat_,
   const Kernel& Stilde_,
   const double eps_
   )
    : lat(lat_)
    , Stilde(Stilde_)
    , eps(eps_)
  {
    assert(&lat==&Stilde.lat);
    assert(eps<0.5);
  };

  ScalarField operator()
  (
   const ScalarField& theta,
   const bool is_even,
   const uint mu
   ) const& {
    ScalarField res(theta);
    const ScalarField dStilde = Stilde.grad(theta,is_even,mu);

    res += dStilde * eps;
    res.proj_u1();

    return res;
  }

  ScalarField inv
  (
   const ScalarField& theta,
   const bool is_even,
   const uint mu
   ) const& {
    ScalarField X(Stilde.grad(theta, is_even, mu));
    ScalarField res(theta);

    for(int i=0; i<ITER_MAX; ++i){
      res = theta;
      res -= X * eps;

      ScalarField Xold(X);
      X = Stilde.grad(res, is_even, mu);
      Xold -= X;

      if( std::sqrt( Xold.squaredNorm()/Xold.size )<1.0e-12 ) break;
    }
    res = theta;
    res -= X * eps;
    res.proj_u1();

    return res;
  }


  void star
  (
   ScalarField& dS_Fstar,
   const ScalarField& dS,
   const ScalarField& theta,
   const bool is_even,
   const uint mu
   ) const& {
    dS_Fstar = dS;

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
              sum += Stilde.hess(theta,y,nu,x,mu)*dS(x,mu);
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
   const ScalarField& theta,
   const bool is_even,
   const uint mu
   ) const& {
    double res = 1.0;
    for(Idx gi=0; gi<lat.vol; ++gi) { // y
      const Coord x(lat,gi);
      if(x.is_even()==is_even) res *= 1.0+eps*Stilde.hess(theta,x,mu,x,mu);
    }
    return res;
  }


  double dlndet
  (
   const ScalarField& theta,
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

        double tmp = eps * Stilde.dd_d(theta,x,mu,y,nu);
        tmp /= 1.0 + eps * Stilde.hess(theta,x,mu,x,mu);

        res += tmp;
      }
    }

    return res;
  }

  void dS
  (
   ScalarField& dSp,
   const ScalarField& dS,
   const ScalarField& theta,
   const bool is_even,
   const uint mu
   ) const& {
    star(dSp, dS, theta, is_even,mu);
    for(Idx gi=0; gi<lat.vol; ++gi) {
      const Coord y(lat,gi);
      for(uint nu=0; nu<2; ++nu){
        const double mdlndet = -dlndet(theta,is_even,mu,y,nu);
        dSp(y,nu) += mdlndet;
      }
    }
  }


};
