class Kernel {
public:
  const Lattice& lat;
  const std::string description = "wilson flow";

  Kernel
  (
   const Lattice& lat_
   )
    : lat(lat_)
  {}

  double kappa( const ScalarField& phi, const Coord& x ) const& {
    Coord xp0(x);
    xp0.shift(0,1);
    Coord xp1(x);
    xp1.shift(1,1);

    double res = 0.0;
    res += phi(x,0) + phi(xp0,1) - phi(xp1,0) - phi(x,1);
    return res;
  }

  double operator()( const ScalarField& phi ) const& {
    double res = 0.0;
    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += std::cos( kappa(phi,x) );
    }
    return res;
  }

  double grad( const ScalarField& phi,
               const Coord& x, const uint mu ) const& {
    double res = 0.0;
    if(mu==0){
      Coord xm1(x);
      xm1.shift(1,-1);
      res -= std::sin(kappa(phi,x)) - std::sin(kappa(phi,xm1));
    }
    else if(mu==1){
      Coord xm0(x);
      xm0.shift(0,-1);
      res -= std::sin(kappa(phi,xm0)) - std::sin(kappa(phi,x));
    }
    else assert(false);

    return res;
  }

  ScalarField grad( const ScalarField& phi,
                    const bool is_even,
                    const uint mu ) const& {
    ScalarField res(lat,2);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      if(x.is_even()==is_even) res(gi,mu) = grad( phi, x, mu );
    }
    return res;
  }


  double hess( const ScalarField& phi,
               const Coord& x0, const uint mu0,
               const Coord& y0, const uint nu0 ) const& {
    double res = 0.0;

    Coord x(lat), y(lat);
    uint mu, nu;
    if(x0<=y0){
      x.set( x0() );
      y.set( y0() );
      mu=mu0;
      nu=nu0;
    }
    else{
      x.set( y0() );
      y.set( x0() );
      mu=nu0;
      nu=mu0;
    }

    if(mu==nu && y==x){
      if(mu==0){
        Coord xm1(x);
        xm1.shift(1,-1);
        res -= std::cos(kappa(phi,x)) + std::cos(kappa(phi,xm1));
      }
      else if(mu==1){
        Coord xm0(x);
        xm0.shift(0,-1);
        res -= std::cos(kappa(phi,xm0)) + std::cos(kappa(phi,x));
      }
    }
    else if(y==x) res = std::cos(kappa(phi,x));
    else if(mu==0){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      Coord xp0m1(xp0);
      xp0m1.shift(1,-1);
      if(y==xp0 && nu==1) res = -std::cos(kappa(phi,x));
      else if(y==xp1 && nu==0) res = std::cos(kappa(phi,x));
      else if(y==xp0m1 && nu==1) {
        Coord xm1(x);
        xm1.shift(1,-1);
        res = std::cos(kappa(phi,xm1));
      }
    }
    else if(mu==1){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      if(y==xp0 && nu==1) res = std::cos(kappa(phi,x));
      else if(y==xp1 && nu==0) res = -std::cos(kappa(phi,x));
    }

    return res;
  }

  double dd_d( const ScalarField& phi,
               const Coord& x0, const uint mu0,
               const Coord& y0, const uint nu0 ) const& {
    double res = 0.0;
    double sign = 1.0;

    Coord x(lat), y(lat);
    uint mu, nu;
    bool is_swapped = false;
    if(x0<=y0){
      x.set( x0() );
      y.set( y0() );
      mu=mu0;
      nu=nu0;
    }
    else{
      x.set( y0() );
      y.set( x0() );
      mu=nu0;
      nu=mu0;
      is_swapped = true;
    }

    if(mu==nu && y==x){
      if(mu==0){
        Coord xm1(x);
        xm1.shift(1,-1);
        res += std::sin(kappa(phi,x)) - std::sin(kappa(phi,xm1));
      }
      else if(mu==1){
        Coord xm0(x);
        xm0.shift(0,-1);
        res += std::sin(kappa(phi,xm0)) - std::sin(kappa(phi,x));
      }
    }
    else if(y==x) {
      if(mu==1) sign = -1.0;
      res = -std::sin(kappa(phi,x));
    }
    else if(mu==0){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      Coord xp0m1(xp0);
      xp0m1.shift(1,-1);
      if(y==xp0 && nu==1) res = std::sin(kappa(phi,x));
      else if(y==xp1 && nu==0) {
        if(is_swapped) sign = -1.0;
        res = -std::sin(kappa(phi,x));
      }
      else if(y==xp0m1 && nu==1) {
        Coord xm1(x);
        xm1.shift(1,-1);
        if(is_swapped) sign = -1.0;
        res = std::sin(kappa(phi,xm1));
      }
    }
    else if(mu==1){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      if(y==xp0 && nu==1) {
        if(is_swapped) sign = -1.0;
        res = std::sin(kappa(phi,x));
      }
      else if(y==xp1 && nu==0) res = -std::sin(kappa(phi,x));
    }

    return sign*res;
  }


};
