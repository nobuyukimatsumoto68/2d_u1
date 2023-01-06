class WilsonAction {
public:
  const Lattice& lat;
  const double beta;

  WilsonAction
  (
   const Lattice& lat_,
   const double beta_
   )
    : lat(lat_)
    , beta(beta_)
  {}

  std::string info() const& {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    ss << "--- action info ---" << std::endl;
    ss << "Wilson action" << std::endl
       << "beta = " << beta << std::endl;
    ss << "--------------------";

    return ss.str();
  }

  double kappa( const ScalarField& theta, const Coord& x ) const& {
    Coord xp0(x);
    xp0.shift(0,1);
    Coord xp1(x);
    xp1.shift(1,1);

    double res = 0.0;
    res += theta(x,0) + theta(xp0,1) - theta(xp1,0) - theta(x,1);
    return res;
  }

  double w0( const ScalarField& theta ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += std::cos( kappa(theta,x) );
    }

    return res;
  }

  double operator()( const ScalarField& theta ) const& {
    return -beta*w0(theta);
  }

  double Q( const ScalarField& theta ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += proj_mPi_Pi(kappa(theta,x));
    }
    res /= -2.0 * M_PI;

    return res;
  }

  double grad_w0( const ScalarField& theta,
                 const Coord& x, const uint mu ) const& {
    double res = 0.0;
    if(mu==0){
      Coord xm1(x);
      xm1.shift(1,-1);
      res -= std::sin(kappa(theta,x)) - std::sin(kappa(theta,xm1));
    }
    else if(mu==1){
      Coord xm0(x);
      xm0.shift(0,-1);
      res -= std::sin(kappa(theta,xm0)) - std::sin(kappa(theta,x));
    }
    else assert(false);

    return res;
  }

  double grad( const ScalarField& theta,
               const Coord& x, const uint mu ) const& {
    return - beta * grad_w0(theta,x,mu);
  }

  ScalarField grad_w0( const ScalarField& theta ) const& {
    ScalarField res(lat,2);
    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      for(uint mu=0; mu<2; ++mu) res(gi,mu) = grad_w0( theta, x, mu );
    }
    return res;
  }

  ScalarField grad( const ScalarField& theta ) const& { return grad_w0(theta) * (-beta); }

  double hess_w0( const ScalarField& theta,
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
        res -= std::cos(kappa(theta,x)) + std::cos(kappa(theta,xm1));
      }
      else if(mu==1){
        Coord xm0(x);
        xm0.shift(0,-1);
        res -= std::cos(kappa(theta,xm0)) + std::cos(kappa(theta,x));
      }
    }
    else if(y==x) res = std::cos(kappa(theta,x));
    else if(mu==0){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      if(y==xp0 && nu==1) res = -std::cos(kappa(theta,x));
      else if(y==xp1 && nu==0) res = std::cos(kappa(theta,x));
    }
    else if(mu==1){
      Coord xp0(x);
      xp0.shift(0);
      Coord xp1(x);
      xp1.shift(1);
      if(y==xp0 && nu==1) res = std::cos(kappa(theta,x));
      else if(y==xp1 && nu==0) res = -std::cos(kappa(theta,x));
    }

    return res;
  }


};
