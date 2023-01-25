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

  double kappa( const ScalarField& phi, const Coord& x ) const& {
    Coord xp0(x);
    xp0.shift(0);
    Coord xp1(x);
    xp1.shift(1);

    double res = 0.0;
    res += phi(x,0) + phi(xp0,1) - phi(xp1,0) - phi(x,1);
    return res;
  }

  double w0( const ScalarField& phi ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += std::cos( kappa(phi,x) );
    }

    return res;
  }

  double operator()( const ScalarField& phi ) const& {
    return -beta*w0(phi);
  }

  double grad_w0( const ScalarField& phi,
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

  double grad( const ScalarField& phi,
               const Coord& x, const uint mu ) const& {
    return - beta * grad_w0(phi,x,mu);
  }

  ScalarField grad_w0( const ScalarField& phi ) const& {
    ScalarField res(lat,2);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      for(uint mu=0; mu<2; ++mu) res(gi,mu) = grad_w0( phi, x, mu );
    }
    return res;
  }

  ScalarField grad( const ScalarField& phi ) const& { return grad_w0(phi) * (-beta); }

};
