class Corr {
public:
  const Lattice& lat;

  Corr
  (
   const Lattice& lat_
   )
    : lat(lat_)
  {}

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

  double operator()( const ScalarField& theta, const uint dx ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      Coord y(lat, gi);
      y.shift(0,dx);
      res += std::cos( kappa(theta,x) )*std::cos( kappa(theta,y) );
    }
    res /= lat.vol;

    return res;
  }

  double Q( const ScalarField& theta ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += proj_mPi_Pi(kappa(theta,x));
    }
    res /= 2.0 * M_PI;

    return res;
  }


};
