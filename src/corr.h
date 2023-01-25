class Corr {
public:
  const Lattice& lat;

  Corr
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

  double w0( const ScalarField& phi ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += std::cos( kappa(phi,x) );
    }

    return res;
  }

  double operator()( const ScalarField& phi, const uint dx ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      Coord y(lat, gi);
      y.shift(0,dx);
      res += std::cos( kappa(phi,x) )*std::cos( kappa(phi,y) );
    }
    res /= lat.vol;

    return res;
  }

  double Q( const ScalarField& phi ) const& {
    double res = 0.0;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      res += proj_mPi_Pi(kappa(phi,x));
    }
    res /= 2.0 * M_PI;

    return res;
  }


};
