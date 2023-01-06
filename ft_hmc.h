class FT_HMC {
public:
  const Lattice& lat;
  const WilsonAction& S;
  const FieldTrsf& trsf;
  const Kernel& kernel;
  Rnd& rnd_theta;
  const int seed;
  const int nsteps;
  std::uniform_real_distribution<> rand01;
  std::mt19937_64 rnd;

  int get_seed(const int seed_, const bool is_random){
    int res = seed_;
    if(is_random) {
      std::random_device r;
      res = r();
    }
    return res;
  }

  FT_HMC
  (
   const Lattice& lat_,
   const WilsonAction& S_,
   const FieldTrsf& trsf_,
   const Kernel& kernel_,
   Rnd& rnd_theta_,
   const int seed_,
   const int nsteps_ = 1,
   const bool is_random_ = false
   )
    : lat(lat_)
    , S(S_)
    , trsf(trsf_)
    , kernel(kernel_)
    , rnd_theta(rnd_theta_)
    , seed(get_seed(seed_,is_random_))
    , nsteps(nsteps_)
    , rand01(0.0,1.0)
  {
    assert(&lat==&S.lat);
    assert(&lat==&rnd_theta_.lat);

    std::mt19937_64 gen;
    gen.seed( seed );
    rnd.seed( gen() );
  };


  void leapfrog
  (
   ScalarField& theta_U,
   ScalarField& pi,
   const double eps
   ) const& {
    ScalarField dS_U = S.grad(theta_U);
    ScalarField dS_V(dS_U);

    ScalarField theta_V(theta_U);
    for(int i=nsteps-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        theta_V = trsf.inv(theta_V,is_even,mu);
        trsf.dS(dS_V, dS_U, theta_V, is_even, mu);
        dS_U = dS_V;

        is_even = true;
        theta_V = trsf.inv(theta_V,is_even,mu);
        trsf.dS(dS_V, dS_U, theta_V, is_even, mu);
        dS_U = dS_V;
      }
    }

    pi += dS_V * (-0.5*eps);

    theta_V += pi * eps;
    theta_V.proj_u1();

    theta_U = theta_V;
    for(int i=0; i<nsteps; ++i){
      for(uint mu=0; mu<lat.dim; ++mu) {
        bool is_even = true;
        theta_U = trsf(theta_U,is_even,mu);
        is_even = false;
        theta_U = trsf(theta_U,is_even,mu);
      }
    }

    dS_U = S.grad(theta_U);

    theta_V = theta_U;
    for(int i=nsteps-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        theta_V = trsf.inv(theta_V,is_even,mu);
        trsf.dS(dS_V,dS_U,theta_V,is_even,mu);
        dS_U = dS_V;

        is_even = true;
        theta_V = trsf.inv(theta_V,is_even,mu);
        trsf.dS(dS_V,dS_U,theta_V,is_even,mu);
        dS_U = dS_V;
      }
    }

    pi += dS_V * (-0.5*eps);
  }

  ScalarField gen_gauss_pi() const& {
    ScalarField pi(rnd_theta.lat, rnd_theta.mult);
    for(Idx idx=0; idx<rnd_theta.size; ++idx) pi[idx] = rnd_theta.gauss(idx);
    return pi;
  }

  double H
  (
   const ScalarField& theta_U,
   const ScalarField& theta_V,
   const ScalarField& pi,
   const bool is_even,
   const uint mu
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += S(theta_U);
    res -= std::log(trsf.star_det(theta_V,is_even,mu));
    return res;
  }

  double Seff
  (
   const ScalarField& theta_U
   ) const& {
    double res = 0.0;
    res += S(theta_U);

    ScalarField theta_V(theta_U);
    for(int i=nsteps-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        theta_V = trsf.inv(theta_V,is_even,mu);
        res -= std::log(trsf.star_det(theta_V,is_even,mu));

        is_even = true;
        theta_V = trsf.inv(theta_V,is_even,mu);
        res -= std::log(trsf.star_det(theta_V,is_even,mu));
      }
    }

    return res;
  }

  double H
  (
   const ScalarField& theta_U,
   const ScalarField& pi
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += Seff( theta_U );
    return res;
  }

  void evolve(ScalarField& theta_U0,
              bool& is_accept, double& dH,
              const double eps_MD, const int nsteps_MD
              ) & {
    ScalarField theta_U( theta_U0 );
    ScalarField pi = gen_gauss_pi();

    const double Hin = H(theta_U,pi);
    for(int k=0; k<nsteps_MD; ++k) {
      leapfrog(theta_U,pi,eps_MD);
    }
    const double Hfi = H(theta_U,pi);

    const double r = rand01(rnd);
    dH = Hfi-Hin;
    const double a = std::min( 1.0, std::exp(-dH) );
    if( r<a ){
      theta_U0 = theta_U;
      is_accept=true;
      theta_U0.proj_u1();
    }
    else is_accept=false;
  }

};
