class FT_HMC {
public:
  const Lattice& lat;
  const WilsonAction& Sw;
  const FieldTrsf& trsf;
  const Kernel& kernel;
  Rnd& rnd_phi;
  const double eps_MD;
  const int nsteps_MD;
  const int seed;
  const int nsteps_flow;
  std::uniform_real_distribution<> rand01;
  std::mt19937_64 rnd;
  const std::string integrator = "leapfrog";

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
   const WilsonAction& Sw_,
   const FieldTrsf& trsf_,
   const Kernel& kernel_,
   Rnd& rnd_phi_,
   const double eps_MD_,
   const int nsteps_MD_,
   const int seed_,
   const int nsteps_flow_ = 1,
   const bool is_random_ = false
   )
    : lat(lat_)
    , Sw(Sw_)
    , trsf(trsf_)
    , kernel(kernel_)
    , rnd_phi(rnd_phi_)
    , eps_MD(eps_MD_)
    , nsteps_MD(nsteps_MD_)
    , seed(get_seed(seed_,is_random_))
    , nsteps_flow(nsteps_flow_)
    , rand01(0.0,1.0)
  {
    assert(&lat==&Sw.lat);
    assert(&lat==&rnd_phi_.lat);

    std::mt19937_64 gen;
    gen.seed( seed );
    rnd.seed( gen() );
  };

  std::string info() const& {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    ss << "--- ft-hmc info ---" << std::endl;
    ss << "kernel = " << kernel.description << std::endl;
    ss << "eps_flow = " << trsf.eps << std::endl;
    ss << "nsteps_flow = " << nsteps_flow << std::endl;
    ss << "eps_MD = " << eps_MD  << std::endl;
    ss << "nsteps_MD = " << nsteps_MD  << std::endl;
    ss << "integrator = " << integrator << std::endl;
    ss << "seed = " << seed << std::endl;
    ss << "--------------------";

    return ss.str();
  }

  void leapfrog
  (
   ScalarField& phi_U,
   ScalarField& pi,
   const double eps
   ) const& {
    ScalarField dS_U = Sw.grad(phi_U);
    ScalarField dS_V(dS_U);

    ScalarField phi_V(phi_U);
    for(int i=nsteps_flow-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        phi_V = trsf.inv(phi_V,is_even,mu);
        trsf.force(dS_V, dS_U, phi_V, is_even, mu);
        dS_U = dS_V;

        is_even = true;
        phi_V = trsf.inv(phi_V,is_even,mu);
        trsf.force(dS_V, dS_U, phi_V, is_even, mu);
        dS_U = dS_V;
      }
    }

    pi += dS_V * (-0.5*eps);

    phi_V += pi * eps;
    phi_V.proj_u1();

    phi_U = phi_V;
    for(int i=0; i<nsteps_flow; ++i){
      for(uint mu=0; mu<lat.dim; ++mu) {
        bool is_even = true;
        phi_U = trsf(phi_U,is_even,mu);
        is_even = false;
        phi_U = trsf(phi_U,is_even,mu);
      }
    }

    dS_U = Sw.grad(phi_U);

    phi_V = phi_U;
    for(int i=nsteps_flow-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        phi_V = trsf.inv(phi_V,is_even,mu);
        trsf.force(dS_V,dS_U,phi_V,is_even,mu);
        dS_U = dS_V;

        is_even = true;
        phi_V = trsf.inv(phi_V,is_even,mu);
        trsf.force(dS_V,dS_U,phi_V,is_even,mu);
        dS_U = dS_V;
      }
    }

    pi += dS_V * (-0.5*eps);
  }

  ScalarField gen_gauss_pi() const& {
    ScalarField pi(rnd_phi.lat, rnd_phi.mult);
    for(Idx idx=0; idx<rnd_phi.size; ++idx) pi[idx] = rnd_phi.gauss(idx);
    return pi;
  }

  double H
  (
   const ScalarField& phi_U,
   const ScalarField& phi_V,
   const ScalarField& pi,
   const bool is_even,
   const uint mu
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += Sw(phi_U);
    res -= std::log(trsf.star_det(phi_V,is_even,mu));
    return res;
  }

  double Seff
  (
   const ScalarField& phi_U
   ) const& {
    double res = 0.0;
    res += Sw(phi_U);

    ScalarField phi_V(phi_U);
    for(int i=nsteps_flow-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        phi_V = trsf.inv(phi_V,is_even,mu);
        res -= std::log(trsf.star_det(phi_V,is_even,mu));

        is_even = true;
        phi_V = trsf.inv(phi_V,is_even,mu);
        res -= std::log(trsf.star_det(phi_V,is_even,mu));
      }
    }

    return res;
  }

  double H
  (
   const ScalarField& phi_U,
   const ScalarField& pi
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += Seff( phi_U );
    return res;
  }

  void evolve(ScalarField& phi_U0,
              bool& is_accept, double& dH
              ) & {
    ScalarField phi_U( phi_U0 );
    ScalarField pi = gen_gauss_pi();

    const double Hin = H(phi_U,pi);
    for(int k=0; k<nsteps_MD; ++k) {
      leapfrog(phi_U,pi,eps_MD);
    }
    const double Hfi = H(phi_U,pi);

    const double r = rand01(rnd);
    dH = Hfi-Hin;
    const double a = std::min( 1.0, std::exp(-dH) );
    if( r<a ){
      phi_U0 = phi_U;
      is_accept=true;
      phi_U0.proj_u1();
    }
    else is_accept=false;
  }

};
