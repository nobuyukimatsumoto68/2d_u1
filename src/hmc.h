class HMC {
public:
  const Lattice& lat;
  const WilsonAction& Sw;
  Rnd& rnd_phi;
  const double eps;
  const int nsteps;
  const int seed;
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

  HMC
  (
   const Lattice& lat_,
   const WilsonAction& Sw_,
   Rnd& rnd_phi_,
   const double eps_,
   const int nsteps_,
   const int seed_,
   const bool is_random_ = false
   )
    : lat(lat_)
    , Sw(Sw_)
    , rnd_phi(rnd_phi_)
    , eps(eps_)
    , nsteps(nsteps_)
    , seed(get_seed(seed_,is_random_))
    , rand01(0.0,1.0)
  {
    assert(&lat==&Sw.lat);
    assert(&lat==&rnd_phi_.lat);

    std::mt19937_64 gen;
    gen.seed(seed);
    rnd.seed( gen() );
  };

  std::string info() const& {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    ss << "--- hmc info ---" << std::endl;
    ss << "eps = " << eps  << std::endl;
    ss << "nsteps = " << nsteps  << std::endl;
    ss << "integrator = " << integrator << std::endl;
    ss << "seed = " << seed << std::endl;
    ss << "--------------------";

    return ss.str();
  }

  void leapfrog
  (
   ScalarField& phi,
   ScalarField& pi,
   const double eps
   ) const& {
    ScalarField dS = Sw.grad(phi);
    pi += dS * (-0.5*eps);

    phi += pi * eps;
    phi.proj_u1();

    dS = Sw.grad(phi);
    pi += dS * (-0.5*eps);
  }

  ScalarField gen_gauss_pi() const& {
    ScalarField pi(rnd_phi.lat, rnd_phi.mult);
    for(Idx gi=0; gi<rnd_phi.size; ++gi) pi[gi] = rnd_phi.gauss(gi);
    return pi;
  }

  double H
  (
   ScalarField& phi,
   ScalarField& pi
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += Sw(phi);
    return res;
  }

  void evolve(ScalarField& phi0,
              bool& is_accept, double& dH) & {
    ScalarField phi(phi0);
    ScalarField pi = gen_gauss_pi();

    const double Hin = H(phi, pi);
    for(int k=0; k<nsteps; ++k) leapfrog(phi,pi,eps);
    const double Hfi = H(phi, pi);

    const double r = rand01(rnd);
    dH = Hfi-Hin;
    const double a = std::min( 1.0, std::exp(-dH) );
    if( r<a ){
      phi0 = phi;
      is_accept=true;
      phi0.proj_u1();
    }
    else is_accept=false;
  }


};
