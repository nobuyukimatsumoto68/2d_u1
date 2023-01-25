class HMC {
public:
  const Lattice& lat;
  const WilsonAction& S;
  Rnd& rnd_theta;
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
   const WilsonAction& S_,
   Rnd& rnd_theta_,
   const double eps_,
   const int nsteps_,
   const int seed_,
   const bool is_random_ = false
   )
    : lat(lat_)
    , S(S_)
    , rnd_theta(rnd_theta_)
    , eps(eps_)
    , nsteps(nsteps_)
    , seed(get_seed(seed_,is_random_))
    , rand01(0.0,1.0)
  {
    assert(&lat==&S.lat);
    assert(&lat==&rnd_theta_.lat);

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
   ScalarField& theta,
   ScalarField& pi,
   const double eps
   ) const& {
    ScalarField dS = S.grad(theta);
    pi += dS * (-0.5*eps);

    theta += pi * eps;
    theta.proj_u1();

    dS = S.grad(theta);
    pi += dS * (-0.5*eps);
  }

  ScalarField gen_gauss_pi() const& {
    ScalarField pi(rnd_theta.lat, rnd_theta.mult);
    for(Idx gi=0; gi<rnd_theta.size; ++gi) pi[gi] = rnd_theta.gauss(gi);
    return pi;
  }

  double H
  (
   ScalarField& theta,
   ScalarField& pi
   ) const& {
    double res = 0.0;
    res += 0.5 * pi.squaredNorm();
    res += S(theta);
    return res;
  }

  void evolve(ScalarField& theta0,
              bool& is_accept, double& dH) & {
    ScalarField theta(theta0);
    ScalarField pi = gen_gauss_pi();

    const double Hin = H(theta, pi);
    for(int k=0; k<nsteps; ++k) leapfrog(theta,pi,eps);
    const double Hfi = H(theta, pi);

    const double r = rand01(rnd);
    dH = Hfi-Hin;
    const double a = std::min( 1.0, std::exp(-dH) );
    if( r<a ){
      theta0 = theta;
      is_accept=true;
      theta0.proj_u1();
    }
    else is_accept=false;
  }


};
