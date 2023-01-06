class Rnd {
public:
  const Lattice& lat;
  const uint mult;
  const Idx size;
  const int seed;
  std::vector<std::mt19937_64> rnd;
  std::normal_distribution<double> normal;

  int get_seed(const int seed_, const bool is_random){
    int res = seed_;
    if(is_random) {
      std::random_device r;
      res = r();
    }
    return res;
  }

  Rnd
  (
   const Lattice& lat_,
   const uint mult_,
   const int seed_,
   const bool is_random=false
   )
    : lat(lat_)
    , mult(mult_)
    , size(lat.vol*mult)
    , seed(get_seed(seed_,is_random))
    , rnd(size)
  {
    std::mt19937_64 gen;
    gen.seed(seed);
    for(Idx gi=0; gi<size; ++gi) rnd[gi].seed( gen() );
  };

  std::string info() const& {
    std::stringstream ss;

    ss << "--- Rnd info ---" << std::endl;
    ss << "seed = " << seed << std::endl;
    ss << "mult = " << mult << std::endl;
    ss << "--------------------";

    return ss.str();
  }

  inline Idx idx(const Idx gi, const uint i) const& { return mult*gi + i; }
  inline Idx idx(const Coord& x, const uint i) const& { return idx(lat.co2gi(x()), i); }

  double gauss(const Idx idx) & { return normal( rnd[idx] ); }

  void set( ScalarField& theta ) & {
    for(Idx gi=0; gi<size; ++gi) theta[gi] = gauss(gi);
    theta.proj_u1();
  }

};
