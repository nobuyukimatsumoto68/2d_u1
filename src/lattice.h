class Lattice {
public:
  const uint dim;
  const std::vector<uint> size;
  const Idx vol;

  uint calculate_volume() const& {
    Idx res = 1;
    for(const uint elem : size) res *= elem;
    return res;
  }

  Lattice
  (
   const uint dim_,
   const std::vector<uint>& size_
   )
    : dim(dim_)
    , size(size_)
    , vol(calculate_volume())
  {
    assert( size.size() == dim );
  }

  std::string info() const& {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    ss << "--- lattice info ---" << std::endl;
    ss << "dimension = " << dim << std::endl
       << "size = ";
    ss << print(size) << std::endl;

    ss << "volume = " << vol << std::endl;
    ss << "--------------------";

    return ss.str();
  }

  uint mod(const int n, const uint m) const& {
    const int tmp = n + 10*m;
    assert(tmp>0);
    return tmp%m;
  }

  uint fund(const uint mu, const int n) const& {
    return mod(n, size[mu]);
  }

  Idx co2gi( const std::vector<int>& co ) const& {
    assert(co.size()==dim);
    Idx res = 0;
    for(uint mu=0; mu<dim; ++mu){
      const uint tmp = fund(mu, co[mu]);
      res = res*size[mu] + tmp;
    }
    return res;
  }

  Idx co2gi( const std::initializer_list<int>& co ) const& {
    assert(co.size()==dim);
    Idx res = 0;
    auto itr = co.begin();
    for(uint mu=0; mu<dim; ++mu){
      const uint tmp = fund(mu, *itr);
      res = res*size[mu] + tmp;
      ++itr;
    }
    return res;
  }

  std::vector<uint> gi2co(const Idx gi) const& {
    std::vector<uint> v1(dim);
    Idx tmp = gi;
    Idx res = 0;
    for(int mu=dim-1; mu>=0; --mu){
      res = tmp%size[mu];
      tmp -= res;
      tmp /= size[mu];
      v1[mu] = res;
    }
    return v1;
  }


};
