class Coord {
// private:
public:
  const Lattice& lat;
  std::vector<int> co;

  Coord
  (
   const Lattice& lat_
   )
    : lat(lat_)
    , co(std::vector<int>(lat.dim))
  {
  }

  Coord
  (
   const Lattice& lat_,
   const std::vector<int> co_
   )
    : lat(lat_)
    , co(co_)
  {
  }

  Coord
  (
   const Lattice& lat_,
   const std::initializer_list<int>& co_
   )
    : lat(lat_)
    , co(lat.dim)
  {
    assert(co_.size()==lat.dim);

    auto itr = co_.begin();
    for(uint mu=0; mu<lat.dim; ++mu){
      co[mu] = *itr;
      ++itr;
    }
  }

  Coord
  (
   const Lattice& lat_,
   const Idx gi_
   )
    : lat(lat_)
    , co(std::vector<int>(lat.dim))
  {
    set( gi_ );
  }

  std::vector<int> operator()() const& {
    return co;
  }

  int operator[](const uint mu) const& {
    assert(mu<lat.dim);
    return lat.fund(mu,co[mu]);
  }

  int& operator[](const uint mu) & {
    return co[mu];
  }

  bool operator<(const Coord& y) const& {
    bool res=false;
    Coord xp0(*this);
    xp0.shift(0);
    Coord xp1(*this);
    xp1.shift(1);
    Coord xp0p1(xp0);
    xp0p1.shift(1);
    Coord xp0m1(xp0);
    xp0m1.shift(1,-1);
    if( y==xp0 || y==xp1 || y==xp0p1 || y==xp0m1 ) res = true;

    return res;
  }

  bool operator>(const Coord& y) const& {
    bool res=false;
    Coord xm0(*this);
    xm0.shift(0,-1);
    Coord xm1(*this);
    xm1.shift(1,-1);
    Coord xm0m1(xm0);
    xm0m1.shift(1,-1);
    Coord xm0p1(xm0);
    xm0p1.shift(1);
    if( y==xm0 || y==xm1 || y==xm0m1 || y==xm0p1 ) res = true;

    return res;
  }

  bool operator<=(const Coord& y) const& { return !( (*this)>y ); }
  bool operator>=(const Coord& y) const& { return !( (*this)<y ); }

  bool operator==(const Coord& y) const& {
    bool res=true;
    for(uint mu=0; mu<lat.dim; ++mu) {
      if( (*this)[mu]!=y[mu] ) {
        res = false;
        break;
      }
    }
    return res;
  }

  bool is_even() const& {
    int sum = 0;
    for(uint mu=0; mu<lat.dim; ++mu) sum += (*this)[mu];
    sum = (sum+2024)%2;

    const bool res = sum==0 ? true : false;
    return res;
  }

  void set(const std::vector<int>& co_) & {
    assert(co_.size()==lat.dim);
    for(uint mu=0; mu<lat.dim; ++mu) co[mu] = lat.fund(mu, co_[mu]);
  }

  void set(const std::vector<uint>& co_) & {
    assert(co_.size()==lat.dim);
    for(uint mu=0; mu<lat.dim; ++mu) co[mu] = lat.fund(mu, co_[mu]);
  }

  void set(const Idx gi) & {
    set( lat.gi2co(gi) );
  }

  void reshape() & { set(co); }

  uint size() const& { return lat.dim; }

  void shift(const uint mu, const int dn=1){
    assert(mu<lat.dim);
    co[mu] = lat.fund(mu, co[mu]+dn);
  }

  std::string show() const& {
    std::stringstream ss;
    ss << '(' << print(co) << ") ";
    return ss.str();
  }

};



