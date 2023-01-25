class ScalarField {
  // private:
public:
  const Lattice& lat;
  const uint mult;
  const Idx size;
  std::vector<double> field;

  ScalarField
  (
   const Lattice& lat_,
   const uint mult_
   )
    : lat(lat_)
    , mult(mult_)
    , size(lat.vol*mult)
    , field( std::vector<double>(size, 0.0) )
  {
    assert(field.size()==size);
  }

  ScalarField
  (
   const ScalarField& other
   )
    : lat(other.lat)
    , mult(other.mult)
    , size(lat.vol*mult)
    , field( other.field )
  {
    assert(field.size()==size);
  }

  inline Idx idx(const Idx gi, const uint i) const& { return mult*gi + i; }
  inline Idx idx(const Coord& x, const uint i) const& { return idx(lat.co2gi(x()), i); }

  inline double operator[](const Idx idx) const& {
    assert(idx<size);
    return field[idx];
  }
  inline double& operator[](const Idx idx) & {
    assert(idx<size);
    return field[idx];
  }

  inline double operator()(const Idx gi, const uint i) const& { return (*this)[idx(gi,i)]; }
  inline double& operator()(const Idx gi, const uint i) & { return (*this)[idx(gi,i)]; }

  inline double operator()(const Coord& x, const uint i) const& { return (*this)[idx(x,i)]; }
  inline double& operator()(const Coord& x, const uint i) & { return (*this)[idx(x,i)]; }

  inline double operator()(const std::initializer_list<int>& x, const uint i) const& {
    return (*this)[idx(lat.co2gi(x),i)];
  }
  inline double& operator()(const std::initializer_list<int>& x, const uint i) & {
    return (*this)[idx(lat.co2gi(x),i)];
  }

  ScalarField& operator=( const ScalarField& other ) {
    if(this != &other) {
      assert(&lat==&other.lat && size==other.size);
      field = other.field;
    }
    return *this;
  }

  std::string show() const& {
    std::stringstream ss;

    for(Idx gi=0; gi<lat.vol; ++gi){
      const Coord x(lat, gi);
      ss << "x = " << x.show() << std::endl;

      for(uint i=0; i<mult; ++i){
        if(i!=0) ss << ", ";
        ss << (*this)(gi,i);
      }
      if(gi!=lat.vol-1) ss << std::endl;
    }

    return ss.str();
  }

  void proj_u1(){ // see util.h
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(Idx idx=0; idx<size; ++idx) (*this)[idx] = proj_mPi_Pi( (*this)[idx] );
  }

  ScalarField& operator+=( const ScalarField& other ){
    for(Idx gi=0; gi<size; ++gi) (*this)[gi] += other[gi];
    return *this;
  }

  ScalarField& operator*=( const double coeff){
    for(Idx gi=0; gi<size; ++gi) (*this)[gi] *= coeff;
    return *this;
  }

  ScalarField operator*( const double coeff) const& {
    ScalarField res(*this);
    res *= coeff;
    return res;
  }

  ScalarField& operator-=( const ScalarField& other ){
    *this += other * (-1.0);
    return *this;
  }

  double squaredNorm() const& {
    double res = 0.0;
    for( Idx gi=0; gi<size; ++gi ) res += (*this)[gi] * (*this)[gi];
    return res;
  }

};


