template <typename T>
std::string print( const std::vector<T>& v ){
  std::stringstream ss;

  bool is_initial = true;
  for(const T elem : v){
    if(!is_initial) ss << ", ";
    ss << elem;
    is_initial = false;
  }

  return ss.str();
}

double proj_mPi_Pi( const double phi ){
  double res = phi;

  while(res>M_PI) res -= 2.0 * M_PI;
  while(res<-M_PI) res += 2.0 * M_PI;

  return res;
}
