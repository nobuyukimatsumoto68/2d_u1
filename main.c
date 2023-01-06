#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <complex>
#include <array>
#include <cassert>
#include <vector>
#include <string>
#include <filesystem>

using Complex = std::complex<double>;
constexpr Complex I = Complex(0,1);
using uint = unsigned int;
using Idx = std::size_t;

#include "util.h"
#include "lattice.h"
#include "coord.h"
#include "scalar_field.h"
#include "action.h"
#include "rnd.h"
#include "hmc.h"
#include "kernel.h"
#include "field_trsf.h"
#include "ft_hmc.h"

// ---
const uint dimension = 2;
// ---
const double a = 1.0/1.6;
const double beta = 1.0/(a*a);
const std::vector<uint> lattice_size = std::vector<uint>{11,11};
// ---
const std::string dir = "./results/";
// ----
const double stot = 1.0;
const int nsteps = 10;
const double eps = stot/nsteps;
const int ninit = 20;
const int nconf = 20;
// ----
const double eps_W = 0.1;
const int nsteps_W = 2;
// ----
const int seed_1 = 1;
const int seed_2 = 2;


int main(){

  std::cout << std::scientific << std::setprecision(15);
  std::cerr << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);

  std::filesystem::create_directory(dir);
  std::string d_w0 = dir+"w0";
  std::filesystem::create_directory(d_w0);
  std::string d_w0_V = dir+"w0_V";
  std::filesystem::create_directory(d_w0_V);
  std::string d_top_ch = dir+"top_ch";
  std::filesystem::create_directory(d_top_ch);
  std::string d_top_ch_V = dir+"top_ch_V";
  std::filesystem::create_directory(d_top_ch_V);

  // --------------------------

  const Lattice lat(dimension, lattice_size);
  std::cout << lat.info() << std::endl;

  const WilsonAction S(lat,beta);
  std::cout << S.info() << std::endl;

  Rnd rnd(lat, dimension, seed_1);
  std::cout << rnd.info() << std::endl;

  Kernel Stilde(lat);
  // HMC hmc(lat, S, rnd, seed_2);
  FieldTrsf trsf(lat,Stilde,eps_W);
  FT_HMC ft_hmc(lat,S,trsf,Stilde,rnd,seed_2,nsteps_W);

  ScalarField theta(lat,dimension);
  rnd.set(theta); // initialized with Gaussian

  {
    bool is_accept = false;
    double dH = 0.0;
    double dH_mean = 0.0;

    double r_mean0 = 0.0;

    for(int k=0; k<ninit; ++k){
      ft_hmc.evolve( theta, is_accept, dH, eps, nsteps);
      r_mean0 += is_accept;
      std::clog << "dH = " << dH << std::endl;
    }

    r_mean0 /= ninit;
    std::clog << "mean0 = " << r_mean0 << std::endl;
    r_mean0 = 0.0;

    for(int k=0; k<nconf; ++k){
      std::ofstream of_w0(d_w0+"/"+std::to_string(k)+".dat");
      of_w0 << std::scientific << std::setprecision(15);
      std::ofstream of_w0_V(d_w0_V+"/"+std::to_string(k)+".dat");
      of_w0_V << std::scientific << std::setprecision(15);
      std::ofstream of_top_ch(d_top_ch+"/"+std::to_string(k)+".dat");
      of_top_ch << std::scientific << std::setprecision(15);
      std::ofstream of_top_ch_V(d_top_ch_V+"/"+std::to_string(k)+".dat");
      of_top_ch_V << std::scientific << std::setprecision(15);

      ft_hmc.evolve( theta, is_accept, dH, eps, nsteps);
      r_mean0 += is_accept;
      dH_mean += std::abs(dH);
      std::clog << "dH = " << dH << std::endl;

      {
        ScalarField theta_V( theta );
        for(int i=nsteps_W-1; i>=0; --i){
          for(int mu=lat.dim-1; mu>=0; --mu) {
            bool is_even = false;
            theta_V = trsf.inv(theta_V,is_even,mu);
            is_even = true;
            theta_V = trsf.inv(theta_V,is_even,mu);
          }
        }

        of_w0 << S.w0(theta)/lat.vol;
        of_w0_V << S.w0(theta_V)/lat.vol;
        of_top_ch << static_cast<int>(S.Q(theta));
        of_top_ch_V << static_cast<int>(S.Q(theta_V));
      }
    }
  }

  return 0;
}

