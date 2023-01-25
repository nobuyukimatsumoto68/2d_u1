#ifndef IS_FTHMC
#define IS_FTHMC 0
#endif

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
#ifdef _OPENMP
#include <omp.h>
#endif

using Complex = std::complex<double>;
using uint = unsigned int;
using Idx = std::size_t;

#include "util.h"
#include "lattice.h"
#include "coord.h"
#include "scalar_field.h"
#include "action.h"
#include "corr.h"
#include "rnd.h"
#include "hmc.h"
#include "kernel.h"
#include "field_trsf.h"
#include "ft_hmc.h"

// ---
const uint dimension = 2; // FIX
// ---
const double a = 1.0/1.8;
const double beta = 1.0/(a*a);
const std::vector<uint> lattice_size = std::vector<uint>{8,8};
// ---
const std::string dir = "./results/";
// ----
const double stot = 1.0;
const int nsteps = 6;
const double eps = stot/nsteps;
const int nconf = 500;
// ----
const int seed_1 = 1;
const int seed_2 = 2;

#if IS_FTHMC
const double eps_W = 0.2;
const int nsteps_W = 1;
#endif


int main(){

  std::cout << std::scientific << std::setprecision(15);
  std::clog << std::scientific << std::setprecision(15);
  std::cerr << std::scientific << std::setprecision(15);
  std::cout << std::endl << std::endl;

#ifdef _OPENMP
  omp_set_num_threads(4);
#endif

  std::filesystem::create_directory(dir);
  std::string d_w0 = dir+"w0";
  std::filesystem::create_directory(d_w0);
  std::string d_top_ch = dir+"top_ch";
  std::filesystem::create_directory(d_top_ch);
  std::string d_corr = dir+"corr";
  std::filesystem::create_directory(d_corr);

#if IS_FTHMC
  std::string d_w0_V = dir+"w0_V";
  std::filesystem::create_directory(d_w0_V);
  std::string d_top_ch_V = dir+"top_ch_V";
  std::filesystem::create_directory(d_top_ch_V);
#endif

  // --------------------------

  const Lattice lat(dimension, lattice_size);
  std::clog << lat.info() << std::endl;

  const WilsonAction S(lat,beta);
  std::clog << S.info() << std::endl;
  const Corr corr(lat);

  Rnd rnd(lat, dimension, seed_1);
  std::clog << rnd.info() << std::endl;


#if IS_FTHMC
  Kernel Stilde(lat);
  FieldTrsf trsf(lat,Stilde,eps_W);
  FT_HMC ft_hmc(lat,S,trsf,Stilde,rnd,eps,nsteps,seed_2,nsteps_W);
  std::clog << ft_hmc.info() << std::endl;
#else
  HMC hmc(lat, S, rnd, eps, nsteps, seed_2);
  std::clog << hmc.info() << std::endl;
#endif

  /* ------------- */

  {
    std::clog << "--- misc ---" << std::endl
              << "nconf = " << nconf  << std::endl
              << "--------------------" << std::endl;
  }

  ScalarField phi(lat,dimension);
  rnd.set(phi); // initialized with Gaussian

  double r_mean=0.0;
  bool is_accept=false;
  double dH=0.0;

  std::ofstream of_acc(dir+"accept_reject.dat");
  std::ofstream of_dH(dir+"dH.dat");
  of_dH << std::scientific << std::setprecision(15);

  for(int k=0; k<nconf; ++k){
    std::ofstream of_w0(d_w0+"/"+std::to_string(k)+".dat");
    of_w0 << std::scientific << std::setprecision(15);
    std::ofstream of_top_ch(d_top_ch+"/"+std::to_string(k)+".dat");
    of_top_ch << std::scientific << std::setprecision(15);
    std::ofstream of_corr(d_corr+"/"+std::to_string(k)+".dat");
    of_corr << std::scientific << std::setprecision(15);

#if IS_FTHMC
    std::ofstream of_w0_V(d_w0_V+"/"+std::to_string(k)+".dat");
    of_w0_V << std::scientific << std::setprecision(15);
    std::ofstream of_top_ch_V(d_top_ch_V+"/"+std::to_string(k)+".dat");
    of_top_ch_V << std::scientific << std::setprecision(15);
#endif


#if IS_FTHMC
    ft_hmc.evolve( phi, is_accept, dH);
#else
    hmc.evolve( phi, is_accept, dH);
#endif

    /* ------------- */

    r_mean += is_accept;
    of_acc << is_accept << ' ';
    of_dH << dH << ' ';

    of_w0 << corr.w0(phi)/lat.vol;
    of_top_ch << corr.Q(phi);

    for(uint dx=0; dx<lat.size[0]; ++dx){
      of_corr << corr(phi, dx) << ' ';
    }

#if IS_FTHMC
    ScalarField phi_V( phi );
    for(int i=nsteps_W-1; i>=0; --i){
      for(int mu=lat.dim-1; mu>=0; --mu) {
        bool is_even = false;
        phi_V = trsf.inv(phi_V,is_even,mu);
        is_even = true;
        phi_V = trsf.inv(phi_V,is_even,mu);
      }
    }
    of_w0_V << corr.w0(phi_V)/lat.vol;
    of_top_ch_V << corr.Q(phi_V);
#endif

  }
  r_mean /= nconf;
  std::clog << "acceptance = " << r_mean << std::endl;

  return 0;
}

