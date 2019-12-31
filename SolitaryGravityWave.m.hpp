// Automatically translated using m2cpp 2.0 on 2019-12-30 22:40:32

#ifndef SOLITARYGRAVITYWAVE_M_HPP
#define SOLITARYGRAVITYWAVE_M_HPP

#include <iostream>
#include <cmath>
#include "matlab2cpp.h"
#include "SolitaryTools.hpp"
#include <cstdio>
#include <complex>
#include <armadillo>

// linker option -llapack -lblas -lgfortran

/*
% SOLITARYGRAVITYWAVE: Computes the steady irrotational surface solitary gravity
%    wave solution of the Euler equations (homogeneous, incompressible and
%    perfect fluids). The wave is defined by its Froude number Fr and the result
%    is about fifteen digits accurate. The method works for all but the highest
%    waves, i.e. for all amplitude/depth ratio less than 0.796.
%
% SYNOPSIS:
% [zs,ws,fs,SWP] = SolitaryGravityWave(Fr);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z,1);
%
% INPUT:
% Fr : Froude number (i.e., c/sqrt(gd), c.f. Note 1 below).
% Z  : Complex abscissa where fields are desired inside the fluid (default Z = []).
%      Z should be strictly below the surface, i.e., -1 <= imag(Z) < eta(real(Z))
%      y = eta(x) being the equation of the free surface.
%
% OUTPUT (dimensionless quantities):
% zs   : Complex abscissa at the surface, i.e., x + i*eta.
% ws   : Complex velocity at the surface, i.e., u - i*v.
% fs   : Complex potential at the surface, i.e., phi + i*psi.
% SWP  : Solitary Wave Parameters, i.e.
%        SWP(1) = wave amplitude, max(eta)
%        SWP(2) = wave mass
%        SWP(3) = circulation
%        SWP(4) = impulse
%        SWP(5) = kinetic energy
%        SWP(6) = potential energy
% W    : Complex velocity in the bulk at abscissas Z.
% F    : Complex potential in the bulk at abscissas Z.
% P    : Pressure in the bulk at abscissas Z.
% A    : Complex acceleration in the bulk at abscissas Z (A = dW / dt).
%
% NOTES:
% 1- Computations are performed in dimensionless units such that rho = g = d = 1,
%    where rho is the fluid density, g is the acceleration due to gravity and d
%    is the constant water depth. It means that all lengths are scaled by d,
%    accelerations by g, speeds by sqrt(g*d), times by sqrt(d/g) and stresses
%    by rho*g*d.
% 2- The numerical scheme is based on the Petviashvili's iterations of the
%    Babenko equation using a pseudo-spectral method.
% 3- The solution is obtained in parametric form resulting from a conformal
%    mapping into a strip.
% 4- The algorithm is efficient for all but the highest waves.
% 5- W, F, P and A are column vectors corresponding to the shape of Z(). The
%    velocities are computed in the frame of reference where the fluid is
%    at rest in the far field.
% 6- For Z very close to the free surface, the fields may be inaccurate.
%
% Have a look at the m-file for more details.

% WARNING: This program was written to illustrate the Petviashvili method
%          for the Babenko equation. It is not designed to compute the wave
%          of limiting height. The code is freely distributed since it may be
%          useful to someone, but it is given 'as such' without any guarantee
%          of any kind. Suggestions for improvements and bug reports are
%          nonetheless welcome.

% Author 1: Denys Dutykh, CNRS & UCD
% E-mail  : Denys.Dutykh@ucd.ie
% URL     : http://www.denys-dutykh.com/
%
% Author 2: Didier Clamond, University of Nice - Sophia Antipolis
*/
using namespace arma ;

void SolitaryGravityWave(double Fr, cx_vec Z, cx_vec& zs, cx_vec& ws, cx_vec& fs, vec& SWP, cx_vec& W, cx_vec& F, vec& P, cx_vec& A)
{
  cx_double S, S1, S2 ;
  cx_vec IOP, Leta_hat, Nl_hat, dW, dzsm1, eta2_hat, eta_hat ;
  double F2, L, dxi, err, etaMean, gam, kd, tol ;
  int N, n ;
  uword LZ, Niter ;
  vec COP, Ceta, Dif_err, LOI, LOP, ONES, Res_err, dexi, eta, eta1, k, phis, qs, us, vs, xi, xs ;
  if (Fr>1.29421 || Fr<=1)
  {
    std::cerr << "The Froude number must be in the interval 1 < Fr <= 1.29421" << std::endl ;
  }
  F2 = pow(Fr, 2) ;
  N = 16384 ;
  tol = 1e-15 ;
  gam = 2.0 ;
  kd = Getkd(F2,0.,1.5,tol) ;
  L = 17*std::log(10)/kd ;
  dxi = 2*L/N ;
  xi = arma::trans((m2cpp::fspan(1-N/2.0, 1, N/2.0)))*dxi ;
  k = arma::trans(arma::join_rows(m2cpp::fspan(0, 1, N/2.0-1), m2cpp::fspan(-N/2.0, 1, -1)))*datum::pi/L ;
  COP = k/tanh(k) ;
  COP(0) = double(1) ;
  LOP = F2*COP-1 ;
  LOI = 1/LOP ;
  ONES = arma::ones<vec>(k.n_rows, k.n_cols) ;
  IOP = -cx_double(0, 1)*(ONES/tanh(k)) ;
  IOP(0) = cx_double(0) ;
  eta = (F2-1)*arma::square((ONES/arma::cosh(0.5*kd*xi))) ;
  err = datum::inf ;
  Res_err.reset() ;
  Dif_err.reset() ;
  while (err>tol)
  {
    eta_hat = arma::fft(eta) ;
    eta2_hat = arma::fft(arma::square(eta)) ;
    Ceta = arma::real(arma::ifft<cx_vec>(COP%eta_hat)) ;
    Leta_hat = LOP%eta_hat ;
    Nl_hat = 0.5*COP%eta2_hat+arma::fft(eta%Ceta) ;
    S1 = arma::as_scalar(arma::trans(eta_hat)*Leta_hat) ;
    S2 = arma::as_scalar(arma::trans(eta_hat)*Nl_hat) ;
    S = pow((S1/S2), gam) ;
    eta1 = arma::real(arma::ifft<cx_vec>(S*LOI%Nl_hat)) ;
    err = norm(eta1-eta, "inf") ;
    eta = eta1 ;
    Dif_err = arma::join_cols(Dif_err, m2cpp::scol<double>(err)) ;
    Res_err = arma::join_cols(Res_err, m2cpp::scol<double>(norm(arma::real(arma::ifft<cx_vec>(Leta_hat-Nl_hat)), "inf"))) ;
  }
  Ceta = arma::real(arma::ifft<cx_vec>(COP%arma::fft(eta))) ;
  dexi = arma::real(arma::ifft<cx_vec>(cx_double(0, 1)*k%arma::fft(eta))) ;
  SWP=vec(6);
  SWP(0) = arma::max(eta) ;
  SWP(1) = arma::as_scalar(dxi*arma::trans(eta)*(1+Ceta)) ;
  SWP(2) = Fr*dxi*double(arma::as_scalar(arma::sum(Ceta))) ;
  SWP(3) = Fr*SWP(1) ;
  SWP(4) = arma::as_scalar(0.5*F2*dxi*arma::trans(eta)*Ceta) ;
  SWP(5) = arma::as_scalar(0.5*dxi*arma::trans((arma::square(eta)))*(1+Ceta)) ;
  etaMean = mean(eta) ;
  xs = (1+etaMean)*xi+arma::real(arma::ifft<cx_vec>(IOP%arma::fft(eta-etaMean))) ;
  phis = Fr*(xs-xi) ;
  qs = arma::square((1+Ceta))+arma::square(dexi) ;
  us = Fr-Fr*(1+Ceta)/qs ;
  vs = -Fr*dexi/qs ;
  zs = xs+cx_double(0, 1)*eta ;
  ws = us-cx_double(0, 1)*vs ;
  fs = phis+cx_double(0, 1)*eta*Fr ;
  Niter = m2cpp::length(Dif_err) ;
  std::printf("+-------------------------------------------------+\n") ;
  std::printf("| Froude Number           = %15.14lf      |\n", Fr) ;
  std::printf("| Amplitude               = %15.14lf      |\n", SWP(0)) ;
  std::printf("| Mass                    = %15.14lf      |\n", SWP(1)) ;
  std::printf("| Circulation             = %15.14lf      |\n", SWP(2)) ;
  std::printf("| Impulse                 = %15.14lf      |\n", SWP(3)) ;
  std::printf("| Kinetic Energy          = %15.14lf      |\n", SWP(4)) ;
  std::printf("| Potential Energy        = %15.14lf      |\n", SWP(5)) ;
  std::printf("|                                                 |\n") ;
  std::printf("| Convergence achieved in %5d iterations.       |\n", Niter) ;
  std::printf("| Error between two latest iterations: %5.4g |\n", err) ;
  std::printf("| Residual                           : %5.4g |\n", Res_err(Res_err.n_rows-1)) ;
  std::printf("+-------------------------------------------------+\n") ;

  if (m2cpp::isempty(Z)==0)
  {
    LZ = m2cpp::length(Z) ;
    W = arma::zeros<cx_vec>(LZ) ;
    F = W ;
    dW = W ;
    P = arma::zeros<vec>(LZ) ;
    dzsm1 = Ceta+cx_double(0, 1)*dexi ;
    for (n=1; n<=LZ; n++)
    {
      W(n-1) = cx_double(arma::as_scalar(arma::sum(dzsm1/(zs-Z(n-1))-arma::conj(dzsm1)/(arma::conj(zs)-cx_double(0, 2)-Z(n-1))))) ;
      W(n-1) = cx_double(0, 1)*(cx_double) 0.5*(cx_double) Fr/datum::pi*(cx_double) W(n-1)*(cx_double) dxi ;
      dW(n-1) = cx_double(arma::as_scalar(arma::sum(dzsm1/arma::square((zs-Z(n-1)))-arma::conj(dzsm1)/arma::square((arma::conj(zs)-cx_double(0, 2)-Z(n-1)))))) ;
      dW(n-1) = cx_double(0, 1)*(cx_double) 0.5*(cx_double) Fr/datum::pi*(cx_double) dW(n-1)*(cx_double) dxi ;
      F(n-1) = cx_double(arma::as_scalar(arma::sum(dzsm1%arma::log((zs+cx_double(0, 1))/(zs-Z(n-1)))-arma::conj(dzsm1%arma::log((zs+cx_double(0, 1))/(zs+cx_double(0, 2)-std::conj(Z(n-1)))))))) ;
      F(n-1) = cx_double(0, 1)*(cx_double) 0.5*(cx_double) Fr/datum::pi*(cx_double) F(n-1)*(cx_double) dxi ;
    }
    P = Fr*arma::real(W)-0.5*arma::square(abs(W))-arma::imag(Z) ;
    A = dW%(arma::conj(W)-Fr) ;
  }
}
#endif
