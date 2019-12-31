
 // in matlab try  zs = SolitaryGravityAmplitude(0.6,[-0.5*i,-0.25*i]); plot(real(zs), imag(zs)); for comparison

#include"SolitaryGravityAmplitude.m.hpp"
#include<armadillo>
#include<vector>
#include<iomanip>

int main()
{
 double Fr=1.25;
 cx_vec Z;
 Z=cx_vec(2);
 Z(0)=cx_double(0,-0.5);
 Z(1)=cx_double(0,-0.25);

 cx_vec W,F,A;
 cx_vec zs, ws, fs;
 vec P;
 vec SWP;

 SolitaryGravityAmplitude(0.6,Z,zs, ws, fs,SWP,W,F,P,A);

 cout<<std::setprecision(15);
 cout<<P<<endl;
 cout<<F<<endl;
 cout<<W<<endl;
 cout<<A<<endl;
 system("pause");
 return 0;
}
