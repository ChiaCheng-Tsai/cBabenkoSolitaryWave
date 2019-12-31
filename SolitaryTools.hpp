// Automatically translated using m2cpp 2.0 on 2019-12-31 11:31:24

// by Tsai 2019, 12, 31

#ifndef SolitaryTools_HPP
#define SolitaryTools_HPP

double SolitarySinc(double x)
{
  if(x==0)
   return 1.;
  else
   return sin(x)/x;
}

double Getkd(double F2,double a,double b,double tol)
{
 double c;
 if(a>b)
 {
  c=a;
  a=b;
  b=c;
 }

 while(b-a>tol)
 {
  c=(a+b)/2.;
  if((SolitarySinc(a)-F2*cos(a))*(SolitarySinc(c)-F2*cos(c))<0.)
   b=c;
  else
   a=c;
 }
 return c;
}

#endif
