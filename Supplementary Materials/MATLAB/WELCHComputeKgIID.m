function [Kg]=WELCHComputeKgIID(n1,n2,h,tol,method)
%Compute Kg for the Welch test, i.i.d. case
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> WELCHComputeKgIID(2,2,@(x) pdf('Normal',x,0,1),1e-6,@ quadl)
%ans =
%    1.0000

n=n1+n2;
alpha=1/n1;
beta =1/n2;

C=2*pi^((n-1)/2)*...
  ((n1-1)/alpha)^((n1-1)/2)*...
  ((n2-1)/beta)^((n2-1)/2)*...
  (1/n1+1/n2)^((n-2)/2)/...
  ((n-2)^((n-2)/2)*gamma((n-1)/2));

w0=acos(sqrt(n2/n));

R=sqrt(2*gammaincinv((2*pi*tol)/2^((n-2)/2),n/2,'upper'));

F=@(r,w) cos(w).^(n-2).*r.^(n-1).*...
         h(r*cos(w-w0)/sqrt(n1)).^n1.*...
         h(r*sin(w-w0)/sqrt(n2)).^n2;
FR=@(t,w) R./(t-1).^2.*F(R./(1-t),w); 

Kg=C*(dblquad(F,0,R,-pi/2,w0,tol,method)+...
      dblquad(F,0,R,w0,pi/2,tol,method)+...
      dblquad(FR,0,1-tol,-pi/2,w0,tol,method)+...
      dblquad(FR,0,1-tol,w0,pi/2,tol,method));