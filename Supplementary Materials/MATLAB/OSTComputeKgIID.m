function [Kg]=OSTComputeKgIID(n,h,tol,method)
%Compute Kg for the Student one-sample t-test, i.i.d. case
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> OSTComputeKgIID(2,@(x) pdf('Normal',x,0,1),1.0e-6,@ quadl)
%ans =
%    1.0000

C=2*pi^(n/2)/gamma(n/2);

R=sqrt(2*gammaincinv((2*pi*tol)/2^((n-2)/2),n/2,'upper'));

F=@(r) r.^(n-1).*h(r/sqrt(n)).^n;
FR=@(t) R./(t-1).^2.*F(R./(1-t)); 

Kg=C*(method(F,0,R,tol)+method(FR,0,1-tol,tol));