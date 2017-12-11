function [Kg]=TSTComputeKgIS(n1,n2,h1,h2,tol,method)
%Compute Kg for the Student two-sample t-test, independent samples
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> TSTComputeKgIS(4,6,@(x) mvnpdf(x,[0 0 0 0],diag([2 2 2 2])),@(x) mvnpdf(x,[0 0 0 0 0 0],diag([4 4 4 4 4 4])),1e-6,@ quadl)
%ans = 7^4/(625*2^(5/2))
%    = 0.6791

n=n1+n2;
alpha=(1/n1+1/n2)*(n1-1)/(n-2);
beta =(1/n1+1/n2)*(n2-1)/(n-2);

C=2*pi^((n-1)/2)*...
  ((n1-1)/alpha)^((n1-1)/2)*...
  ((n2-1)/beta)^((n2-1)/2)*...
  (1/n1+1/n2)^((n-2)/2)/...
  ((n-2)^((n-2)/2)*gamma((n-1)/2));

w0=acos(sqrt(n2/n));

R=sqrt(2*gammaincinv((2*pi*tol)/2^((n-2)/2),n/2,'upper'));

F=@(r,w) cos(w).^(n-2).*r.^(n-1).*...
  arrayfun(@(x) h1(x*cos(w-w0)/sqrt(n1)*ones(1,n1)),r).*...
  arrayfun(@(x) h2(x*sin(w-w0)/sqrt(n2)*ones(1,n2)),r);
FR=@(t,w) R./(t-1).^2.*F(R./(1-t),w); 

Kg=C*(dblquad(F,0,R,-pi/2,w0,tol,method)+...
      dblquad(F,0,R,w0,pi/2,tol,method)+...
      dblquad(FR,0,1-tol,-pi/2,w0,tol,method)+...
      dblquad(FR,0,1-tol,w0,pi/2,tol,method));