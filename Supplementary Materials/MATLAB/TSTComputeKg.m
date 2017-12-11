function [Kg]=TSTComputeKg(n1,n2,g,tol,method)
%Compute Kg for the Student two-sample t-test
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> TSTComputeKg(3,2,@(x) mvnpdf(x,[0 0 0 0 0],[1 1/2 -1/3 0 1/3; 1/2 2 1/6 0 0; -1/3 1/6 3 0 1/5; 0 0 0 4 1/2; 1/3 0 1/5 1/2 5;]),1e-6,@ quadl)
%ans = 13585329*sqrt(1509481/5)/3993677192
%    = 1.8691

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
  arrayfun(@(x) g(x*cos(w-w0)/sqrt(n1)*[ones(1,n1) zeros(1,n2)]+...
    x*sin(w-w0)/sqrt(n2)*[zeros(1,n1) ones(1,n2)]),r);
FR=@(t,w) R./(t-1).^2.*F(R./(1-t),w); 

Kg=C*(dblquad(F,0,R,-pi/2,w0,tol,method)+...
      dblquad(F,0,R,w0,pi/2,tol,method)+...
      dblquad(FR,0,1-tol,-pi/2,w0,tol,method)+...
      dblquad(FR,0,1-tol,w0,pi/2,tol,method));