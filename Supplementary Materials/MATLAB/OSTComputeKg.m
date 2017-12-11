function [Kg]=OSTComputeKg(n,g,tol,method)
%Compute Kg for the Student one-sample t-test
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> OSTComputeKg(3,@(x) mvnpdf(x,[0 0 0],[1 1/2 -1/3; 1/2 2 1/6; -1/3 1/6 3;]),1.0e-6,@ quadl)
%
%ans = 267/250
%    = 1.068

C=2*pi^(n/2)/gamma(n/2);

R=sqrt(2*gammaincinv((2*pi*tol)/2^((n-2)/2),n/2,'upper'));

F=@(r) r.^(n-1).*arrayfun(@(x) g(x/sqrt(n).*ones(1,n)),r);
FR=@(t) R./(t-1).^2.*F(R./(1-t)); 

Kg=C*(method(F,0,R,tol)+method(FR,0,1-tol,tol));