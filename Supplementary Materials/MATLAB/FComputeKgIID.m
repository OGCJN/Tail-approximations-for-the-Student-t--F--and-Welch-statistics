function [Kg]=FComputeKgIID(n1,n2,h,tol,method)
%Compute Kg for the F-test, i.i.d. case
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> FComputeKgIID(2,3,@(x) pdf('Normal',x,0,1),1.0e-6,@ quadl)
%ans =
%    1.0000

C=gamma((n1-1)/2)*(pi*(n1-1))^((n2-1)/2)/...
  gamma((n1+n2-2)/2);

if (n1==2)
    R1=...
        sqrt(2*gammaincinv((2*pi*tol)/2^((n2-1)/2),(n2+1)/2,'upper'));
    F1=@(r,w)  r.^n2.*...
        std([cos(w); sin(w)]).^(n2-1).*...
        h(r.*cos(w)).*...
        h(r.*sin(w));
    F1R=@(t,w) R1./(t-1).^2.*F1(R1./(1-t),w); 
    C1=dblquad(F1,0,R1,0,2*pi,tol,method)+...
       dblquad(F1R,0,1-tol,0,2*pi,tol,method);
elseif (n1==3)
    R1=...
        sqrt(2*gammaincinv(((2*pi)^(3/2)*tol)/2^(n2/2),(n2+2)/2,'upper'));
    F1=@(r,w1,w2) r.^(n2+1).*sin(w1).*...
        std([cos(w1); sin(w1).*cos(w2); sin(w1).*sin(w2)]).^(n2-1).*...
        h(r.*cos(w1)).*...
        h(r.*sin(w1).*cos(w2)).*...
        h(r.*sin(w1).*sin(w2));    
    F1R=@(t,w1,w2) R1./(t-1).^2.*F1(R1./(1-t),w1,w2);    
    C1=triplequad(F1,0,R1,0,pi,0,2*pi,tol,method)+...
       triplequad(F1R,0,1-tol,0,pi,0,2*pi,tol,method);    
else
    buf=10000;
    C1=-ones(1,100);
    k=1;
    F1=@ (x) std(x).^(n2-1).*prod(h(x))./prod(pdf('Normal',x,0,1));
    x=random('Normal',0,1,n1,buf);
    V=F1(x);
    C1(1)=mean(V); 
    sigma=std(V);
    handle = waitbar(k*buf*tol^2/sigma,['Estimated standard deviation: ' num2str(sigma/sqrt(k*buf)) '.'],...
            'Name','Kg is computed using Monte Carlo simulations. This may take a while...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handle,'canceling',0);
    while sigma/sqrt(k*buf)>tol
        k=k+1;
        if mod(k,100)==1
           if getappdata(handle,'canceling')
               delete(handle);
               break;
           else
               waitbar(k*buf*tol^2/sigma,handle,['Estimated standard deviation: ' num2str((sigma/sqrt(k*buf))) '.']);
               drawnow;
           end 
           C1=[C1 -ones(1,100)]; %#ok<AGROW>
        end
        x=random('Normal',0,1,n1,buf);
        V=F1(x);
        C1(1)=mean(V); 
        sigma=sqrt(((k-1)*(buf-1)*sigma^2+(buf-1)*std(V)^2)/(k*(buf-1)));
    end
    delete(handle);
    C1(C1==-1)=[];
    C1=mean(C1);
end

R2=norminv(1-tol,0,1)/sqrt(n2);

F2=@(r) sqrt(n2)*h(r).^n2;
F2Rp=@(t) R2./(t-1).^2.*h(R2./(1-t)).^n2;
F2Rm=@(t) R2./(t-1).^2.*h(R2./(t-1)).^n2;

C2=method(F2,-R2,R2)+...
   method(F2Rp,0,1-tol,tol)+...
   method(F2Rm,0,1-tol,tol);

Kg=C*C1*C2;