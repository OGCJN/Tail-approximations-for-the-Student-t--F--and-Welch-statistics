function [Kg]=FComputeKg(n1,n2,g,tol)
%Compute Kg for the F-test, 
%(c) Dmitrii Zholud 2011
%
%Example: 
%>> FComputeKg(3,3,@(x) mvnpdf(x,[0 0 0 0 0 0],1/6*[1 1/2 -1/3 0 1/3 0; 1/2 2 1/6 0 0 0; -1/3 1/6 3 0 1/5 0; 0 0 0 4 1/2 1/5; 1/3 0 1/5 1/2 5 0; 0 0 0 1/5 0 6]),1e-4)
%ans = 1118038/(771147*sqrt(13))
%    = 0.4021

C=gamma((n1-1)/2)*(pi*(n1-1))^((n2-1)/2)/...
  gamma((n1+n2-2)/2);

buf=10000;
C1=-ones(1,100);
k=1;
F1=@ (x) std(x(1:(end-1),:)).^(n2-1).*cellfun(@(z)g([z(1:(end-1)); 1/sqrt(n2)*repmat(z(end),n2,1)]'),mat2cell(x,n1+1,ones(1,buf)))./prod(pdf('Normal',x,0,1));
x=random('Normal',0,1,n1+1,buf);
V=F1(x);
C1(1)=mean(V); 
sigma=std(V);
handle = waitbar(k*buf*tol^2/sigma^2,['Estimated standard deviation: ' num2str(sigma/sqrt(k*buf)) '.'],...
        'Name','Kg is computed using Monte Carlo simulations. This may take a while...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
setappdata(handle,'canceling',0);
while sigma/sqrt(k*buf)>tol
    k=k+1;
    if mod(k,10)==1
       if getappdata(handle,'canceling')
           delete(handle);
           break;
       else
           waitbar(k*buf*tol^2/sigma^2,handle,['Estimated standard deviation: ' num2str((sigma/sqrt(k*buf))) '.']);
           drawnow;
       end 
       C1=[C1 -ones(1,10)]; %#ok<AGROW>
    end
    x=random('Normal',0,1,n1+1,buf);
    V=F1(x);
    C1(k)=mean(V); 
    sigma=sqrt(((k-1)*(buf-1)*sigma^2+(buf-1)*std(V)^2)/(k*(buf-1)));
end
delete(handle);
C1(C1==-1)=[];
C1=mean(C1);

Kg=C*C1;