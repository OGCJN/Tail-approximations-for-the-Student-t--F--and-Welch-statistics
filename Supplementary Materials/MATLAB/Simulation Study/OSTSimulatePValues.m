function [PValuesRaw,PValuesCorrected]=OSTSimulatePValues(dist,n,r,Kg)
% dist    - population distribution: 
%           000-Uniform(-1,1);
%           001-Normal(0,1);
%           002-Centered exponential;
%           003-Cauchy;
%           004-t2;
%           005-t3;
%           ...
% n>1     - integer, sample size.
% r>10    - integer, only p-values less than 1/r will be considered
%
% PValuesRaw           - vector of p-values computed using t-distribution.
% PValuesCorrected     - vector of p-values computed using 
%                        the asymptotic approximation formula.
% df>0    - integer, degrees of freedom in the asymptotic formula,
%           as well as the degrees of freedom for the actual test.
% u=1/r   - threshold u. 
% res>100 - integer, approximate number of simulated p-values
%           that fall in the interval [0,1/r]. 
% Kg>0    - asymptotic constant. 
% N>1000  - total number of simulated p-values. N=r*res.
% buf     - number of simulation per cycle. Provides more effecient
%           use of memory.
% k       - number of simulation cycles done so far.

df =n-1;
u=1/r;
res=10000;

if dist==0
    Simulate_RV=@(n,m) random('Uniform',-1,1,n,m);  
elseif dist==1
    Simulate_RV=@(n,m) random('Normal',0,1,n,m); 
elseif dist==2
    Simulate_RV=@(n,m) random('Exponential',1,n,m)-1; 
elseif dist>2
    Simulate_RV=@(n,m) random('t',dist-2,n,m); 
elseif (dist==-1)&&(n==3)
    MU=[0 0 0];
    rho=0.2;
    sigma1=1;
    sigma2=1;
    sigma3=1;
    SIGMA=[sigma1^2           rho*sigma1*sigma2  0  ;... 
           rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
           0                  rho*sigma2*sigma3  sigma3^2  ;];
    Simulate_RV=@(n,m) mvnrnd(MU,SIGMA,m)'; 
elseif (dist==-2)&&(n==3)
    MU=[0 0 0];
    rho=-0.2;
    sigma1=1;
    sigma2=1;
    sigma3=1;
    SIGMA=[sigma1^2           rho*sigma1*sigma2  0  ;... 
           rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
           0                  rho*sigma2*sigma3  sigma3^2  ;];
    Simulate_RV=@(n,m) mvnrnd(MU,SIGMA,m)'; 
elseif (dist==-3)&&(n==3)
    MU=[0 0 0];
    rho=0;
    sigma1=1;
    sigma2=2;
    sigma3=3;
    SIGMA=[sigma1^2           rho*sigma1*sigma2  0  ;... 
           rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
           0                  rho*sigma2*sigma3  sigma3^2  ;];
    Simulate_RV=@(n,m) mvnrnd(MU,SIGMA,m)';     
else 
    display('This density is not supported');
    return;
end
     
N=r*res;
buf=min(N,100000);
k=0;

PValuesRaw=-1*ones(1,res);
PValuesCorrected=-1*ones(1,res);
NRaw=0;               
NCorrected=0;

while (k*buf<N)
    
    x=Simulate_RV(n,buf);
    xbar=mean(x);
    S=std(x);
    T=sqrt(n)*xbar./S;
    
    PValuesCorrected_tmp=Kg*(1-tcdf(T,df));
    PValuesCorrected_tmp=PValuesCorrected_tmp(PValuesCorrected_tmp<=u);

    if (NCorrected+length(PValuesCorrected_tmp)>length(PValuesCorrected))
        PValuesCorrected=[PValuesCorrected -1*ones(1,res)];         %#ok<AGROW>
    end
    PValuesCorrected((NCorrected+1):(NCorrected+length(PValuesCorrected_tmp)))=PValuesCorrected_tmp;
    NCorrected=NCorrected+length(PValuesCorrected_tmp);

    PValuesRaw_tmp=1-tcdf(T,df);
    PValuesRaw_tmp=PValuesRaw_tmp(PValuesRaw_tmp<=u);

    if (NRaw+length(PValuesRaw_tmp)>length(PValuesRaw))
        PValuesRaw=[PValuesRaw -1*ones(1,res)];                     %#ok<AGROW>
    end
    PValuesRaw((NRaw+1):(NRaw+length(PValuesRaw_tmp)))=PValuesRaw_tmp;
    NRaw=NRaw+length(PValuesRaw_tmp);

    k=k+1;
end

PValuesCorrected(PValuesCorrected==-1)=[];
PValuesRaw(PValuesRaw==-1)=[];

SavePValues('OST',dist,n,0,0,r,res,Kg,PValuesRaw,PValuesCorrected);