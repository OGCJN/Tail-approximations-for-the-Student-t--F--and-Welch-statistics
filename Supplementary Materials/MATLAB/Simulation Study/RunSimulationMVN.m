addpath('..');

%% One-sample t-test
test='OST';
n=3;
for dist=[-1 -2 -3]
 if dist==-1
     rho=0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g=@(x) mvnpdf(x,[0 0 0],[sigma1^2           rho*sigma1*sigma2  0  ;... 
                              rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
                              0                  rho*sigma2*sigma3  sigma3^2  ;]);
     Kg=OSTComputeKg(n,g,1e-6,@quadl);    
 elseif dist==-2
     rho=-0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g=@(x) mvnpdf(x,[0 0 0],[sigma1^2           rho*sigma1*sigma2  0  ;... 
                              rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
                              0                  rho*sigma2*sigma3  sigma3^2  ;]);
     Kg=OSTComputeKg(n,g,1e-6,@quadl);
 elseif dist==-3
     rho=0;
     sigma1=1;
     sigma2=2;
     sigma3=3;
     g=@(x) mvnpdf(x,[0 0 0],[sigma1^2           rho*sigma1*sigma2  0  ;... 
                              rho*sigma1*sigma2  sigma2^2           rho*sigma2*sigma3;... 
                              0                  rho*sigma2*sigma3  sigma3^2  ;]);
     Kg=OSTComputeKg(n,g,1e-6,@quadl);
 else
    return;
 end  
 for r=[1 20 100 1000 10000]
    display(['test=' test...
             ' n=' num2str(n)...
             ' dist=' num2str(dist)...
             ' r=' num2str(r)]);   
        tic;   
        OSTSimulatePValues(dist,n,r,Kg);
        toc
 end
end

%% Two-sample t-test
test='TST';
n1=2;
n2=3;
for dist=[-1 -2 -3]
 if dist==-1
     rho=0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5));        
     Kg=TSTComputeKg(n1,n2,g,1e-6,@quadl);    
 elseif dist==-2
     rho=-0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5));      
     Kg=TSTComputeKg(n1,n2,g,1e-6,@quadl); 
 elseif dist==-3
     rho=0;
     sigma1=1;
     sigma2=2;
     sigma3=3;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=TSTComputeKg(n1,n2,g,1e-6,@quadl); 
 else
    return;
 end   
 for r=[1 20 100 1000 10000]
        display(['test=' test...
                 ' n1=' num2str(n1)...
                 ' n2=' num2str(n2)...
                 ' dist=' num2str(dist)...
                 ' r=' num2str(r)]);   
            tic;   
            TSTSimulatePValues(dist,n1,n2,r,Kg);
            toc
 end
end

%% Welch test
test='WELCH';
n1=2;
n2=3;
for dist=[-1 -2 -3]
 if dist==-1
     rho=0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5));      
     Kg=WELCHComputeKg(n1,n2,g,1e-6,@quadl);    
 elseif dist==-2
     rho=-0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=WELCHComputeKg(n1,n2,g,1e-6,@quadl); 
 elseif dist==-3
     rho=0;
     sigma1=1;
     sigma2=2;
     sigma3=3;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=WELCHComputeKg(n1,n2,g,1e-6,@quadl); 
 else
    return;
 end   
 for r=[1 20 100 1000 10000]
        display(['test=' test...
                 ' n1=' num2str(n1)...
                 ' n2=' num2str(n2)...
                 ' dist=' num2str(dist)...
                 ' r=' num2str(r)]);   
            tic;   
            WELCHSimulatePValues(dist,n1,n2,r,Kg);
            toc
 end
end

%% F-test
test='F';
n1=2;
n2=3;
for dist=[-1 -2 -3]
 if dist==-1
     rho=0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=0.934199;  %FComputeKg(n1,n2,g,1e-4); % Takes too long, so we replaced it by the exact value computed in Mathematica    
 elseif dist==-2
     rho=-0.2;
     sigma1=1;
     sigma2=1;
     sigma3=1;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=1.06623;  %FComputeKg(n1,n2,g,1e-4); % Takes too long, so we replaced it by the exact value computed in Mathematica     
 elseif dist==-3
     rho=0;
     sigma1=1;
     sigma2=2;
     sigma3=3;
     g1=@(x) mvnpdf(x,[0 0],   [sigma1^2           rho*sigma1*sigma2;...
                                rho*sigma1*sigma2  sigma2^2]);
     g2=@(x) mvnpdf(x,[0 0 0], [sigma1^2             rho*sigma1*sigma2    0;...
                                rho*sigma1*sigma2    sigma2^2             rho*sigma2*sigma3;...
                                0                    rho*sigma2*sigma3    sigma3^2]);                                                       
     g=@(x)  g1(x(1:2))*g2(x(3:5)); 
     Kg=5*sqrt(3)/14;%FComputeKg(n1,n2,g,1e-4); % Takes too long, so we replaced it by the exact value computed in Mathematica     
 else
    return;
 end   
 for r=[1 20 100 1000 10000]
        display(['test=' test...
                 ' n1=' num2str(n1)...
                 ' n2=' num2str(n2)...
                 ' dist=' num2str(dist)...
                 ' r=' num2str(r)]);   
            tic;   
            FSimulatePValues(dist,n1,n2,r,Kg);
            toc
 end
end