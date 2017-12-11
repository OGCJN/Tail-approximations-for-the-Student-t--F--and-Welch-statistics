%% OBS: It took more than a month to run the script on a regular PC.
addpath('..');

%% One-sample t-test
test='OST';
for n=[2 3 5]  
  for dist=[0 2 3]
     if dist==0
         h=@(x) pdf('Uniform',x,-1,1);
         Kg=OSTComputeKgIID(n,h,1e-6,@quad);
     elseif dist==2
         h=@(x) pdf('Exponential',x+1,1);
         Kg=OSTComputeKgIID(n,h,1e-6,@quadl);
     elseif dist==3
         h=@(x) pdf('t',x,dist-2);
         Kg=OSTComputeKgIID(n,h,1e-6,@quadl);
     else
         return;
     end  
    for r=[1 20 100 1000 10000 100000 1000000]
        display(['test=' test...
                 ' n=' num2str(n)...
                 ' dist=' num2str(dist)...
                 ' r=' num2str(r)]);   
            tic;   
            OSTSimulatePValues(dist,n,r,Kg);
            toc
    end
  end
end

%% Two-sample t-test
test='TST';
for n2=[2 3 5]
  if (n2==2)  
      n1=2;
  elseif (n2==3)
      n1=2;
  elseif (n2==5)
      n1=3;
  else
      return;
  end
  for dist=[0 2 4]
     if dist==0
         h=@(x) pdf('Uniform',x,-1,1);
         Kg=TSTComputeKgIID(n1,n2,h,1e-6,@quad);
     elseif dist==2
         h=@(x) pdf('Exponential',x+1,1);
         Kg=TSTComputeKgIID(n1,n2,h,1e-6,@quadl);
     elseif dist==4
         h=@(x) pdf('t',x,dist-2);
         Kg=TSTComputeKgIID(n1,n2,h,1e-6,@quadl);
     else
         return;
     end  
    for r=[1 20 100 1000 10000 100000 1000000]
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
end

%% Welch t-test
test='WELCH';
for n2=[2 3 5]
  if (n2==2)  
      n1=2;
  elseif (n2==3)
      n1=2;
  elseif (n2==5)
      n1=3;
  else
      return;
  end
  for dist=1
     if dist==0
         h=@(x) pdf('Uniform',x,-1,1);
         Kg=WELCHComputeKgIID(n1,n2,h,1e-6,@quad);
     elseif dist==1
         h=@(x) pdf('Normal',x,0,1);
         Kg=WELCHComputeKgIID(n1,n2,h,1e-6,@quadl);
     elseif dist==2
         h=@(x) pdf('Exponential',x+1,1);
         Kg=WELCHComputeKgIID(n1,n2,h,1e-6,@quadl);
     elseif dist==4
         h=@(x) pdf('t',x,dist-2);
         Kg=WELCHComputeKgIID(n1,n2,h,1e-6,@quadl);
     else
         return;
     end  
    for r=[1 20 100 1000 10000 100000 1000000]
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
end

%% F-test
test='F';
for n2=[5 3 2 ]
  if (n2==2)  
      n1=2;
  elseif (n2==3)
      n1=2;
  elseif (n2==5)
      n1=3;
  else
      return;
  end
  for dist=[0 2 7]
     if dist==0
         h=@(x) pdf('Uniform',x,-1,1);
         Kg=FComputeKgIID(n1,n2,h,1e-6,@quad);
     elseif dist==2
         h=@(x) pdf('Exponential',x,1);
         Kg=FComputeKgIID(n1,n2,h,1e-6,@quadl);
     elseif dist==7
         h=@(x) pdf('t',x,dist-2);
         Kg=FComputeKgIID(n1,n2,h,1e-6,@quadl);
     else
         return;
     end  
    for r=[1 20 100 1000 10000 100000 1000000]
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
end

