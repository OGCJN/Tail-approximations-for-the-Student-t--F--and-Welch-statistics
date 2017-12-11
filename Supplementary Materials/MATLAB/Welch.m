function p = Welch(X,Y)

    n1 = size(X,2);
    n2 = size(Y,2);
        
    sigma1_hat = std(X,0,2)';
    sigma2_hat = std(Y,0,2)';

    nu =...
     (1/n1*sigma1_hat.^2+1/n2*sigma2_hat.^2).^2./...
    ((1/n1^2*sigma1_hat.^4)/(n1-1)+(1/n2^2*sigma2_hat.^4)/(n2-1));
    
    xbar = mean(X,2)';
    ybar = mean(Y,2)';
    
    Sp = sqrt(1/n1*sigma1_hat.^2+1/n2*sigma2_hat.^2);
    T  = (xbar-ybar)./Sp;
    
    p = 1-tcdf(T,nu);
    
end