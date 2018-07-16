function s_est = WSOAV(y,H,ProjMat1,ProjMat2,gamma,K,w)
%WSOAV Weighted Sum of Absolute Values optimization
    
    [~,n]=size(H);
        
    % parameter of Douglas-Rachford algorithm
    rho=1.9;
    
    % Weighted Sum of Absolute Value Optimization via Douglas-Rachford Algorithm
    foo=ProjMat2*y;
    r=zeros(n,1);
    for k=1:K
       z=proxfw(r,gamma*ones(n,1),w);
       r=r+rho*(ProjMat1*(2*z-r+foo)-z);
    end
    
    s_est=z;
    
end