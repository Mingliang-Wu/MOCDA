function y = MMF1_e(x)            
    [N,D] = size(X);
    M     = 2;   
    S     = ceil(0.1*(D-M));
    g     = zeros(N,4);
    for i = 1 : 4
        g(:,i) = sum(g1(X(:,M+(i-1)*S:M+i*S-1),pi/3),2)+sum(g2(X(:,[M:M+(i-1)*S-1,M+i*S:end]),0),2);
    end
    PopObj = repmat(1+min(g,[],2)/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),X(:,1:M-1)],2)).*[ones(N,1),1-X(:,M-1:-1:1)];
end