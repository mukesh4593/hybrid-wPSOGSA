function [X]=initialization(dim,N,up,down)

if size(up,2)==1
    X=rand(N,dim).*(up-down)+down; % Initialization of each parameter with random number
end
if size(up,2)>1
    for k=1:dim
    high=up(k);low=down(k);
    X(:,k)=rand(N,1).*(high-low)+low;
    end
end
end