                              %-------------------------------------------%
                              %         Evaluate the population           %           
                              %-------------------------------------------%                               
function [gBestScore,gBest,GlobalBestCost,gbest1,r_cal1,p_cal1]=PSOGSA(n,iteration,dataFrequencies,r_obs,low,up,dim,p_obs)

current_fitness =zeros(n,1);
gBest=zeros(1,dim);
gBestScore=inf;

for k=1:n
        pBestScore(k)=inf;
end
        pBest=zeros(n,dim);
G0=100;                                          % gravitational constant
% Boundary_no = size(up',1);
if dim==1
    current_position = rand()*(up-low)+low; %initial positions in the problem's boundary
end
% If each variable has a different up and low
if dim>1
    for k=1:dim
        up_k=up(1,k);
        low_k=low(1,k);
       current_position(1:n,k)=rand(n,1).*(up_k-low_k)+low_k; %initial positions in the problem's boundary
    end
end
velocity = 0.1*randn(n,dim) ;
acceleration=zeros(n,dim);
mass(n)=0;
force=zeros(n,dim);
C1=1.2; %C1 in Equation (9)
C2=1.2;%C2 in Equation (9)
w=0.67;%weight factor
%%main loop
iter = 0 ;                  % Iterations? counter
while  ( iter < iteration )
G=G0*exp(-20*iter/iteration); %Equation (4)
iter = iter + 1;
% acceleration=zeros(n,dim);
% mass(n)=0;
% force=zeros(n,dim);
for k = 1:n
    for j=1:dim
    fitness=0;
    %///Bound the search Space///
     Flag4up=current_position(k,j)>up(1,j);
    Flag4low=current_position(k,j)<low(1,j);
    current_position(k,j)=(current_position(k,j).*(~(Flag4up+Flag4low)))+up(1,j).*Flag4up+low(1,j).*Flag4low;
  % current_position(i,j)=(current_position(i,j).*(~(Flag4up+Flag4low)))+up.*Flag4up+low'.*Flag4low; 
    %////////////////////////////
    end
end                               %-------------------------------------------%
                                 %         Evaluate the population           %           
for k= 1:n                               %-------------------------------------------%      
    [fitness,r_cal1,p_cal1]=benchmark_functions(current_position(k,:),dataFrequencies,r_obs,p_obs);
    current_fitness(k)=fitness;    
        
    if(pBestScore(k)>fitness)
        pBestScore(k)=fitness;
        pBest(k,:)=current_fitness(k,:);
    end
    if(gBestScore>fitness)
        gBestScore=fitness;
        gBest=current_position(k,:);
    end
    
end
best=min(current_fitness);
worst=max(current_fitness);
        GlobalBestCost(iter)=gBestScore;
        GlobalBestCost(iter);
        best;
    for pp=1:n
        if current_fitness(pp)==best
            break;
        end
        
    end
    
    bestIndex=pp;
            
    for pp=1:dim
        best_fit_position(iter,1)=best;
        best_fit_position(iter,pp+1)=current_position(bestIndex,pp);   
    end
                                               %-------------------%
                                               %   Calculate Mass  %
                                               %-------------------%
    for k=1:n
    mass(k)=(current_fitness(k)-0.99*worst)/(best-worst);    
end
for k=1:n
    mass(k)=mass(k)*5/sum(mass);    
    
end
                                               %-------------------%
                                               %  Force    update  %
                                               %-------------------%
for k=1:n
    for j=1:dim
        for k1=1:n
            if(current_position(k1,j)~=current_position(k,j))
                % Equation (3)
                force(k,j)=force(k,j)+ rand()*G*mass(k1)*mass(k)*(current_position(k1,j)-current_position(k,j))/abs(current_position(k1,j)-current_position(k,j));
                
            end
        end
    end
end
                                               %------------------------------------%
                                               %  Accelations $ Velocities  UPDATE  %
                                               %------------------------------------%
for k=1:n
       for j=1:dim
            if(mass(k)~=0)
                %Equation (6)
                acceleration(k,j)=force(k,j)/mass(k);
            end
       end
end   
for k=1:n
        for j=1:dim
            %Equation(9)
            velocity(k,j)=w*velocity(k,j)+C1*rand()*acceleration(k,j) + C2*rand()*(gBest(j)-current_position(k,j));
        end
end
                                               %--------------------------%
                                               %   positions   UPDATE     %
                                               %--------------------------%
                                                        
%Equation (10) 
for k=1:n
        for j=1:dim
        current_position(k,j) = current_position(k,j) + velocity(k,j) ;
        end
end
gbest1(iter,:)=gBest;
end
end
