function [gBestScore,gBest,cg_curve]=PSO(N,iter,lb,ub,dim,y)

%PSO Infotmation
Vmax=ones(1,dim).*(ub-lb).*0.15;           %速度最大值
noP=N;
w=0.7;
c1=1.5;
c2=1.5;

% Initializations
vel=zeros(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
cg_curve=zeros(1,iter);

% Random initialization for agents.

lb=ones(1,dim).*(lb);                              % Lower limit for variables
ub=ones(1,dim).*(ub);                              % Upper limit for variables

for i=1:dim
    pos(:,i) = lb(i)+rand(N,1).*(ub(i) - lb(i));                          % Initial population
end


for i=1:noP
    pBestScore(i)=inf;
end

% Initialize gBestScore for a minimization problem
gBestScore=inf;


for l=1:iter
    
    % Return back the particles that go beyond the boundaries of the search space
    for i=1:size(pos,1)
        Flag4ub=pos(i,:)>ub;
        Flag4lb=pos(i,:)<lb;
        pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    
    for i=1:size(pos,1)
        %Calculate objective function for each particle
        fitness= y(pos(i,:) );
        if(pBestScore(i)>fitness)
            pBestScore(i)=fitness;
            pBest(i,:)=pos(i,:);
        end
        if(gBestScore>fitness)
            gBestScore=fitness;
            gBest=pos(i,:);
        end
    end
    
    %Update the W of PSO
  
    %Update the Velocity and Position of particles
    for i=1:size(pos,1)
        for j=1:size(pos,2)
            vel(i,j)=w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+c2*rand()*(gBest(j)-pos(i,j));
            
            if(vel(i,j)>Vmax(j))
                vel(i,j)=Vmax(j);
            end
            if(vel(i,j)<-Vmax(j))
                vel(i,j)=-Vmax(j);
            end
            pos(i,j)=pos(i,j)+vel(i,j);
        end
    end
    cg_curve(l)=gBestScore;
end

end