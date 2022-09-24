%Code for the paper "Optimal control of energy storage under random
% operation permissions, IISE Transactions, Volume 50, 2018 - Issue 8: 
% Energy Systems Modeling and Analytics ", Somayeh Moazeni
%This is the main code to get the optimal value function and actions. 
%The code calls data from "PermissionArrivalRate.m".
%This function can be called as:
% [V, aGrid]=vGrid_TimeVarying('maxP','pt')
% or
% [V, aGrid]=vGrid_TimeVarying('maxP','p')
% [V, aGrid]=vGrid_TimeVarying('max','pt')
% [V, aGrid]=vGrid_TimeVarying('max','p')

function [V, aGrid]=vGrid_TimeVarying(ch,EPch)
close all;
K=4;      %Capacity of the Storage Unit
T=16;     %Planning Time horizon 
dt=1/12;  % Time discritization
dk=1/8;   % Capacity discritization
Tsteps=T/dt;
Ksteps=K/dk;
lambda=24.9355/16;   %Permission Arrival Rate
EndOfHorizonF=0;
V=zeros(Ksteps+1,Tsteps+1);
aGrid=zeros(Ksteps,Tsteps);

if strcmp(EPch,'pt')==1
    EP=PermissionArrivalRate(50);
elseif strcmp(EPch,'p')==1
    EP=PermissionArrivalRate(50);
    EP=mean(EP)*ones(size(EP));
end
MEP=length(EP);
if strcmp(ch,'maxP')==1
    F=@(x,p) log(x*p+1);   %log utility
    %F=@(x,p) x-0.5*p*x*x; 
    %F=@(x,p) p-(1/(1+x)); %power reward function
    %F=@(x,p) x-p*(x-2)*(x-2)/2; %power reward function
    %boundary values
    for ik=1:(Ksteps+1)
       V(ik,1)=F((Ksteps+1-ik)*dk,EndOfHorizonF); %corresponding to t=0;
    end
    for it=2:(Tsteps+1)
       V(Ksteps+1,it)=0;  %corresponding to k=0;
       aGrid(Ksteps+1,it)=0;
       V(Ksteps,it)=(1-dt*lambda)*V(Ksteps,it-1)+dt*lambda*max(F(1*dk,EP(MEP-(it-2)))+V(Ksteps+1,it-1), V(Ksteps,it-1)); %corresponding to k=1;
       maxval=V(Ksteps,it-1);
       aGrid(Ksteps,it)=0;
       if (F(1*dk,EP(MEP-(it-2)))+V(Ksteps+1,it-1))>V(Ksteps,it-1)
         maxval=(F(1*dk,EP(MEP-(it-2)))+V(Ksteps+1,it-1));
         aGrid(Ksteps,it)=1*dk;
       end       
    end
    for ik=(Ksteps-1):-1:1 
       for it=2:Tsteps+1 % at each row we move from left to right
           maxVal=V(ik,it-1);
           aGrid(ik,it)=0;
           for a=1:(Ksteps-ik)
               if (V(ik+a,it-1)+F(a*dk,EP(MEP-it+2)))>maxVal
                   maxVal=V(ik+a,it-1)+F(a*dk,EP(MEP-it+2));
                   aGrid(ik,it)=a*dk;
               end
           end
           if F((Ksteps+1-ik)*dk,EP(MEP-(it-2)))>maxVal
               maxVal=F((Ksteps+1-ik)*dk,EP(MEP-(it-2)));
               aGrid(ik,it)=(Ksteps+1-ik)*dk;
           end
           V(ik,it)=(1-dt*lambda)*V(ik,it-1)+dt*lambda*maxVal;
       end
    end 
elseif strcmp(ch,'max')==1
    F=@(x) log(x+1);
    % boundary values
    for ix=1:(Ksteps+1)
       V(ix,1)=F((Ksteps-ix+1)*dk); %corresponding to t=0;
    end
    for iy=2:(Tsteps+1)
       V(Ksteps+1,iy)=0;  %corresponding to k=0;
       V(Ksteps,iy)=max(F(1*dk),V(Ksteps,iy-1)); %corresponding to k=1;
    end
    for ix=(Ksteps-1):-1:1 
       for iy=1:Tsteps % at each row we move from left to right
           maxVal=F((Ksteps+1-ix)*dk);  
           for a=1:(Ksteps-ix)
               maxVal=max(V(ix+a,iy)+F(a*dk),maxVal);
           end
           V(ix,iy+1)=(1-dt*lambda)*V(ix,iy)+dt*lambda*max(V(ix,iy),maxVal);
       end
    end 
    for k=1:Ksteps
        for t=1:Tsteps
            vec=zeros(k,1);                     
            for a=1:k
                vec(a)=V(Ksteps+1-(k-a),t)+F(a*dk);
            end
         [w1,akT]=max(vec);
         aGrid(k,t)=akT;
        end
    end
elseif strcmp(ch,'min')==1
    F=@(x) 0.5*x*x; 
    % boundary values
    for ix=1:(Ksteps+1)
       V(ix,1)=F((Ksteps-ix+1)*dk); %corresponding to t=0;
    end
    for iy=2:(Tsteps+1)
       V(Ksteps+1,iy)=0;  %corresponding to k=0;
       V(Ksteps,iy)=min(F(1*dk),V(Ksteps,iy-1)); %corresponding to k=1;
    end
   for ix=(Ksteps-1):-1:1 
       for iy=1:Tsteps 
           minVal=F((Ksteps+1-ix)*dk);  
           for a=1:(Ksteps-ix)
               minVal=min(V(ix+a,iy)+F(a),minVal);
           end
           V(ix,iy+1)=(1-dt*lambda)*V(ix,iy)+dt*lambda*min(V(ix,iy),minVal);
       end
   end 
   for k=1:Ksteps
        for t=1:Tsteps
            vec=zeros(k,1);                     
            for a=1:k
                vec(a)=V(Ksteps+1-(k-a),t)+F(a*dk);
            end
         [w1,akT]=min(vec);
         aGrid(k,t)=akT;
        end 
   end
end


VKT=zeros(Ksteps+1,Tsteps+1);
for ix=1:(Ksteps+1)
    VKT(ix,:)=V(Ksteps+2-ix,:);
end
for it=1:(Tsteps+1)
    V(:,it)=VKT(:,Tsteps+2-it);
end
V;
mesh(V)
ylabel('inventory level k','FontSize',14,'FontWeight','bold');
xlabel('time t','FontSize',14,'FontWeight','bold');
zlabel('value function','FontSize',14,'FontWeight','bold');
h = get(gca,'xlabel'); set(h,'Rotation',20);
h = get(gca,'ylabel'); set(h,'Rotation',-30);

figure;
aKT=zeros(Ksteps+1,Tsteps+1);
for ix=1:(Ksteps+1)
    aKT(ix,:)=aGrid(Ksteps+2-ix,:);
end
for it=1:(Tsteps+1)
    aGrid(:,it)=aKT(:,Tsteps+2-it);
end
mesh(aGrid)
zlabel('discharge action','FontSize',14,'FontWeight','bold');
ylabel('inventory level k','FontSize',14,'FontWeight','bold');
xlabel('time t','FontSize',14,'FontWeight','bold');
h = get(gca,'xlabel'); set(h,'Rotation',20);
h = get(gca,'ylabel'); set(h,'Rotation',-30);
aGrid;




