% Performs RWMH algorithm for baseline results in paper
function [theta,acceptances,h,Ws]= RWMH(S0,S1,P,F,K,L,theta0)
Plag=getLags(P,K+L+1);
Flag=getLags(F,K+L+1);
T=length(Plag);
P=P(end-T+1:end);
thetascale=[0.0162 0.0025 0.0165 0.0264 0.1951].^2; % Manually tuned to be proportional to the diagonal variance matrix of theta 
S=S0+S1;
theta=zeros(S,5);
Ws=zeros(T,S);

h=zeros(S,1);
theta(1,:)=theta0;
u=rand(S,1);
acceptances=zeros(S,1);
lambda=.0075; % Manually tuned to get rule of thumb acceptance rates
stepsize=lambda*diag(thetascale);
posttheta=@(theta,h) (-0.5*h*(P-Housing(theta,Plag,Flag,K,L))'*(P-Housing(theta,Plag,Flag,K,L)));

for s=2:S
    % Draw h
    nu1=T;
    ssq1=((P-Housing(theta(s-1,:),Plag,Flag,K,L))'*(P-Housing(theta(s-1,:),Plag,Flag,K,L)))/nu1;
    h(s)=gamrnd(nu1/2,(ssq1*nu1/2)^-1,1);

    % Draw theta
    thetastar=mvnrnd(theta(s-1,:),stepsize); % draw proposal 

    acc=exp(posttheta(thetastar',h(s))-posttheta(theta(s-1,:)',h(s))); % calculate R

    % Acceptance or rejection decision
    if u(s)<min(acc,1)
        theta(s,:)=thetastar;
    else 
        theta(s,:)=theta(s-1,:);
    end

    acceptances(s)=min(acc,1);

    Ws(:,s)=getW(theta(s,:),Plag,Flag,K,L);

if mod(s,10000)==0
fprintf('%1.0f Draws \n',s)
end

end