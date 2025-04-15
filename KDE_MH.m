% Perfoms RWMH with KDE estimate of Likelihood. For the most part identical
% to RWMH code except when calculating the Likelihood
% (iter=N in the paper)
function [theta,acceptances,h]= KDE_MH(S0,S1,P,F,K,L,theta0,iter,C)
Plag=getLags(P,K+L+1);
Flag=getLags(F,K+L+1);
T=length(Plag);

thetascale=[0.0162 0.0025 0.0165 0.0264 0.1951].^2;
S=S0+S1;
theta=zeros(S,5);
h=zeros(S,1);
theta(1,:)=theta0;
u=rand(S,1);
acceptances=zeros(S,1);
lambda=.0075;
stepsize=lambda*diag(thetascale);
seeds=randi([1 2^32-1],S,1);
for s=2:S
    % Draw h (Gibbs step)
    nu1=T;
    ssq1=((P(end-T+1:end)-Housing(theta(s-1,:),Plag,Flag,K,L))'*(P(end-T+1:end)-Housing(theta(s-1,:),Plag,Flag,K,L)))/nu1;
    h(s)=gamrnd(nu1/2,(ssq1*nu1/2)^-1,1);

    % Draw theta 
    thetastar=mvnrnd([theta(s-1,:)],stepsize);

    rng(seeds(s)) % Fix seeds before simulation 
    postthetaold= SimulLik(theta(s-1,:),h(s),P,F,K,L,iter,C); % Log-Likelihood under proposed draw
    rng(seeds(s))
    postthetaprop= SimulLik(thetastar,h(s),P,F,K,L,iter,C); % Log-Likelihood under old draw
    acc=exp(postthetaprop-postthetaold);


    if u(s,1)<min(acc,1)
        theta(s,:)=thetastar;
    else 
        theta(s,:)=theta(s-1,:);
    end
    acceptances(s,1)=min(acc,1);


if mod(s,100)==0
    fprintf('%1.0f draws \n',s)
end

end