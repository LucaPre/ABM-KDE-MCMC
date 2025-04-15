function LogLikelihood = SimulLik(theta,h,Pfix,Ffix,K,L,iter,C)
c=theta(1);
d=theta(2);
alpha=theta(3);
beta=theta(4);
gamma_param=theta(5);
sigma=sqrt(h^-1);
Plag=getLags(Pfix,K+L+1);
Flag=getLags(Ffix,K+L+1);
T=length(Plag);
Pfix=Pfix(end-T+1:end);
LogLikelihood=0;
% Generate data
for t=1:T
Psimul=zeros(iter,1);

% Simulate the price for one period fixed on the past in each period
for i=1:iter
pi_f=0;
for k=1:K
pi_f=pi_f+abs(alpha*(Plag(t,k+1)-Flag(t,k+1))-(Plag(t,k)-Plag(t,k+1))); % (A8) for fundamentalist rule
end

pi_c=0;
for k=1:K
    sum1=0;
    for l=1:L
       sum1=sum1+Plag(t,k+l)-Plag(t,k+l+1);
    end
    pi_c=pi_c+abs(beta*sum1-Plag(t,k)+Plag(t,k+1)); % (A8) for chartist rule
end

W=(1+exp(gamma_param*((pi_f-pi_c)./(pi_f+pi_c)))).^-1; % (A7)
sum2=0;
for l=1:L
    sum2=sum2+Plag(t,l)-Plag(t,l+1); % (A6) with chartist expectation (without beta)
end

Psimul(i)=c+d*Plag(t,1)+Plag(t,1)+W.*(alpha*(Plag(t,1)-Flag(t,1)))+(1-W).*(beta*sum2)+sigma*randn(1); % (A10)

end

bw=C*1.059*std(Psimul)*iter^-0.2;

LogLikelihood=LogLikelihood+log(KDE(Psimul,Pfix(t),bw)); % Estimate the log Likelihood as sums of conditional log Likelihoods estimated with parametric pdf

end

