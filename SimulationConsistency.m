% Code performs simulations calibrated around MLE 
clear
Tsimul=230; % Choose length of simulated series (confirm that estimation gets closer as Tsimul increases)
theta=[-0.0749 0.0117 -0.1484 0.4986 1.2838]; K=1; L=3; sigma=0.0063; % Set true parameters and lag lengths (here to MLE)
c=theta(1);
d=theta(2);
alpha=theta(3);
beta=theta(4);
gamma=theta(5);
Psimul=zeros(Tsimul+10,1);
Psimul(1:10)=6.1972; % Set presample
Fsimul=ones(10+Tsimul,1)*6.2440; % Set F (here to some constant)
% Recursively generate data
for t=11:Tsimul+10
pi_f=0;
for k=1:K
pi_f=pi_f+abs(alpha*(Psimul(t-(k+1))-Fsimul(t-(k+1)))-(Psimul(t-k)-Psimul(t-(k+1)))); % (A8) for fundamentalist rule
end
pi_c=0;
for k=1:K
    sum1=0;
    for l=1:L
       sum1=sum1+Psimul(t-(k+l))-Psimul(t-(k+l+1));
    end
    pi_c=pi_c+abs(beta*sum1-Psimul(t-k)+Psimul(t-(k+1))); % (A8) for chartist rule
end

W=(1+exp(gamma*((pi_f-pi_c)./(pi_f+pi_c)))).^-1; % (A7)
sum2=0;
for l=1:L
    sum2=sum2+Psimul(t-l)-Psimul(t-(l+1)); % (A6) with chartist expectation (without beta)
end

Psimul(t)=c+d*Psimul(t-1)+Psimul(t-1)+W.*(alpha*(Psimul(t-1)-Fsimul(t-1)))+(1-W).*(beta*sum2)+sigma*randn(1); % (A10) with noise
end

P=Psimul(11:end);
F=Fsimul(11:end);

% Apply estimation techniques (starting value chosen to the true value for
% simplification)
thetahat=MLE(P,F,K,L,theta,0)

S0=100000; S1=1000000;
[thetadraws acceptances h]=RWMH(S0,S1,P,F,K,L,theta); % may not work withouth enough time for high Tsimul without adjusting the function (acceptance rate will decrease as Tsimul increases, proposal distribution would need to be adjusted)

names = {'c','d','$\alpha$','$\beta$','$\gamma$'};
figure 
for i=1:5
    subplot(3,2,i)
    histogram(thetadraws(1:end,i),100,'Normalization','probability')
    title(names{i},'Interpreter','Latex','FontSize',12)
end

