%Replicate Figure 2 Panel A from Kouwenberg & Zwinkels to ensure the model
%was implemented correctly
Tsimul=200;
theta=[-0.2307 0.0236 -0.6329 0.3032 2.1818]; K=1; L=4; 
c=theta(1);
d=theta(2);
alpha=theta(3);
beta=theta(4);
gamma=theta(5);
Psimul=10*ones(10,1)+0.0000001;
Psimul=[Psimul;zeros(Tsimul,1)];
Fsimul=10*ones(Tsimul+10,1);
Ws=zeros(Tsimul+10,1);
count=0;
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

W=(1+exp(gamma*((pi_f-pi_c)./(pi_f+pi_c)))).^-1;
Ws(t)=W;
sum2=0;
for l=1:L
    sum2=sum2+Psimul(t-l)-Psimul(t-(l+1)); % (A6) with chartist expectation (without beta)
end

Psimul(t)=c+d*Psimul(t-1)+Psimul(t-1)+W.*(alpha*(Psimul(t-1)-Fsimul(t-1)))+(1-W).*(beta*sum2); % (A10)
end

figure 
plot(Psimul(11:end),'DisplayName', 'P_t')
ylim([9.95 10.5])
hold on
plot(Fsimul(11:end),'DisplayName', 'F_t')
hold on
yyaxis right
plot(Ws(11:end),'DisplayName', 'W_t')
ylim([-0.8 0.8])
legend('Location','best')


