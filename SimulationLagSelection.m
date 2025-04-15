% Code performs Monte Carlo simulations calibrated around MLE and checks
% how often the correct Lag number is selected, see Chapter 5. (Warning takes a lot of
% time, about 2 hours on my laptop)
clear
% Data
data=xlsread("RENT-PRICE-RATIO.2018q2.xlsx");
lastdata=234;
rent=data(1:lastdata,2);
price=data(1:lastdata,3);
dates=data(1:lastdata,1);
CPI=xlsread("CPIAUCSL.xls");
CPI=CPI(1:lastdata);

% Seed and Monte Carlo setup
rng(666)
MC=1000;
Ks=zeros(MC,1);
Ls=zeros(MC,1);
for m=1:MC
% Simulation of time series
Ffix=log(calcfund(rent./CPI,price./CPI)); % Get actual series of F
theta=[-0.0749 0.0117 -0.1484 0.4986 1.2838]; K=1; L=3;
c=theta(1);
d=theta(2);
alpha=theta(3);
beta=theta(4);
gamma=theta(5);
Tsimul=229;
Psimul=zeros(Tsimul+10,1);
Psimul(1:10)=6.1972; % Set presample
Fsimul=[ones(10,1)*6.2440; Ffix];
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

Psimul(t)=c+d*Psimul(t-1)+Psimul(t-1)+W.*(alpha*(Psimul(t-1)-Fsimul(t-1)))+(1-W).*(beta*sum2)+0.0063*randn(1); % (A10) with noise
end

P=Psimul(11:end);
F=Fsimul(11:end);

Likelihoods=zeros(12,12);
for k=1:12
    for l=1:12
K=k;
L=l;
initial=theta;
[thetahat exitflag Likelihood] = MLE(P,F,K,L,initial,1);
Likelihoods(k,l)=Likelihood;
    end
end
[MaxPerColumn, RowIndices] = max(Likelihoods);
[MaxValue, ColumnIndex] = max(MaxPerColumn);
RowIndex = RowIndices(ColumnIndex);
Ks(m)=RowIndex; % Save selected K
Ls(m)=ColumnIndex; % Save selected L
m 
end

% 897 correct