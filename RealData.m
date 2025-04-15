clear
% Import data
data=xlsread("RENT-PRICE-RATIO.2018q2.xlsx");
lastdata=234;
rent=data(1:lastdata,2);
price=data(1:lastdata,3);
dates=data(1:lastdata,1);
CPI=xlsread("CPIAUCSL.xls");
CPI=CPI(1:lastdata);

% Calc fundamental price
F=log(calcfund(rent./CPI,price./CPI));
P=log(price(6:end)./CPI(6:end));
dates=dates(6:end);

% Figure A.1
figure 
plot(dates,P,'DisplayName', 'P_t')
hold on
plot(dates,F,'DisplayName', 'F_t')
legend('Location','best')

figure 
plot(dates,P-F,'DisplayName', 'P_t - F_t')
legend('Location','best')
yline(0,'k--','HandleVisibility', 'off')


% MLE and Lag selection (chapter 5)
Likelihoods=zeros(12,12);
for k=1:12
    for l=1:12
K=k;
L=l;
thetahat=MLEres(P,F,K,L,[0 0 0 0 0]); % Estimate restricted model to get starting values
initial=thetahat;
initial(5)=2; % Set starting value for gamma to something reasonable
[thetahat, exitflag, Likelihood] = MLE(P,F,K,L,initial,1);
Likelihoods(k,l)=Likelihood;
    end
end

% Choose values K and L that maximize the Likelihood (here Likelihoods(1,3))
K=1;
L=3;

thetahat=MLEres(P,F,K,L,[0 0 0 0 0]); 
initial=thetahat;
initial(5)=2; % 
[thetahat, exitflag, Likelihood, sigmaMLE] = MLE(P,F,K,L,initial,0); % get the final MLE estimates
 
% Code to generate baseline results (Sampling of Wchains may lead to memory
% problems for a large chain depending on the working station, decrease S1 if thats the case)
rng(666)
seeds=randi([1 2^32-1],4,1);
S0=10000; % Burn in
S1=1000000; % Number of draws
thetachains=zeros(S1,5,4);
hchains=zeros(S1,1,4);
Wchains=zeros(224,S1,4);
parfor i=1:4
rng(seeds(i))
[theta,acceptances,h,Ws]=RWMH(S0,S1,P,F,K,L,[-0.0749 0.0117 -0.1484 0.4986 1.2838]);
mean(acceptances)
thetachains(:,:,i)=theta(S0+1:end,:);
hchains(:,:,i)=h(S0+1:end,:);
Wchains(:,:,i)=Ws(:,S0+1:end);
end
theta=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)]; % Mix all chains
% Running time (at my Laptop) 408 seconds, acceptance rates for chains between
% 24.46%-24.52%

% Posterior statistics Table 1
thetamean=mean(theta);
thetamode=histmode(theta);
thetamedian=median(theta);
thetastd=std(theta);

% Figure A.2
names = {'c','d','$\alpha$','$\beta$','$\gamma$'};
figure 
for j=1:4
for i=1:5
    subplot(2,3,i)
    plot(thetachains(:,i,j))
    title(names{i},'Interpreter','Latex','FontSize',12)
    xlabel('S')
    hold on
end
end

% Convergence tests for A.4 (Warning: takes a lot of time, 25 minutes for me)
CD=zeros(5,4,1);
parfor j=1:4
    for i=1:5
        CD(i,j)=(mean(thetachains(1:100000,i,j))-mean(thetachains(500001:end,i,j)))/sqrt(NeweyWest(thetachains(1:100000,i,j),20000)/100000+NeweyWest(thetachains(500001:end,i,j),20000)/500000);
    i
    end
end

% Figure 1
names = {'c','d','$\alpha$','$\beta$','$\gamma$'};
figure 
for i=1:5
    subplot(3,2,i)
    histogram(theta(1:end,i),100,'Normalization','probability')
    ax = gca;
    ax.XRuler.Exponent = 0;
    title(names{i},'Interpreter','Latex','FontSize',12)
end

% Figure 2
Ws=[Wchains(:,:,1) Wchains(:,:,2) Wchains(:,:,3) Wchains(:,:,4)];
LowerW=zeros(224,1);
MedianW=zeros(224,1);
UpperW=zeros(224,1);
for t=1:224
LowerW(t)=quantile(Ws(t,:),0.025);
MedianW(t)=quantile(Ws(t,:),0.5);
UpperW(t)=quantile(Ws(t,:),0.975);
end

figure
plot(dates(6:end),MedianW,'DisplayName', 'Median','LineWidth',1)
ylabel('W_t','Rotation',0)
hold on
plot(dates(6:end),LowerW,'LineStyle','--','Color','red','LineWidth',.05,'DisplayName', '2.5% and 97.5% quantile')
hold on
plot(dates(6:end),UpperW,'LineStyle','--','Color','red','LineWidth',.05,'HandleVisibility', 'off')
hold on
plot((min(dates):max(dates)),0.5*ones(1,57),'Color','black','HandleVisibility', 'off','LineWidth',1)
legend('Location','best')


% Simulation limit cycle
load("baseline.mat")
theta=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)]; 
rng(666)
simulations=10000;
Tsimul=500;
paths=zeros(Tsimul,simulations);
mins=zeros(simulations,1);
maxs=zeros(simulations,1);
diffs=zeros(simulations,1);
lengths=zeros(simulations,1);
for i=1:simulations
if i==1
    thetadraw=thetahat; K=1; L=3; % Use MLE for first simulation
else
ind=randi([1 length(theta)],1); % Draw randomly from theta
thetadraw=theta(ind,:); K=1; L=3;
end
c=thetadraw(1);
d=thetadraw(2);
alpha=thetadraw(3);
beta=thetadraw(4);
gamma=thetadraw(5);
Psimul=6.6029*ones(10,1)+0.0000001;
Psimul=[Psimul;zeros(Tsimul,1)];
Fsimul=6.6029*ones(Tsimul+10,1);
count=0;
for t=11:Tsimul+10
pi_f=0;
for k=1:K
pi_f=pi_f+abs(alpha*(Psimul(t-(k+1))-Fsimul(t-(k+1)))-(Psimul(t-k)-Psimul(t-(k+1))));
end
pi_c=0;
for k=1:K
    sum1=0;
    for l=1:L
       sum1=sum1+Psimul(t-(k+l))-Psimul(t-(k+l+1));
    end
    pi_c=pi_c+abs(beta*sum1-Psimul(t-k)+Psimul(t-(k+1)));
end
W=(1+exp(gamma*((pi_f-pi_c)./(pi_f+pi_c)))).^-1;
sum2=0;
for l=1:L
    sum2=sum2+Psimul(t-l)-Psimul(t-(l+1));
end
Psimul(t)=c+d*Psimul(t-1)+Psimul(t-1)+W.*(alpha*(Psimul(t-1)-Fsimul(t-1)))+(1-W).*(beta*sum2);
if count==0
if Psimul(t)<Psimul(t-1) 
    j1=t-1;
    count=count+1;
end
end
if count==1
    if Psimul(t)>Psimul(t-1) 
    j2=t-1;
    count=count+1;
    end
end
if count==2
    if Psimul(t)<Psimul(t-1) 
    j3=t-1;
    count=count+1;
    end
end
end
paths(:,i)=Psimul(11:end); 
mins(i)=min(Psimul(100:end)); % Save minimum of cycle
maxs(i)=max(Psimul(100:end)); % Save maximum of cycle
diffs(i)=maxs(i)-mins(i); % Save amplitude
lengths(i)=j3-j1; % Save cycle length
end

% Behaviour of simulated cycles (Table 2)
prctile(lengths(2:end),[2.5 50 97.5])
prctile(diffs(2:end),[2.5 50 97.5])
prctile(mins(2:end),[2.5 50 97.5])-6.6029
prctile(maxs(2:end),[2.5 50 97.5])-6.6029

% Figure 3
figure
for i=1:20
    if i==1
      plot(paths(:,i),'--','linewidth',2,'DisplayName', 'MLE')
      legend('location','best')
      hold on
    else
    plot(paths(:,i),'HandleVisibility', 'off')
    xlabel('t')
    ylabel('$P_t$','Interpreter','Latex','rotation',0)
    hold on
    end
end

%% KDE estimate of likelihood  
% Code to generate results for Table 3 (Warning: Takes a lot of time)
tic
rng(666)
seeds=randi([1 2^32-1],4,1);
S0=10000;
S1=100000;
iter=1000; % Number of simulations at each t
thetachains=zeros(S1,5,4);
parfor i=1:4
rng(seeds(i))
[theta,acceptances]=KDE_MH(S0,S1,P,F,K,L,[-0.0749 0.0117 -0.1484 0.4986 1.2838],iter,2); % Replace last input with the desired constant C equation (6)
mean(acceptances)
thetachains(:,:,i)=theta(S0+1:end,:);
end
toc
theta=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
% Running time 20225 seconds (can vary by a couple thousand seconds)

thetamean=mean(theta);
thetastd=std(theta);

% KLD (here I loaded data collected over several days from code before)
load("baseline.mat")
theta1=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
load("C=2.mat") % Insert here data from draws with bandwidth of choice
theta2=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
rng(666)
KLD(theta1,theta2,1000)


% Figure 3
load("baseline.mat")
theta1=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
load("C=1.mat") 
theta2=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
load("C=2.mat")
theta3=[thetachains(:,:,1); thetachains(:,:,2); thetachains(:,:,3); thetachains(:,:,4)];
thetamean=mean(theta1);

names = {'c','d','$\alpha$','$\beta$','$\gamma$'};
figure
count=0;
for j=1:3
if j==1
    theta=theta1;
end
if j==2
    theta=theta2;
end
if j==3
    theta=theta3;
end
for i=1:5
    count=count+1;
    subplot(3,5,count)
    histogram(theta(1:end,i),100,'Normalization','probability')
    hold on
    plot(thetamean(i)*ones(1,5),0:0.01:0.04)
    ax = gca;
    ax.XRuler.Exponent = 0;
    if j==1
    title(names{i},'Interpreter','Latex','FontSize',12)
    end
    if i==1
    if j==1
ylabel("RWMH",'Rotation',0);
    end
    if j==2
 ylabel({'Silverman'},'Rotation',0);
    end
    if j==3
ylabel({'C=2'},'Rotation',0);
    end
    end
end
end

% BW selection A.5 (Warning: Takes a lot of time)
rng(666)
iter=1000;
S0=1000;
S1=10000;
[theta,acceptances,h]=KDE_MH(S0,S1,P,F,K,L,[-0.0749 0.0117 -0.1484 0.4986 1.2838],iter,1); % Approximate sampling from posterior with silvermans rule
Plag=getLags(P,K+L+1);
Flag=getLags(F,K+L+1);
T=length(Plag);
posttheta=@(theta,h) (-0.5*h*(P(end-T+1:end)-Housing(theta,Plag,Flag,K,L))'*(P(end-T+1:end)-Housing(theta,Plag,Flag,K,L))); % Function for true acceptance ratio
M=100; % Number of draws of posterior for averaging of acceptance ratios
seeds=randi([1 2^32-1],M,1);
bws=0.25:0.25:5; % Grid of bandwidths
Nstar=50000; % Simulation size of pseudo-true ratio
pseudotrueratios=zeros(M,1);
ratiokdes=zeros(M,length(bws));
trueratios=zeros(M,1);

% Calculate true and pseudotrue acceptance ratios over the (fixed seeds)
% draws of the posterior
for i=1:M
    rng(seeds(i))
    randnum=randi([S0 S0+S1]);
    thetaold=theta(randnum-1,:);
    thetastar=mvnrnd(thetaold,0.0075*diag([0.0162 0.0025 0.0165 0.0264 0.1951].^2));
    rng(seeds(i)) % Fix seeds before simulation 
    postthetaold= SimulLik(thetaold,h(randnum),P,F,K,L,Nstar,1); % Log-Likelihood under proposed draw
    rng(seeds(i))
    postthetaprop= SimulLik(thetastar,h(randnum),P,F,K,L,Nstar,1); % Log-Likelihood under old draw
    pseudotrueratio=exp(postthetaprop-postthetaold);
    trueratio=exp(posttheta(thetastar',h(randnum))-posttheta(thetaold',h(randnum)));
    pseudotrueratios(i)=min(pseudotrueratio,1);
    trueratios(i)=min(trueratio,1);
    i
end

% Calculate the estimated ratio for each bandwidth 
for j=1:length(bws)
MAE=zeros(M,1);
for i=1:M
    rng(seeds(i))
    randnum=randi([S0 S0+S1]);
    thetaold=theta(randnum-1,:);
    thetastar=mvnrnd(thetaold,0.0075*diag([0.0162 0.0025 0.0165 0.0264 0.1951].^2));
    C=bws(j);
    rng(seeds(i)) % Fix seeds before simulation 
    postthetaold= SimulLik(thetaold,h(randnum),P,F,K,L,iter,C); % Log-Likelihood under proposed draw
    rng(seeds(i))
    postthetaprop= SimulLik(thetastar,h(randnum),P,F,K,L,iter,C); % Log-Likelihood under old draw
    ratiokde=exp(postthetaprop-postthetaold);
    diff=(min(ratiokde,1)-min(pseudotrueratios(i),1));
    MAE(i)=diff;
    ratiokdes(i,j)=min(ratiokde,1);
end
j
end

MAE_true=mean(abs(ratiokdes-trueratios));
MAE_est=mean(abs(ratiokdes-pseudotrueratios));

% Run time of 50,000 simulations with M=100 for pseudotrue ratio around 900 seconds

% Figure A.3 (with data collected by the code above)
bws=0.25:0.25:5;
load("MAE_true.mat")
plot(bws(3:end),MAE_true(3:end),'-','DisplayName', 'True Likelihood')
xlabel("C")
ylabel("MAE",'rotation',0)
hold on 
load("MAE_500000.mat")
plot(bws(3:end),MAE_est(3:end),'--','DisplayName', 'N* = 500,000')
hold on 
load("MAE_50000.mat")
plot(bws(3:end),MAE_est(3:end),':','DisplayName', 'N* = 50,000','Color','black','LineWidth',1.5)
hold on 
load("MAE_25000.mat")
plot(bws(3:end),MAE_est(3:end),'-.','DisplayName', 'N* = 25,000')
legend('Location','best')



