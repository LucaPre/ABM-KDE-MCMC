% Computes MLE numerically for some initial value
function [thetahat, exitflag, Likelihood, sigmaMLE] = MLE(P,F,K,L,initial,lagsearch)
if lagsearch==0
Plag=getLags(P,K+L+1);
Flag=getLags(F,K+L+1);
T=length(Plag);
P=P(end-T+1:end);
else % drop observations when searching K and L so that the sample size stays constant for all pairs of lags
Plag=getLags(P,25);
Flag=getLags(F,25);
T=length(Plag);
P=P(end-T+1:end);
end

x0=initial; % Set starting values for optimization
opts = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunctionEvaluations',1000,'Display','off');
ms = MultiStart('UseParallel',false,'Display','off','XTolerance',0.001,'FunctionTolerance',0.001);
fixedFunction = @(x) LogLik(x,P,Plag,Flag,K,L); % Set function that computes the Log Likelihood for values of theta
problem = createOptimProblem('fmincon','x0',x0,'objective',fixedFunction,'options',opts);
[thetahat fval exitflag] = run(ms,problem,1); % Adjust last input to test several local minima (didn't change anything here and increases computation time)

sigmaMLE=sqrt((P-Housing(thetahat,Plag,Flag,K,L))'*(P-Housing(thetahat,Plag,Flag,K,L))/T); % Calculate MLE of error variance
hMLE=sigmaMLE^-2;
Likelihood=-0.5*hMLE*fval+0.5*T*log(hMLE); % Return the Log-Likelihood 


