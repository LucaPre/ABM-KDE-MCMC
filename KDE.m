% Kernel estimator for given bandwidth, sample, point of evaluation with
% Gaussian Kernel
function density = KDE(sample,x,bw)
K=@(u) 1/sqrt(2*pi)*exp(-0.5*u.^2);
N=length(sample);
density=1/(N*bw)*sum(K((sample-x)./bw))+1e-50; % Add small value for numerical stability (otherwise we might take logs of 0)