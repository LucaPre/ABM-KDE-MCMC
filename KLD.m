% Estimates KLD based MC draws and KDE calculated with samples from two
% pdfs (theta1 is the sample for the reference pdf, note that KLD is not
% symmetric)
function KL = KLD(theta1,theta2,MC)
d=width(theta1);
n1=length(theta1);
n2=length(theta2);
H1=0.94*cov(theta1)^0.5*n1^-(1/(d+4)); % Rule of thumb bandwidth matrix proportional to covariance
H2=0.94*cov(theta2)^0.5*n2^-(1/(d+4));
K=@(u) 1/sqrt((2*pi)^d)*exp(-0.5*sum(u.^2)'); 
seeds=randi([1 2^32-1],MC,1);
ratios=zeros(MC,1);
parfor i=1:MC
    rng(seeds(i))
    ind=randi([1 n1],1);
    x=theta1(ind,:); % Random draw from the posterior 
    p=1/(n1*det(H1))*sum(K(H1^-1*(theta1-x)')); % Multivariate KDE evaluated at random draw
    q=1/(n2*det(H2))*sum(K(H2^-1*(theta2-x)'));
    ratios(i)=log(p/q); % Log ratio of pdfs
    if mod(i,10)==0
        i
    end
end
KL=mean(ratios); % Mean of log-ratios


