% Same as function Housing but returns the W as a vector
function W = getW(theta,Plag,Flag,K,L)
c=theta(1);
d=theta(2);
alpha=theta(3);
beta=theta(4);
gamma=theta(5);
T=length(Plag);

pi_f=zeros(T,1);
for k=1:K
pi_f=pi_f+abs(alpha*(Plag(:,k+1)-Flag(:,k+1))-(Plag(:,k)-Plag(:,k+1)));
end

pi_c=zeros(T,1);
for k=1:K
    sum1=zeros(T,1);
    for l=1:L
       sum1=sum1+Plag(:,k+l)-Plag(:,k+l+1);
    end
    pi_c=pi_c+abs(beta*sum1-Plag(:,k)+Plag(:,k+1));
end

W=(1+exp(gamma*((pi_f-pi_c)./(pi_f+pi_c)))).^-1;
sum2=zeros(T,1);
for l=1:L
    sum2=sum2+Plag(:,l)-Plag(:,l+1);
end

Pt=c*ones(T,1)+d*Plag(:,1)+Plag(:,1)+W.*(alpha*(Plag(:,1)-Flag(:,1)))+(1-W).*(beta*sum2);

end
