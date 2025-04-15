% Calcs Fundamental price with given price and rent series following
% Kouwenberg & Zwinkels (2014)
function F = calcfund(rent,P)
T=length(rent);
F=zeros(T,1);
for i=6:T
    g=mean(rent(5:i-1)./rent(1:i-4-1)-1); % Real time estimate of rent growth rate
    HP=mean(rent(1:i-1)./P(1:i-1)); % Real time estimate of expected rent yield
    F(i)=(1+g)/HP*rent(i);
end
F=F(6:end);


