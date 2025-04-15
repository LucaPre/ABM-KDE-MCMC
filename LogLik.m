% Returns the negative of the Log-Likelihood (ignoring constants when
% maximizing with respect to theta)
function L = LogLik(theta,P,Plag,Flag,K,L)
L=-(P-Housing(theta,Plag,Flag,K,L))'*(P-Housing(theta,Plag,Flag,K,L));
L=-L;
end