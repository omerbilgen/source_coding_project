function [c,ceq] = nonlinear_inner(w,gammas,g)

%This function computes nonlinear constraints gamma_i-I(X_i; U1, U2, U3) given noise
%variances w and correlation coefficients g for the inner bound's
%optimization problem

%There is no nonlinear constraint
ceq=[];

%This expressions occur frequently in c(i)'s so it is preferable to compute
%them earlier.
d=g(1)^2 *(1+w(3)-g(2)^2+w(2)*g(2)^2)/((1+w(2))*(1+w(3))-g(2)^2);
d2=g(2)^2 *(1+w(1)-g(1)^2+w(2)* g(1)^2) / ((1+w(1))*(1+w(2))-g(1)^2);

c(1)=gammas(1)-0.5*log2((1+w(1)-d)/(w(1)*(1-d)));
c(2)=gammas(2)-0.5*log2((1+w(1)-d)*((1+w(2))*(1+w(3))-g(2)^2)/((1+w(1)-g(1)^2)*(1+w(3)-g(2)^2)*w(2)));
c(3)=gammas(3)-0.5*log2((1+w(3)-d2)/(w(3)*(1-d2)));


end

