function [ c,ceq ] = nonlinear_outer( x,g )
%This function computes nonlinear inequalities of the outer bound optimization
%problem.

%g=[correlation btwn X1&X2, correlation btwn X2&X3, correlation btwn X1&X3]

%We denote ith encoder's message by Mi

%Decision variables of the optimization problem on outer bound are the
%following

%I(X1^n;M1)/n=x(1)
%I(X1^n;M2)/n=x(2)
%I(X1^n;M3)/n=x(3)
%I(X2^n;M1)/n=x(4)
%I(X2^n;M2)/n=x(5)
%I(X2^n;M3)/n=x(6)
%I(X3^n;M1)/n=x(7)
%I(X3^n;M2)/n=x(8)
%I(X3^n;M3)/n=x(9)
%I(X1^n;M2,M3)/n=x(10)
%I(X3^n;M1,M2)/n=x(11)
%I(M1;M2)/n=x(12)
%I(M1;M3)/n=x(13)
%I(M2;M3)/n=x(14)
%I(M1;M2, M3)/n=x(15)

%One may also expect to see I(M3;M1, M2)/n=x(16) because of the symmetry of
%the problem. Since we have I(M3;M1, M2)=I(M1;M2, M3)+I(M2;M3) -I(M1;M2), 
%we choose not to include it among decision variables of the optimization 
%problem. Whenever it appears in teh inequalities, we will simply use 
%I(M3;M1, M2)/n=x(15)+x(14)-x(12)

%We have no nonlinear equality constraint in the optimization problem
ceq=[];

%Consider the Markov chain M1-X1^n-X2^n-M2 and the associated Courtade
%inequality. Taking log of both sides, then gathering terms on one side
%yields 
%1/n[I(X1^n;M2)+I(X2^n;M1)-I(M1;M2)]+
%0.5*log[1-g1^2+g1^2 * 2^(-2/n^(I(X1^n;M1)+I(X2^n;M2)-I(M1;M2)))]<=0
c(1)=x(2)+x(4)-x(12)+0.5 * log2( 1-g(1)^2 + g(1)^2 *2^(-2*(x(1)+x(5)-x(12))));

%Proceeding similarly, Courtade's inequality on M1-X1^n-X3^n-M3 yields
c(2)=x(3)+x(7)-x(13)+0.5 * log2( 1-g(3)^2 + g(3)^2 *2^(-2*(x(1)+x(9)-x(13))));

%Courtade's inequality on M1-X1^n-X2^n-M3 yields 
c(3)=x(3)+x(4)-x(13)+0.5 * log2( 1-g(1)^2 + g(1)^2 *2^(-2*(x(1)+x(6)-x(13))));

%Courtade's inequality on M1-X1^n-X2^n-(M2,M3) yields 
c(4)=x(10)+x(4)-x(15)+0.5 * log2( 1-g(1)^2 + g(1)^2 *2^(-2*(x(1)+x(5)+x(6)-x(14)-x(15))));

%Courtade's inequality on M1-X1^n-X2^n-X3^n yields 
c(5)=-0.5*log2(1-g(3)^2)+x(4)-x(7)+0.5*log2(1-g(1)^2+g(1)^2*(1-g(2)^2)*2^(-2*(x(1)-x(7))));

%Courtade's inequality on X1^n-X2^n-X3^n-M3 yields 
c(6)=-0.5*log2(1-g(3)^2)+x(6)-x(3)+0.5*log2(1-g(2)^2 + g(2)^2*(1-g(1)^2)*2^(-2*(x(9)-x(3))));

%Courtade's inequality on M2-X2^n-X3^n-M3 yields 
c(7)=x(6)+x(8)-x(14)+0.5 * log2( 1-g(2)^2 + g(2)^2 *2^(-2*(x(5)+x(9)-x(14))));

%Courtade's inequality on M1-X2^n-X3^n-M3 yields 
c(8)=x(6)+x(7)-x(13)+0.5 * log2( 1-g(2)^2 + g(2)^2 *2^(-2*(x(4)+x(9)-x(13))));

%Courtade's inequality on (M1,M2)-X2^n-X3^n-M3 after expressing I(M3;M1,M2)
%as I(M1;M2,M3)+I(M2;M3)-I(M1;M2) gives
c(9)=x(6)+x(11)-x(15)-x(14)+x(12)+0.5*log2(1-g(2)^2+g(2)^2 *2^(-2*(x(4)+x(5)+x(9)-x(15)-x(14))));

%Consider the Markov chain X2^n-X1^n-M1 and apply Oohama's inequality.
%Taking logs of both sides and gathering terms on one side gives
% I(X2^n;M1)+ 0.5 * log2(1-g(1)^2+g(1)^2 * 2^(-2/n * I(X1^n;M1)))
c(10)=x(4)+0.5 * log2(1-g(1)^2+g(1)^2 * 2^(-2*x(1)));

%Oohama's inequality on %X3^n-X2^n-M1
c(11)=x(7)+0.5 * log2(1-g(2)^2+g(2)^2 * 2^(-2*x(4)));

%Oohama's inequality on X3^n-X1^n-M1
c(12)=x(7)+0.5 * log2(1-g(3)^2+g(3)^2 * 2^(-2*x(1)));

%Oohama's inequality on X1^n-X2^n-M2
c(13)=x(2)+0.5 * log2(1-g(1)^2+g(1)^2 * 2^(-2*x(5)));

%Oohama's inequality on X3^n-X2^n-M2
c(14)=x(8)+0.5 * log2(1-g(2)^2+g(2)^2 * 2^(-2*x(5)));

%Oohama's inequality on X1^n-X3^n-M3
c(15)=x(3)+0.5 * log2(1-g(3)^2+g(3)^2 * 2^(-2*x(9)));

%Oohama's inequality on X2^n-X3^n-M3
c(16)=x(6)+0.5 * log2(1-g(2)^2+g(2)^2 * 2^(-2*x(9)));

%Oohama's inequality on X1^n-X2^n-M3
c(17)=x(3)+0.5 * log2(1-g(1)^2+g(1)^2 * 2^(-2*x(6)));

%Oohama's inequality on X1^n-X2^n-(M2,M3)
c(18)=x(10)+0.5 * log2(1-g(1)^2+g(1)^2 * 2^(-2*(x(5)+x(6)-x(14))));

%Oohama's inequality on X3^n-X2^n-(M1,M2)
c(19)=x(11)+0.5 * log2(1-g(2)^2+g(2)^2 * 2^(-2*(x(5)+x(4)-x(12))));


end