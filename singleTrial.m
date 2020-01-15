function [inner_min_so_far,outer_value,inner_optimizers_so_far, outer_optimizers,...
    percentage_loss]=singleTrial(mu, gammas, rhos, ...
    use_inital_inner,use_initial_outer,initial_inner, initial_outer)

%Every decision variable is a mutual information so all of them are
%lower bounded by 0.
A=-eye(15);

%We will add 13 more linear constraints
A=[A ; zeros(13,15)];
% gamma1<=I(X1; M1,M2,M3)/n=[I(X1;M1)+I(X1;M2,M3)-I(M1;M2,M3)]/n
A(16, 1)=-1;  A(16,10)=-1;  A(16, 15)=1;
% gamma2<=I(X2; M1,M2,M3)/n=[I(X2^n;M1)+I(X2^n;M2)+I(X2^n;M3)-I(M1;M2,M3)-I(M2;M3)]/n
A(17, 4)=-1;  A(17, 5)=-1;  A(17,6)=-1;A(17, 14)=1; A(17, 15)=1;
%gamma3<=[I(X3^n;M3)+I(X3^n;M1,M2)-I(M1;M2,M3)-I(M2;M3)+I(M1;M2)]/n
A(18, 9)=-1;  A(18,11)=-1;  A(18, 15)=1; A(18, 14)=1; A(18, 12)=-1;

%I(M3;M1,M2)/n >= I(M3;M1)/n
%But express I(M3;M1,M2) as I(M1;M2,M3)+I(M2;M3)-I(M1;M2)
A(19, 13)=1; A(19, 12)=1;A(19, 14)=-1; A(19, 15)=-1;

%I(X3^n;M1,M2)/n >= I(X3^n;M2)/n
A(20, 11)=-1; A(20, 8)=1;
%I(X3^n;M1,M2)/n >= I(X3^n;M1)/n
A(21, 11)=-1; A(21, 7)=1;

%I(X1^n;M2,M3)/n >= I(X1^n;M3)/n
A(22, 10)=-1; A(22,3)=1;
%I(X1^n;M2,M3)/n >= I(X1^n;M2)/n
A(23, 10)=-1; A(23,2)=1;

%I(M1;M2,M3)/n >= I(M1;M3)/n
A(24, 15)=-1; A(24,13)=1;
%I(M1;M2,M3)/n >= I(M1;M2)/n
A(25, 15)=-1; A(25,12)=1;

%Data-processing inequality arising from M1-X1^n-(M2,M3)
A(26, 1)=-1;  A(26,15)=1;

%Data-processing inequalities arising from M3-X3^n-(M1,M2)
%Again we express I(M3;M1,M2) as I(M1;M2,M3)+I(M2;M3)-I(M1;M2)
A(27,9)=-1; A(27,15)=1; A(27,14)=1; A(27,12)=-1;

%Data-processing inequality from M2-X2^n-(M1,M3) to force I(X2;U2|U1U3)=I(X2;U2)-I(U2;U1U3)=
%I(X2;U2)+I(U1;U3)-I(U2;U3)-I(U1;U2U3)to be nonnegative
A(28,5)=-1; A(28,15)=1; A(28,14)=1; A(28,13)=-1;

%Ax+b<=0
b=zeros(28,1);
%Add distortions to linear constraint of the form Ax+b<=0
b(16:18)=-gammas;

%Since the outer bound optimization problem is convex, we can start
%anywhere. Assign the point that the user prvoded as the initial point 
%if flag is raised. Otherwise assign all zero point.
if use_initial_outer==1
    initial_value_outer=initial_outer;
else
    initial_value_outer=zeros(15,1);
end

%Inner bound is expressed in terms of noise variances and noise variances
%are nonnegative so we have the following lower bound on three noise
%variances
lb=[0,0,0];


%We start with relatively tight tolerances for constraints and optimality.
%But in some choice of distortion constraints, weights and correlation
%coefficients, fmincon will give warning of almost singular matrix. When
%that happens, we can't completely trust the result of fmincon. Then we
%will gradually use more relaxed constraint and optimality tolerances.

%Relaxed tolerances may result in outer bound being slightly(around the first 
%optimality tolerance that does not lead to singular matrix or nearly singular matrix warning)
%greater than Berger-Tung inner bound from time to time, which is
%theoretically impossible. In those cases, one can deduce that they
%actually coincide. Since we output the values of the inner and outer bound,
%you can inspect those cases for yourself.
options1=optimoptions(@fmincon,'Display','off','Algorithm','interior-point', ...
    'ConstraintTolerance',1e-7,'MaxFunctionEvaluations', 30000, ...
    'MaxIterations',10000 ,'OptimalityTolerance',10^-7,'StepTolerance',10^-10);


warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');



%find the contrapolymatroid "mode" for the Berger-Tung bound
if mu(1)==mu(2) && mu(2)==mu(3)
    decomposition=0;
else
    if mu(1)>=mu(2)
        if mu(2)>=mu(3)
            decomposition=1;
        else
            decomposition=3;
        end
    else
        decomposition=2;
    end
end
  
%Express the objective of the outer bound's optimization problem in terms of 
%x vector based on the ordering of the weights mu
if decomposition==0
    outer=@(x)mu(1)*(x(1)+x(5)+x(9)-x(14)-x(15));
elseif decomposition==1
    outer=@(x)mu(1)*(x(1)-x(15)) +mu(2)*(x(5)-x(14))+mu(3)*x(9);
elseif decomposition==2
    outer=@(x)mu(1)*(x(1)-x(13)) +mu(2)*(x(5)+x(13)-x(14)-x(15))+mu(3)*x(9);
elseif decomposition==3
    outer=@(x)mu(1)*(x(1)-x(15)) +mu(2)*x(5)+mu(3)*(x(9)-x(14));
end

%Compute the outer bound value and outer bound optimizers
while true
    %Solve the convex optimization problem for the outer bound with the set
    %of the tolerances we specified above
    warning('');
    [outer_optimizers,outer_value,exitflag]...
        =fmincon(outer,initial_value_outer,A,b,[],[],[],[],@(x)nonlinear_outer(x,rhos),options1);
    
    %If there is a warning of an almost singular matrix or singular matrix, we should catch
    %it and we can't trust this outer value and outer optimizers
    [msglast, msgidlast] = lastwarn;
    if ~strcmp(msgidlast,'MATLAB:singularMatrix') && ...
            ~strcmp(msgidlast,'MATLAB:nearlySingularMatrix')&& ...
            (exitflag==1 || exitflag==2)
        %Now that there is no warning, and either local minimum has been
        %found (exitflag==1) or local minimum is possible(exitflag==2),
        % there is no reason to suspect the minima
        
        %Assign tolerances to the our default values for next iterations
        options1.ConstraintTolerance=10^-7;
        options1.OptimalityTolerance=10^-7;
        %Clear the warning
        warning('');
        break
    elseif strcmp(msgidlast,'MATLAB:singularMatrix') || ...
            strcmp(msgidlast,'MATLAB:nearlySingularMatrix')
        %Gradually relax optimality and constraint tolerances to avoid
        %warnings
        options1.OptimalityTolerance=options1.OptimalityTolerance*2;
        options1.ConstraintTolerance=options1.ConstraintTolerance*2;
        warning('');
        
    %If it enters here double the 'MaxIterations' and 'MaxFunctionEvaluations'
    else
        options1.MaxIterations=options1.MaxIterations*2;
        options1.MaxFunctionEvaluations=options1.MaxFunctionEvaluations*2;
    end
    
    %If the tolerances are too muchs relaxed, the program should warn us
    %This never happened in our simulations
    if options1.ConstraintTolerance>10^-3 || options1.OptimalityTolerance>10^-3
        disp("The tolerances have been relaxed too much to avoid singular matrix warning");
    end
end

%We are done with the outer bound. Next, we will try a number of initial
%points for the inner bound's optimization problem and pick the minimum of 
%the inner values. This is simply because the inner bound's optimization
%problem is not convex and sometimes quite sensitive to the initial point.

%Derive an initial noise variance wi by assuming I(Xi^n;Mi)=I(Xi;Ui)
%where Ui=Xi+Wi with Wi's variance is denoted by wi
outer_optimizers_noises=[((2^(2*outer_optimizers(1))-1)^-1),...
    ((2^(2*outer_optimizers(5))-1)^-1),...
    ((2^(2*outer_optimizers(9))-1)^-1)];

%For some large values of distortion constraint, which means very high
%resolution, we need very small step tolerance to find local minima because noise
%variances will be very small
step_tol=min(min(outer_optimizers_noises)*10^-3,10^-17);

inner_min_so_far=Inf;
coefficient=1;
options2=optimoptions(@fmincon,'Display','off','Algorithm','interior-point', ...
    'ConstraintTolerance',1e-9,'MaxFunctionEvaluations', 15000, ...
    'MaxIterations',5000 ,'OptimalityTolerance',1e-9,'StepTolerance',step_tol);


%Since the inner bound problem is non-convex, we try different initial
%noise variance tuples starting with noise tuple (w1,w2,w3) derived from
%outer optimizers, with less noise in each iteration.
%And we stop either after 20 different iteration or inner bound is less
%than 1.02 times the outer value
for m=1:20
    % We use the noise variance triple derived from outer optimizers
    %(w1,w2,w3) as the initial point
    [inner_optimizers,inner_value,exitflag]=...
        fmincon(@(w)inner_func(w,rhos,mu,decomposition),...
        [((2^(2*outer_optimizers(1))-1)^-1)/coefficient, ...
        ((2^(2*outer_optimizers(5))-1)^-1)/coefficient,...
        ((2^(2*outer_optimizers(9))-1)^-1)/coefficient],...
        [],[],[],[],lb,[],@(w)nonlinear_inner(w,gammas,rhos),options2);
    %If there was a warning, clear it
    warning('')
    %Just in case there may be a warning, check the distortion
    %constraints again and recompute the objective with the supposedly
    %optimal inner_optimizers
    [cond,eq]=nonlinear_inner(inner_optimizers,gammas,rhos);
    inner_value=inner_func(inner_optimizers, rhos, mu, decomposition);
    %If the solver has not failed and the constraints are satisfied,
    %and the new objective is less than the minimums of objectives so
    %far, update the minimum inner bound
    if all(cond<=0+1e-9)&&(exitflag==1 || exitflag==2) && inner_value < inner_min_so_far
        inner_min_so_far=inner_value;
        inner_optimizers_so_far=inner_optimizers;
        %To speed up the code, we break the loop if the inner bound found
        %in this iteration is less than 2% higher than the outer bound
        if inner_value< outer_value*1.02
            break
        end
    end
    %decrease the noise variances for the next iteration
    coefficient=coefficient*1.1;
end

%Provide a starting point to the non-convex inner bound that would satisfy
%the distortion constraints as if the problem were composed of three
%point-to-point problems. For numerical stability, add a small positive
%offset to the distortion constraint in case a target distortion is 0
[inner_optimizers,inner_value,exitflag]=...
    fmincon(@(w)inner_func(w,rhos,mu,decomposition),...
    [((2^(2*gammas(1)+10^-7)-1)^-1),...
    ((2^(2*gammas(2)+10^-7)-1)^-1),...
    ((2^(2*gammas(3)+10^-7)-1)^-1)],...
    [],[],[],[],lb,[],@(w)nonlinear_inner(w,gammas,rhos),options2);
%Clear warning if there was any
warning('');
[cond,eq]=nonlinear_inner(inner_optimizers,gammas,rhos);
inner_value=inner_func(inner_optimizers, rhos, mu, decomposition);
if all(cond<=0+1e-9)&& (exitflag==1 || exitflag==2) && inner_value<inner_min_so_far
    inner_min_so_far=inner_value;
    inner_optimizers_so_far=inner_optimizers;
end

%Use the inital point that the user provided if the flag is 1
if use_inital_inner==1
    [inner_optimizers,inner_value,exitflag]=...
        fmincon(@(w)inner_func(w,rhos,mu,decomposition),initial_inner,...
        [],[],[],[],lb,[],@(w)nonlinear_inner(w,gammas,rhos),options2);
    %Clear warning if there was any
    warning('');
    [cond,eq]=nonlinear_inner(inner_optimizers,gammas,rhos);
    inner_value=inner_func(inner_optimizers, rhos, mu, decomposition);
    if all(cond<=0+1e-9)&& (exitflag==1 || exitflag==2) && inner_value<inner_min_so_far
        inner_min_so_far=inner_value;
        inner_optimizers_so_far=inner_optimizers;
    end
end

%Sometimes U2=null initial point work better than the initial points above
%We evaluate the optimization problem only if there has been too much performance loss
%so far
%Start searching with M2=null
if inner_min_so_far > outer_value*1.02
    [inner_optimizers,inner_value,exitflag]=...
        fmincon(@(w)inner_func(w,rhos,mu,decomposition),...
        [((2^(2*gammas(1)+10^-7)-1)^-1),...
        10^20,...
        ((2^(2*gammas(3)+10^-7)-1)^-1)],...
        [],[],[],[],lb,[],@(w)nonlinear_inner(w,gammas,rhos),options2);
    warning('');
    [cond,eq]=nonlinear_inner(inner_optimizers,gammas,rhos);
    inner_value=inner_func(inner_optimizers, rhos, mu, decomposition);
    if all(cond<=0+1e-9)&& (exitflag==1 || exitflag==2) && inner_value<inner_min_so_far
        inner_min_so_far=inner_value;
        inner_optimizers_so_far=inner_optimizers;
    end
end

percentage_loss=max(0,100*(inner_min_so_far-outer_value)/inner_min_so_far);
end