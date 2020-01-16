function []=generatePlots(mu,fixed_rho_12,fixed_rho_23,fixed_gammas)

%Note that if mu(3)>mu(1), this function exchanges them since 
%X1 and X3's roles are interchangeable by the symmetry of the underlying 
%source-coding problem

%fixed_rho_12 and fixed_rho_23 will be used for the plot in which gamma values are varied
%together

%fixed gammas will be used for the plot in which rho_12 and rho_23 values are varied
%together

if mu(3)>mu(1)
    disp("Exchanging \mu_1 and \mu_3...");
    tmp=mu(3);
    mu(3)=mu(1);
    mu(1)=tmp;
end

%vary distortions gamma_1=gamma_2=gamma_3 from 0.1 to 10 with a step size
%of 0.01
gammas=[0.1:0.01:10];

%g holds correlations between X1 and X2, X2 and X3, X1 and X3
%g(3) is simply correlation between X1 and X3, which is basically
%correlation of (X1&X2) times correlation of (X2&X3)
%correlation co
rhos=[fixed_rho_12, fixed_rho_23,fixed_rho_12*fixed_rho_23 ];
%initialize the result arrays to hold inner and outer bound values
outer_bound_array=zeros(1,length(gammas));
inner_bound_array=ones(1,length(gammas))*Inf;

%initial noise vector
initial_w=[[((2^(2*gammas(1)+10^-7)-1)^-1),...
    ((2^(2*gammas(1)+10^-7)-1)^-1),...
    ((2^(2*gammas(1)+10^-7)-1)^-1)]];

for i=1:length(gammas)
    [inner_bound_array(i),outer_bound_array(i),inner_optimizers, ...
        outer_optimizers,percentage_loss]=singleTrial(mu, ...
        [gammas(i), gammas(i),gammas(i)], rhos, 1,0,initial_w, []);
    
%In the next iteration, we might as well use this iteration's inner
%optimizers as one of the starting points given that distortions are increasing by
%only step_size in the next iteration. singleTrial program tries many
%initial points for inner bound's non-convex optimization problem so there
%is no downside except increased running time
    initial_w=inner_optimizers;
end
figure()
plot(gammas,outer_bound_array,'b','LineWidth',1);
hold on;
plot(gammas,inner_bound_array,'r--','LineWidth',2);
xlabel({'\gamma_1=\gamma_2=\gamma_3'})
s=sprintf("%0.2f R_1 + %0.2f R_2 + %0.2f R_3", mu(1),mu(2),mu(3));
ylabel(s)
legend('Outer Bound','Inner Bound')


disp("Maximum performance loss as percentage of the inner bound ")
disp(100*max((inner_bound_array-outer_bound_array)./inner_bound_array));
percent_diff=max(zeros(1,length(gammas)), ...
    100*(inner_bound_array-outer_bound_array)./inner_bound_array );
disp("Mean percentage performance loss");
disp(mean(percent_diff));
disp("Median percentage performance loss");
disp(median(percent_diff));

%now we vary correlations rho_12=rho_23 from 0 to 0.99
rho=[0.00:0.001:0.99];
%initialize the result arrays to hold inner and outer bound values
outer_bound_array=zeros(1,length(rho));
inner_bound_array=ones(1,length(rho))*Inf;
initial_w=[[((2^(2*fixed_gammas(1)+10^-7)-1)^-1),...
    ((2^(2*fixed_gammas(2)+10^-7)-1)^-1),...
    ((2^(2*fixed_gammas(3)+10^-7)-1)^-1)]];

for i=1:length(rho)
    [inner_bound_array(i),outer_bound_array(i),inner_optimizers, ...
        outer_optimizers,percentage_loss]=singleTrial(mu, ...
        fixed_gammas, [rho(i), rho(i), rho(i)^2], 1,0,initial_w, []);
    %In the next iteration, we might as well use this iteration's inner
    %optimizers as one of the starting points given that distortions are increasing by
    %only step_size in the next iteration. singleTrial program tries many
    %initial points for inner bound's non-convex optimization problem so there
    %is no downside except increased running time
    initial_w=inner_optimizers;
    if inner_bound_array(i)==Inf
        dis("It looks like the solver for all of inner bound optimization problems failed!")
        keyboard
    end
end

figure()
plot(rho,outer_bound_array,'b','LineWidth',1);
hold on;
plot(rho,inner_bound_array,'r--','LineWidth',2);
xlabel({'\rho_{12}=\rho_{23}'})
s=sprintf("%0.2f R_1 + %0.2f R_2 + %0.2f R_3", mu(1),mu(2),mu(3));
ylabel(s)
legend('Outer Bound','Inner Bound')
disp("Maximum performance loss as percentage of the inner bound ")
disp(100*max((inner_bound_array-outer_bound_array)./inner_bound_array));
percent_diff=max(zeros(1,length(rho)), ...
    100*(inner_bound_array-outer_bound_array)./inner_bound_array );
disp("Mean percentage performance loss");
disp(mean(percent_diff));
disp("Median percentage performance loss");
disp(median(percent_diff));
end