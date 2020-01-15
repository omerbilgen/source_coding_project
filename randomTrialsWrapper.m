function [results]=randomTrialsWrapper()

%Get the number of trials from the user
numberOfTrials=input('Please enter the number of trials: ');

%We will store the randomly chosen weights, distortions, correlation
%coefficients and the resulting inner and outer bounds and the performance
%loss in the following matrix
results=zeros(numberOfTrials,30);


%Because of tolerances, outer bound sometimes end up being larger than
%inner bound, which is theoretically impossible. We measure this
%difference as absolute and percentage and keep track of the maximum of
%such differences
max_outer_minus_inner_absolute=0;
max_outer_minus_inner_percentage=0;


for i=1:length(results(:,1))
    
    %generate a uniformly random row-vector of 8 elements on (0,1)
    unif_vec=rand(1,8);
    
    %generate weights of encoders uniformly on (1,5)
    results(i,1:3)=4*unif_vec(1:3)+1;
    
    %generate distortions uniformly between 0 and 10
    results(i,4:6)=unif_vec(4:6)*10;
    
    %generate correlation coefficients between 0 and 1
    results(i,7)=unif_vec(7);
    results(i,8)=unif_vec(8);
    %because of Markov chain
    results(i,9)=results(i,7)*results(i,8);
    
    %if mu_3>=mu_1, since the cases are symmetric, replace them
    if results(i,3)>=results(i,1)
        tmp=results(i,3);
        results(i,3)=results(i,1);
        results(i,1)=tmp;
    end
    
    
   
    %Evaluate the outer and inner bounds for this particular choice of
    %weights, distortions and correlation coefficients
    [inner_value,outer_value,inner_optimizers, outer_optimizers,...
        percentage_loss]=singleTrial(results(i,1:3),results(i,4:6) ,...
        results(i,7:9), 0,0,[],[]);
    
    % Update the maximum absolute computational error, which is positive
    %outer_value-inner_value difference. This difference cannot be possible
    %theoretically but possible here, due to tolerances of fmincon. But we
    %will report this value at the very end, which will turn out to be very small indeed.
    if outer_value-inner_value>max_outer_minus_inner_absolute
        max_outer_minus_inner_absolute=outer_value-inner_value;
    end
    
    %Update the maximum percentage error of positive outer_value-inner_value
    %, which cannot exist theoretically but possible here, due to tolerances in
    %optimization problems' settings. But we will report at the very end that it is very
    %small like 0.0001% of inner_value
    if (outer_value-inner_value)/inner_value>max_outer_minus_inner_percentage
        max_outer_minus_inner_percentage=(outer_value-inner_value)/inner_value *100;
    end
    
    %Store the performance loss as a percentage
    results(i,10)=percentage_loss;
    %store the inner and outer bounds as well as corresponding optimizers
    results(i,11)=inner_value;
    results(i,12)=outer_value;
    results(i,13:15)=inner_optimizers;
    results(i,16:30)=outer_optimizers;    
end
%Reporting results
disp("Displaying the worst case performance loss as percentage of the inner bound");
disp(max(results(:,10)));
disp("Displaying the average percentage performance loss");
disp(mean(results(:,10)));
disp("Displaying the median percentage performance loss");
disp(median(results(:,10)));
disp("The maximum of outer_value-inner_value, which is supposed to nonpositive, is :")
disp(max_outer_minus_inner_absolute);
disp("The maximum of outer_value-inner_value as a percentage of inner_value, which is supposed to nonpositive, is :");
disp(max_outer_minus_inner_percentage);

end
