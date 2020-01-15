function [inner_objective]=inner_func(w,g,mu, decomposition)

%This function evaluates inner bound objective, which is mu1*R1+mu2*R2+mu*R3 
%where Ri denotes ith encoder's rate in Berger-Tung scheme and mui's can 
%be interpreted as the cost of sending information from that encoder

%g(1) is the correlation coefficient between X1 and X2 
%g(2) is the correlation coefficient between X2 and X3
%g(3) is the correlation coefficient between X1 and X3 and it is equal to
%g(1)*g(2) due to the Markov chain X1-X2-X3

%The variable decomposition defines the decoding order of messages because of
%contrapolymatroid structure of Berger-Tung achievable coding scheme

%decomposition==0 means mu(1)=mu(2)=mu(3)
%decomposition==1 means mu(1)>=mu(2)>=mu(3)
%decomposition==2 means mu(2)>=mu(1)>=mu(3)
%decomposition==3 means mu(1)>=mu(3)>=mu(2)


if w(1)<Inf
    if w(2)<Inf
        %None of the auxiliary random variables are null
        if w(3)<Inf
            % the order of decoding does not matter
            if decomposition ==0
                inner_objective=mu(1)*0.5 * log2((1+w(1)-g(1)^2 *...
                    (1+w(3)-g(2)^2+w(2)*g(2)^2)/((1+w(2))*(1+w(3))-g(2)^2))*...
                    ((1+w(2))*(1+w(3))-g(2)^2)/(w(1)*w(2)*w(3)));
                
            %decode U3, then U2 and then U1
            elseif decomposition==1
                inner_objective=0.5 * mu(1)*log2((1+w(1)-(g(1)^2 ...
                    *(1+w(3)-g(2)^2+w(2)*g(2)^2)/((1+w(2))*(1+w(3))-g(2)^2)))/(w(1)))...
                    +0.5 * mu(2)*log2(( (1+w(2)) * (1+w(3)) - g(2)^2 )/(w(2)* (1+w(3))))...
                    +0.5 * mu(3)*log2((1+w(3))/w(3));
                
            %decode U3 first, then U1, then U2
            elseif decomposition==2
                inner_objective= 0.5 * mu(2)*log2(((1+w(2))*(1+w(3))-g(2)^2)/...
                    (w(2)*((1+w(1))*(1+w(3))-g(3)^2)))+0.5*mu(2)*log2...
                    ((1+w(1)-(g(1)^2 *(1+w(3)-g(2)^2+w(2)*g(2)^2)/((1+w(2))*(1+w(3))-g(2)^2))))...
                     +0.5 * mu(1)*log2(( (1+w(1)) * (1+w(3)) - g(3)^2 )/(w(1)* (1+w(3))))...
                     +0.5 * mu(3)*log2((1+w(3))/w(3));
            
            %decode U2 first, then U3 and then U1
            elseif decomposition==3
                inner_objective=0.5 * mu(1)*log2((1+w(1)-(g(1)^2 *(1+w(3)-g(2)^2+w(2)*g(2)^2)/...
                    ((1+w(2))*(1+w(3))-g(2)^2)))/(w(1)))...
                    +0.5 * mu(2)*log2( (1+w(2))/ w(2))...
                    +0.5 * mu(3)*log2(((1+w(3))*(1+w(2))-g(2)^2)/(w(3)*(1+w(2))));
            end
        %U3 is null
        else
            % This is sum-rate so the order of decoding does not matter
            if decomposition ==0
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(2))-g(1)^2)/(w(1)*w(2)));
                
            %Since U3 is null, this case reduces to decoding U2 first and then U1
            elseif decomposition==1
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(2))-g(1)^2)/(w(1)*(1+w(2))))...
                    +mu(2)*0.5*log2((1+w(2))/w(2));
                
            %Since U3 is null, this case reduces to decoding U1 first and then U2
            elseif decomposition==2
                inner_objective=mu(1)*0.5*log2((1+w(1))/w(1))+0.5*mu(2)*log2(...
               ((1+w(1))*(1+w(2))-g(1)^2)/((1+w(1))*w(2))); 
           
            %Since U3 is null, this case reduces to decoding U2 first and then U1
            elseif decomposition==3
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(2))-g(1)^2)/...
                    (w(1)*(1+w(2))))...
                    +mu(2)*0.5*log2((1+w(2))/w(2));
            end
        end
    else
        %U2 is null 
        if w(3)<Inf
            % This is sum-rate so the order of decoding does not matter
            if decomposition ==0
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(3))-g(3)^2)/(w(1)*w(3)));
                
            %Since U2 is null, this case reduces to decoding U3 first and
            %then U1
            elseif decomposition==1
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(3))-g(3)^2)/(w(1)*(1+w(3))))...
                    +mu(3)*0.5*log2((1+w(3))/w(3));
                
            %Since U2 is null, this case reduces to decoding U3 first and
            %then U1
            elseif decomposition==2
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(3))-g(3)^2)/(w(1)*(1+w(3))))...
                   +mu(3)*0.5*log2((1+w(3))/w(3));
            %Since U2 is null, this case reduces to decoding U3 first and
            %then U1   
            elseif decomposition==3
                inner_objective=mu(1)*0.5*log2(((1+w(1))*(1+w(3))-g(3)^2)/(w(1)*(1+w(3))))...
                   +mu(3)*0.5*log2((1+w(3))/w(3));
            end
            
        %U2 and U3 are null. Since we have only U1 as non-trivial, then
        %decomposition does not matter
        else
            inner_objective=mu(1)*0.5*log2((1+w(1))/w(1));
            
        end
    end
    
    
else
    
    if w(2)<Inf
        %U1 is null
        if w(3)<Inf
            
            %Decoding order does not matter
            if decomposition==0
                inner_objective=mu(2)*0.5*log2(((1+w(2))*(1+w(3))-g(2)^2)/(w(2)*w(3)));
            %Since U1 is null, all we have to do is to decode U3 first then
            %decode U2
            elseif decomposition==1
                inner_objective=mu(2)*0.5*log2(((1+w(2))*(1+w(3))-g(2)^2)/(w(2)*(1+w(3))))...
                    +mu(3)*0.5*log2((1+w(3))/w(3));
            %Since U1 is null, all we have to do is to decode U3 first then
            %decode U2  
            elseif decomposition==2
                inner_objective=mu(2)*0.5*log2(((1+w(2))*(1+w(3))-g(2)^2)/(w(2)*(1+w(3))))...
                    +mu(3)*0.5*log2((1+w(3))/w(3));
            %Since U1 is null, all we have to do is to decode U2 first then
            %decode U3  
            elseif decomposition==3
                inner_objective=mu(3)*0.5*log2(((1+w(3))*(1+w(2))-g(2)^2)/(w(3)*(1+w(2))))...
                    +mu(2)*0.5*log2((1+w(2))/w(2));
                
            end
        
        else
            %U1 and U2 are null
            if w(3)<Inf
                inner_objective=mu(3)*0.5*log2((1+w(3))/w(3));
        
            %U1, U2 and U3 are null
            else
                inner_objective=0;
            end
        
        end
    end
end

end