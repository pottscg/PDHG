%% SKETCH_PROJECTION
% Build a projection matrix based on the requested type. Note that 
% FJLT only works if m is a power of 2. 
% INPUT:
%   m = original dimension of data
%   d = dimension data is being projected onto
%   type = string indicating which type of projection is desired
%
function PSI = Sketch_Projection(m,d,type)
    
    if(strcmp(type,'Gaussian'))
        %create d x m projection matrix 
        PSI = reshape(rand(d*m,1),[d,m]);

        %normalize over rows
        for j = 1:d
            PSI(j,:) = PSI(j,:)./sum(PSI(j,:));
        end
    else if(strcmp(type, 'Ternary'))
            %create d x m projection matrix 
            seed = reshape(rand(d*m,1),[d,m]);
            PSI = zeros(size(seed));
            PSI(seed<=0.25) = -1;
            PSI(seed>=0.75) = 1;
            
        else if(strcmp(type,'FJLT'))
                %sparsity coeff
                q = 0.25;
                seed = rand(d,m);
                P = zeros(size(seed));
                P(seed > (1-q)) = normrnd(0,q^(-1),1,sum(sum(seed > (1-q))));
                H = (1/sqrt(m))*hadamard(m);
                seed = rand(1,m);
                seed(seed < 0.5) = -1;
                seed(seed >=0.5) = 1;
                D = diag(seed);
                PSI = P*H*D;
            end
        end
    end

end
