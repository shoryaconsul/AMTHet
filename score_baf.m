%% Function to implement alternating minization (as in alt_min.m)
% snv_baf: BAFs of SNPs (row vector)
% Li: Copy number of clones (row vector)
% mu: Tumor fraction (column vector)
% bnd: Range of values around snp_baf that are acceptable
% score: True if all computed bafs are consistent with snv_baf
function [score] = score_baf(snv_baf,Li,mu,bnd)
    n = length(Li); % Number of clones
    Lmax = max(Li); % Max copy number
    
    % Determining number of possible allelic configurations
    num_sol = 1;
    if rem(Lmax,2) == 1 % If Lmax is odd
        for c = 1: floor(Lmax/2)
            num_sol = num_sol*(2*c)^(sum(Li==2*c-1 | Li==2*c));
        end
        num_sol = num_sol*(Lmax+1)^(sum(Li==Lmax));
    else % If Lmax is even
        for c = 1: Lmax/2
            num_sol = num_sol*(2*c)^(sum(Li==2*c-1 | Li==2*c));
        end
    end

    comp_baf = zeros(num_sol,1);  % Compute BAFs for each configuration  
    
    Lcurr = zeros(1,n); % Allelic configuration (num of ref allele in each clone)
    for sol = 1:num_sol
       frac = [(Lcurr*mu) (Li-Lcurr)*mu]; % Temp variable for computation of BAF
       comp_baf(sol) = frac(1)/sum(frac);   % BAF computaton
       Lcurr = nextL(Lcurr,Li); % Next possible solution
    end
    
    score_mat = (snv_baf>comp_baf-bnd)&(snv_baf<comp_baf+bnd); % True if snp_baf values lie withing range around comp_baf values
    score_config = sum(score_mat,2)/length(snv_baf); % Fraction of values that fall in acceptable range
    score = max(score_config);    
    
end

%% Function to generate next possible Li based on current Li and inferred copy number
% Lcurr: Current Li: Single vector
% Li: Vector of copy numbers 
function Lnew = nextL(Lcurr,Lmax)
    if size(Lcurr) ~= size(Lmax)
        error('Incompatible input sizes');
    elseif size(Lcurr,1)~= 1
        error('Input must be row vector')
    else
        Lnew = Lcurr;   % result
        
        c = 1; % Set carry to 1
        i = 1;  % Start with first index
        while c>0     % While values have to be increased
            if i > length(Lcurr) 
                break;    % Exhausting possible solutions
            end
            
            tmp = Lcurr(i)+c;
            if tmp == Lmax(i)/2   % For even copy numbers, we disallow equal copies of both alleles
                tmp = tmp+1;
            end
            
            c = mod(tmp,Lmax(i)+1); % Carry
            s = rem(tmp,Lmax(i)+1); % Modulo sum

            Lnew(i) = s;
            i = i + 1;  % Increment index
                
        end
    end
    
end