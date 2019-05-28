%% Function to implement alternating minization (as in alt_min.m)
% r: Read counts
% m: No of intervals
% n: Number of subclones 
% K: Max number of CNAs
% neta: Average coverage
% num_iter: Max number of iterations 

function [mu_est,L_est] = alt_min_fn_upd(r,m,n,K,neta,num_iter)

% mu_conv = false;      % True is estimates for mu converge
% L_conv = false;       % True if estimates for L converge
% mu_diff = 1e-3;
% L_diff = m*1e-3;

%% Initializations
mu_curr = sort(rand(n,1));
mu_curr = mu_curr/sum(mu_curr);
L_curr = 2*ones(m,n);
iter = 0;
obj_fn = [];

while (iter<num_iter) % && ~(obj_conv)
%     L_prev = L_curr;
    obj_est = sumsqr(r-neta*L_curr*mu_curr); % Initalize objective function with current value
    
    %% Searching for best L given mu
    r_L = r/neta- L_curr(:,1)*mu_curr(1);      % Removing contribution of normal genome
    mu_tum = mu_curr(2:end);             % Fractions corresponding to tumor subclones

    L_b = round(r_L*mu_tum'/(mu_tum'*mu_tum));  % Babai estimate
    L_b = min(L_b,K);                           % Capping max value of CNA in Babai estimate to K
    d = abs(r_L - L_b*mu_tum);                  % Bound for search

    [clone_mu,clone_ord] = sort(mu_tum,'descend');    % Sort clones in order of decreasing fractions
    L_best = zeros(m,n-1);                              % Best estimate of L

    for i = 1:m         % Same process for each region
        d_i = d(i);   % Error bound
        r_i = r_L(i)*ones(1,n-1);  % For purposes of finding bounds at each level of tree
        level = 1;  % Level of tree where search is (1 is top level, max depth is n-1)
        incr = false;                  % True if must update higher level l
        L_i = L_b(i,:);
        L_tmp = zeros(1,n-1);

        while true
            l_max = min(K,floor((r_i(level) + d_i)/clone_mu(level)));
            l_min = max(0,ceil((r_i(level) - d_i)/clone_mu(level))-(n-level)*K);

            if l_min>l_max      % No possible values
                if level~=1
                    level = level - 1;     % Go up one level and increment
                    incr = true;           % Already true
                else                       % Terminate while loop if at top level (exhausted search space)
                    break;
                end

            elseif ~incr
                L_tmp(clone_ord(level)) = l_min;

                if level == n-1   % Reached a leaf node
                    d_red = abs(r_i(level) - clone_mu(level)*l_min);
                    if d_red < d_i  % If solution is better, use it
                        d_i = d_red;    % Improved bound for search
                        L_i = L_tmp;
                    end
                    incr = true;
                else                            % if not at a leaf node
                    r_i(level+1) = r_i(level) - clone_mu(level)*l_min;
                    level = level+1;    % Descend to next level
                end

            elseif L_tmp(clone_ord(level))+1 > l_max

                if level~=1
                    level = level - 1;     % Go up one level and increment
                    incr = true;           % Already true
                else                       % Terminate while loop if at top level (exhausted search space)
                    break;
                end

            else                % Increment value at current level (incr = true and within bounds)
                L_tmp(clone_ord(level)) = L_tmp(clone_ord(level))+1;

                if level == n-1   % Reached a leaf node
                    d_red = abs(r_i(level) - clone_mu(level)*L_tmp(clone_ord(level)));

                    if d_red < d_i   % If solution is better, use it
                        d_i = d_red;    % Improved bound for search
                        L_i = L_tmp;
                    end
                    incr = true;                % Already true
                else                            % if not at a leaf node
                    r_i(level+1) = r_i(level) - clone_mu(level)*L_tmp(clone_ord(level));
                    level = level+1;    % Descend to next level
                    incr = false;
                end

            end

        end

        L_best(i,:) = L_i;
    end

    L_curr = [2*ones(m,1) L_best];
    %     L_conv = (sumsqr(L_curr-L_prev)<L_diff);

    %% Finding mu given L
%     mu_prev = mu_curr;
    mu_curr = proj_simplex(r/neta,L_curr);

    err = (r - neta*L_curr*mu_curr).^2;  % Absolute error

    %% Check objective function 

    iter = iter+1;
%         fprintf('Error in tumor purity at iter %d: %.2f\n', iter, norm(mu-mu_curr));
    obj_fn = [obj_fn, sumsqr(r-neta*L_curr*mu_curr)];
%     mu_arr = [mu_arr mu_curr];
%     cov_arr = [cov_arr var(err)];
%     obj_conv = (obj_fn(end)<obj_thresh);

    % Updating estimates
    if obj_fn(end)<obj_est
        obj_est = obj_fn(end);
        mu_est = mu_curr; % Best estimate of tumor fractions
        L_est = L_curr;   % Best estimates of CNAs
    end

    %% Random mu
    mu_curr = rand(n,1);
    mu_curr = sort(mu_curr/sum(mu_curr));    % Tumor fractions are in ascending order
    %     mu_conv = (norm(mu_curr-mu_prev)<mu_diff);

end    
    
end