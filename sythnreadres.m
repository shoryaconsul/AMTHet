%% Reading THetA and MixClone results and comparing vs AMTHet results for synthetic data
num_rep = 10;       % Experiment repeated num_rep times to evaluate performance
num_iter = 1e5;     % Max no of iterations for Alt_Min
m = 12;     % No of intervals
n = 4;      % No of subclones
K = 4;      % Maximum number of copy number aberrations (CNAs)
theta_read = false; % True if THetA results are to be read

% dir_ip = 'ThetA_06_bootstrap\python\interval_files\amp\';        % Relative path to input folder
% dir_op = 'ThetA_06_bootstrap\python\output\amp\';                % Relative path to output folder
% dir_ip = 'ThetA_06_bootstrap\python\interval_files\rand\';        % Relative path to input folder
% dir_op = 'ThetA_06_bootstrap\python\output\rand\';                % Relative path to output folder
dir_ip = 'ThetA_06_bootstrap\python\interval_files\scnv\';        % Relative path to input folder
dir_op = 'ThetA_06_bootstrap\python\output\scnv\';                % Relative path to output folder
pre = strcat('interval_count_n',int2str(n));  %input prefix

mu_err = zeros(3,num_rep);  % Purity error
mu_diff = zeros(3,num_rep); % L2 norm of diff in tumor fractions
L_err = zeros(3,num_rep);   % L2 norm of difference in CNA matrix
L_diff = zeros(3,num_rep);  % L0 norm of difference in CNA matrix

tic
h = waitbar(0, 'Reading results');

for rep = 0:num_rep-1
    
    % Read in ground truth
    fname = strcat(dir_ip,'muC_n',int2str(n),int2str(rep));
    fid = fopen(fname,'r');
    muC = textscan(fid,'%s',2,'HeaderLines',1);
    mu_str = muC{1,1}{1,1};
    L_str = muC{1,1}{2,1};
    
    mu = str2num(mu_str);   % Used str2num instead of str2double as it allows ',' as a delimiter
    L_split = strsplit(L_str,':');
    L = zeros(m,n-1);
    for i=1:m
        L(i,:) = str2num(L_split{i});
    end
    fclose(fid);
    L = [2*ones(m,1) L];
    
    % Read THetA results
    if theta_read
        fname = strcat(dir_op,pre,int2str(rep),'.n',int2str(n),'.results');
        fid = fopen(fname,'r');
        res_t = textscan(fid,'%s',4,'HeaderLines',1);
        mu_t = str2num(res_t{1,1}{2,1});
        L_split = strsplit(res_t{1,1}{3,1},':');
        L_t = zeros(m,n-1);
        for i=1:m
            L_t(i,:) = str2num(L_split{i});
        end
        fclose(fid);   
        [ord_t, L_ts] = L_perm(L(:,2:end),L_t);
        L_t = [2*ones(m,1) L_ts(:,ord_t)];
        mu_t = [mu_t(1) mu_t(ord_t+1)];

        mu_diff(1,rep+1) = abs(mu(1)-mu_t(1)); % Purity error
        L_err(1,rep+1) = norm(L(:)-L_t(:))/(m*(n-1)); % Tumor fraction error
        mu_err(1,rep+1) = norm(mu-mu_t); % CNV error
        L_diff(1,rep+1) = nnz(L-L_t)/(m*(n-1)); % Fraction of diff CNV entries
    end
    %% Read MixClone results
    L_mix_raw = 2*ones(m,n);
    fname = strcat(dir_op,'n',int2str(n),int2str(rep),'.MixClone.segments');
    fopen(fname,'r');
    for i=1:m
        res_m = textscan(fid,'%s %d %d %d %d %d %f %s %f %d %s %f %d\n','HeaderLines',1);
        if ~isempty(res_m{end}) % if subclone cluster detected
            L_mix_raw(i,res_m{end}) = res_m{10};    
        end
    end
    fclose(fid);   
    
    [ord_mix,L_mix] = L_perm(L,L_mix_raw);
    
    mu_mix = zeros(1,n);
    fname = strcat(dir_op,'n',int2str(n),int2str(rep),'.MixClone.summary');
    fopen(fname,'r');
    res_format = "%s %s %s %s "+strjoin(repmat("%f",1,n)," ")+"\n";
    res_m = textscan(fid,res_format,'HeaderLines',4);
    for i=1:n
        mu_mix(i) = res_m{4+i};
    end
    mu_mix = mu_mix(ord_mix);
    fclose(fid);
    
    mu_diff(2,rep+1) = abs(mu(1)-mu_mix(1)); % Purity error 
    L_err(2,rep+1) = norm(L(:)-L_mix(:))/(m*(n-1)); % Tumor fraction error
    mu_err(2,rep+1) = norm(mu-mu_mix); % CNV error
    L_diff(2,rep+1) = nnz(L-L_mix)/(m*(n-1)); % Fraction of diff CNV entries
    
    %% Run Alt_Min 
    % Determining coverage per interval
    fname = strcat(dir_ip,'interval_count_n',int2str(n),int2str(rep));
    r_tn = dlmread(fname,' ',1,4);
    r = 2*r_tn(:,1)./r_tn(:,2); 
    
    neta = 1;           % Coverage has been normalized
%     tic
    [mu_am,L_am] = alt_min_fn_upd(r,m,n,K,neta,num_iter); 
%     L_am = L_am(:,2:end);
%     toc
    
    [ord_am, L_ams] = L_perm(L(:,2:end),L_am(:,2:end));
    L_am = [2*ones(m,1) L_ams(:,ord_am)];
    mu_am = [mu_am(1); mu_am(ord_am+1)];
    
    mu_diff(3,rep+1) = abs(mu(1)-mu_am(1));
    mu_err(3,rep+1) = norm(mu-mu_am');
    L_err(3,rep+1) = norm(L(:)-L_am(:))/(m*(n-1)); 
    L_diff(3,rep+1) = nnz(L-L_am)/(m*(n-1));
    
        
    waitbar(rep/num_rep);
end
    
close(h);
toc
    
fprintf('Median purity error (Normal,Tumor): %.3f, %.3f %.3f\n',median(mu_diff,2))
fprintf('Mean purity error (Normal,Tumor): %.3f, %.3f %.3f\n',mean(mu_diff,2))
fprintf('Median tf error (Normal,Tumor): %.3f, %.3f %.3f\n',median(mu_err,2))
fprintf('Mean tf error (Normal,Tumor): %.3f, %.3f %.3f\n',mean(mu_err,2))
fprintf('Median CNA error (Normal,Tumor): %.3f, %.3f %.3f\n',median(L_err,2))
fprintf('Mean CNA error (Normal,Tumor): %.3f, %.3f %.3f\n',mean(L_err,2))
fprintf('Median CNA diff (Normal,Tumor): %.3f, %.3f %.3f\n',median(L_diff,2))
fprintf('Mean CNA diff (Normal,Tumor): %.3f, %.3f %.3f\n\n',mean(L_diff,2))