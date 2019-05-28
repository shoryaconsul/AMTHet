%% Generating datasets with CNVs similar to that in SCNVSim
% Interval lengths drawn from exponential distribution

%% Parameters
m = 12;     % No of intervals
n = 4;      % No of subclones
K = 4;      % Maximum number of copy number aberrations (CNAs)
neta = 30;  % Average coverage
mu_min = 1/(2*n-1);
if n == 2   % Minimum fraction for single tumor subclone
    mu_min = 0.05;
end

Lr = 100;   % Length of reads
Lmax = 15e7;    % Length of chromosome

% params from SCNVsim
svrate = 0.2; %0.2;
hethomratio = 2;
lam = 1;    % Lambda for Exp rnd for intervals

num_rep = 10;       % Experiment repeated num_rep times to evaluate performance
dir = 'ThetA_06_bootstrap\python\interval_files\scnv\';  % output directory

for rep=0:num_rep-1
    L = 2*ones(m,n);
    
    Lseg = [1000];  % Breakpoints of intervals
    Lcurr = Lseg(1);      % Current position along genome
    
    while true % While end of chromosome not reached
        intlen = round(1e6*exprnd(lam));
        if Lcurr+intlen<Lmax
            Lseg = [Lseg Lcurr+intlen];
            Lcurr = Lcurr+intlen;
        else
            Lseg = [Lseg Lmax];
            break;
        end
    end       
    Lind = sort(randperm(length(Lseg)-1,m));    % Select m intervals
   
    numcl = 1;  % Number of sublcones generated thus far
    cl = 2*ones(m,1); % New subclone initialized
    
    while numcl<n
        for i=1:m % Each interval
            if rand<svrate
                if randi(4)==1  % introduce tandem duplication
                    cl(i) = min(round(exprnd(3))+1,K);
                elseif randi(4) == 2 % deletion with copy number loss
                    if rand>hethomratio/(hethomratio+1)
                        cl(i)=1;
                    end
                end
            end
        end

        if ~ismember(cl',L(:,1:numcl)','rows') % Comparing previous subclones to new one
            L(:,numcl+1) = cl;
            numcl = numcl+1;    % New subclone has been generated
            cl = 2*ones(m,1);   % For next subclone generation
        end
    end
    
    
    % Tumor fractions are each at least 20%
    mu = rand(n,1); mu = mu_min + (1-n*mu_min)*mu/sum(mu); mu = sort(mu);
    nr = zeros(m,1);
    nr_ref = zeros(m,1);
    r_wn = neta*L*mu;   % Not noisy

    for i=1:m
        Li = Lseg(Lind(i)+1) - Lseg(Lind(i));
        nr_wn = r_wn(i)*Li/Lr+neta;
        phi = 0.02; % Noise level is r_wn*phi
        noise = randn;
        r = r_wn(i) + phi*r_wn.*noise;
        nr(i) = round(nr_wn + phi*nr_wn.*noise);
        nr_ref(i) = round(neta*Li/Lr*2+neta); % Reads from normal genome (used for weights in ThetA)
%     fprintf('Ref count at rep %d is %d\n',rep,nr_ref)
    end
    
    fname = strcat(dir,'interval_count_n',int2str(n),int2str(rep));
    fid = fopen(fname,'w+');
    fprintf(fid,'#ID chrm start end tumorCount normalCount\n');
    for i=1:m
        fprintf(fid,'%d %d %d %d %d %d\n',i,3,Lseg(Lind(i)),Lseg(Lind(i)+1),nr(i),nr_ref(i));
    end
    
    fname = strcat(dir,'muC_n',int2str(n),int2str(rep));
    fid = fopen(fname,'w+');
    fprintf(fid,'mu C\n');
    str = strip(repmat('%.3f,',[1,n]),'right',',');
    fprintf(fid,str,mu);
    fprintf(fid,' ');
    str = strcat(strip(repmat('%d,',[1,n-1]),'right',','),':');
    fprintf(fid,str,L(:,2:end)');

    fclose(fid);
end
               
            