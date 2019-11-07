%% Experiment to check if BAF (SNP) information improves AMTHet

%% Generating datasets with CNVs similar to that in SCNVSim
% Interval lengths drawn from exponential distribution

%% Parameters
num_iter = 1e5; % Max no of iterations for AMTHet
m = 10;     % No of intervals
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

%% Generating CNA matrix
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

%% Generating read counts
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
    

%% Generating heterozygous SNP positions
snvrate = 1e3*3; % Mean of exp RV to draw heterozygous SNP positions

stpos = Lseg(Lind); % Start positions of each selected segment
endpos = Lseg(Lind+1); % End positions of each selected segment
snv_seg = [];   % Segment for each SNP
snv_pos = [];   % SNP positions

for i = 1:length(stpos)
    currpos = stpos(i);
    while currpos < endpos(i)
        incr = exprnd(snvrate); % Distance to next SNP
        if currpos+incr > endpos
            break;
        else
            currpos = currpos + incr;
            snv_seg = [snv_seg i]; % Segment number
            snv_pos = [snv_pos currpos]; % Adding SNP position
        end
    end
end

% Generating allelic configurations for each segment
seg_config = zeros(m,n,2); % Third dimension is copy number of other allele
for i = 1:m
    for j = 1:n
        Lnum = L(i,j); % Copy number
        if Lnum == 1
            ref = randn<0;  % true if reference alele present
            seg_config(i,j,:) = int8([ref ~ref]);
        elseif Lnum == 2
            seg_config(i,j,:) = int8([1 1]);
        else
            while true
                ref = randi(Lnum); % number of reference alleles
                if ref ~=Lnum/2 || rem(Lnum,2)==1 % Avoid AABB or the like
                    seg_config(i,j,:) = int8([ref Lnum-ref]);
                    break;
                end
            end
        end
    end
end

% Generating BAFs for heterozygous SNP positions
BAF_N_MIN = 0.4;
BAF_N_MAX = 0.6;
baf_n = BAF_N_MIN + (BAF_N_MAX-BAF_N_MIN)*rand(1,length(snv_pos));
baf_mat = [baf_n', 1-baf_n'];

snv_baf = zeros(1,length(snv_pos)); % BAFs of reference allele for each SNV
for i = 1:m
    ind = find(snv_seg == i);
    ref_frac = squeeze(seg_config(i,:,:))'*mu; % Fraction of ref and non-ref allele
    baf_cl = baf_mat(ind,:).*repmat(ref_frac',length(ind),1); % BAF for each clone
    snv_baf(ind) = baf_cl(:,1)./sum(baf_cl,2);    
end
    

%% Running AMTHet and checking if solution is consistent with BAF info
baf_t = 0.6;    % Fraction of SNP BAF
baf_bnd = 0.1;  % Range of acceptable BAF error 
num_rep = 5;       % Experiment repeated num_rep times to evaluate performance

mu_err = zeros(2,num_rep); % L2 norm of diff in tumor fractions
L_err = zeros(2,num_rep);   % L2 norm of difference in CNA matrix
comp_time = zeros(2,num_rep);   % Time taken for each run

for rep = 0:num_rep-1
    tic
    [mu_snp,L_snp] = alt_min_fn_snp(r,snv_baf,snv_seg,m,n,K,neta,baf_t,baf_bnd,num_iter);
    comp_time(1,rep+1) = toc;

    [ord_snp, L_snps] = L_perm(L(:,2:end),L_snp(:,2:end));
    L_snp = [2*ones(m,1) L_snps(:,ord_snp)];
    mu_snp = [mu_snp(1); mu_snp(ord_snp+1)];

    mu_err(1,rep+1) = norm(mu-mu_snp');
    L_err(1,rep+1) = norm(L(:)-L_snp(:))/(m*(n-1)); 

    tic
    [mu_am,L_am] = alt_min_fn_upd(r,m,n,K,neta,num_iter);
    comp_time(2,rep+1) = toc;
    
    [ord_am, L_ams] = L_perm(L(:,2:end),L_am(:,2:end));
    L_am = [2*ones(m,1) L_ams(:,ord_am)];
    mu_am = [mu_am(1); mu_am(ord_am+1)];

    mu_err(2,rep+1) = norm(mu-mu_am');
    L_err(2,rep+1) = norm(L(:)-L_am(:))/(m*(n-1)); 
end