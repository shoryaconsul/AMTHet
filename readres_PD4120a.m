%% Reading and running AltMin on PD4120a

num_iter = 1e5;     % Max no of iterations for Alt_Min
m = 335;     % No of intervals
n = 2;      % No of subclones
K = 3;      % Maximum number of copy number aberrations (CNAs)


dir_ip = 'ThetA_06_bootstrap\python\interval_files\';        % Relative path to input folder
dir_op = 'ThetA_06_bootstrap\python\output\';                % Relative path to output folder
pre = 'PD4120a';  %input prefix

%% Read THetA results
fname = strcat(dir_op,pre,'.n',int2str(n),'.results');
fid = fopen(fname,'r');
res_t = textscan(fid,'%s',4,'HeaderLines',1);
mu_t = str2num(res_t{1,1}{2,1});
L_split = strsplit(res_t{1,1}{3,1},':');
L_t = zeros(m,n-1);
for i=1:m
    L_t(i,:) = str2num(L_split{i});
end
fclose(fid);   
L_t = [2*ones(m,1) L_t];

%% Running AltMin
% Determining coverage per interval
fname = strcat(dir_ip,pre,'.intervals');
fid = fopen(fname,'r');
r_tn = textscan(fid,'%s %d %d %d %d %d','HeaderLines',1);
r_n = double(r_tn{1,6});
r_t = double(r_tn{1,5});
r_tscale = r_n(1)/r_t(1);
r = 2*r_tscale*r_t./r_n;
fclose(fid);

neta = 1;           % Coverage has been normalized
[mu_am,L_am] = alt_min_fn_upd(r,m,n,K,neta,num_iter); 

%% Plots
set(groot,'DefaultAxesFontSize',14);

if n==2
    figure(1);
    bar(L_t(:,2));
    xlabel('Interval number')
    ylabel('Copy number')
    set(gca,'FontName','Times New Roman','Color','W')
    set(gcf, 'Color','W')
    
    figure(2);
    bar(L_am(:,2));
    xlabel('Interval number')
    ylabel('Copy number')
    set(gca,'FontName','Times New Roman','Color','W')
    set(gcf, 'Color','W')
end