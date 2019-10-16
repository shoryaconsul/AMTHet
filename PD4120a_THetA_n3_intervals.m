%% Script to determine intervals to be selected for n=3 run of THetA

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

fname = strcat(dir_ip,pre,'.intervals');
fid = fopen(fname,'r');
r_tn = textscan(fid,'%s %d %d %d %d %d','HeaderLines',1);
intlen = r_tn{1,4}-r_tn{1,3}; % Lengts of intervals
fclose(fid);

int_sel = (L_t(:,2)~=2).*(intlen>2e6);

fname = strcat(dir_ip,pre,'.n3.intervals');
fid = fopen(fname,'w');
fprintf(fid,'#ID chrm start end tumorCount normalCount\n');
for i=1:m
    if int_sel(i)
        %s = strcat(r_tn{1,1}{i},' ',int2str(r_tn{1,2}{i}),' ',int2str(r_tn{1,3}{i}),
        s = sprintf('%s %d %d %d %d %d\n',r_tn{1,1}{i},r_tn{1,2}(i),r_tn{1,3}(i), ...
            r_tn{1,4}(i),r_tn{1,5}(i),r_tn{1,6}(i));
        fprintf(fid,s);
    end
end
fclose(fid);
