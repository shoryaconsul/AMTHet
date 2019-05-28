%% Script to make plots for ACM-BCB paper from results on spreadsheet

clear
close all

set(groot,'DefaultAxesFontSize',14);

%% For m=10 intervals, varying number of populations (n=2,3,4,5)
n_arr = 2:5;
pu_err = [0.057,0.466,0.037;
          0.527,0.141,0.089;
          NaN,0.122,0.077;
          NaN,0.067,0.074];
tf_err = [0.08,0.701,0.053;
          0.672,0.289,0.143;
          NaN,0.255,0.165;
          NaN,0.258,0.23];

cna_err = [0.032,0.248,0.0;
          0.248,0.15,0.135;
          NaN,0.11,0.159;
          NaN,0.094,0.16];
      
cpu_time = [1, 12, 7.1;
            7200, 32, 8.4;
            NaN, 41, 16.1;
            NaN, 29, 24];
      
figure(1);
h = plot(n_arr,pu_err);
set(h,{'Marker'},{'o';'s';'^'},{'Color'},{'r';'b';'k'},{'MarkerEdgeColor'},{'r';'b';'k'},...
    {'MarkerFaceColor'},{'r';'b';'k'},{'LineStyle'},{'None';'None';'None'})
xlabel('K')
ylabel('Purity error')
xlim([1.5,5.5])
ylim([0,0.7])
legend({'THetA','MixClone','AMTHet'})
set(gca,'FontName','Times New Roman','Color','W')
set(gcf, 'Color','W')

figure(2);
h = plot(n_arr,tf_err);
set(h,{'Marker'},{'o';'s';'^'},{'Color'},{'r';'b';'k'},{'MarkerEdgeColor'},{'r';'b';'k'},...
    {'MarkerFaceColor'},{'r';'b';'k'},{'LineStyle'},{'None';'None';'None'})
xlabel('K')
ylabel('Tumor fraction error')
xlim([1.5,5.5])
ylim([0,0.7])
legend({'THetA','MixClone','AMTHet'})
set(gca,'FontName','Times New Roman','Color','W')
set(gcf, 'Color','W')

figure(3);
h = plot(n_arr,cna_err);
set(h,{'Marker'},{'o';'s';'^'},{'Color'},{'r';'b';'k'},{'MarkerEdgeColor'},{'r';'b';'k'},...
    {'MarkerFaceColor'},{'r';'b';'k'},{'LineStyle'},{'None';'None';'None'})
xlabel('K')
ylabel('CNA error')
xlim([1.5,5.5])
ylim([0,0.3])
legend({'THetA','MixClone','AMTHet'}')
set(gca,'FontName','Times New Roman','Color','W')
set(gcf, 'Color','W')

figure(4);
h = semilogy(n_arr,cpu_time);
set(h,{'Marker'},{'o';'s';'^'},{'Color'},{'r';'b';'k'},{'MarkerEdgeColor'},{'r';'b';'k'},...
    {'MarkerFaceColor'},{'r';'b';'k'},{'LineStyle'},{'None';'None';'None'})
xlabel('K')
ylabel('Time')
xlim([1.5,5.5])
% ylim([0,0.3])
legend({'THetA','MixClone','AMTHet'}')
set(gca,'FontName','Times New Roman','Color','W')
set(gcf, 'Color','W')
