close all;
clear all;
clc;

Fs = 250; % Sampling frequency
total_dur = 580; % Duration of entire main task in seconds
total_dur_short = 340;
stim_dur = 40; % Duration of stimulus trials in seconds
rest_dur = 20; % Duration of rest trials in seconds
N = Fs*total_dur; % Total number of samples
N_short = Fs*total_dur_short;
trial_num = 10; % Number of stimulus trials

main_data_path = "F:\Documents\SUT M.Sc\Research\Matlab codes\newm data_auditory\";
out_path = "F:\Documents\SUT M.Sc\Research\Matlab codes\newm data_auditory\profiles\";
res_path = "F:\Documents\SUT M.Sc\Research\Matlab codes\newm data_auditory\figs\";


patients = ["nazari", "omidvar", "madannejad",...
            "ostadkarim", "rahchamani", "yousefi",...
            "ghorayshi", "pakdel", "saeedi",...
            "daneshpajoh", "moshtagh", "nokhostin",...
            "abdi", "jabari", "farayei",...
            "rezaei", "shishechi", "mahmoudnasab",...
            "ghasemi", "khalili", "mostafavi",...
            "hosseinnejad", "shafaghi", "souri",...
            "naghdi", "rezvani", "khalili-moslem",...
            "alizade", "teymouri", "mohammadganji",...
            "fathi", "jafari", "karimi",...
            "nematpour", "yousefi-abalfazl"]; 
        
MMSE = [28 23 27 ...
        21 21 23 ...
        19 24 26 ...
        28 20 30 ...
        0 18 19 ...
        25 25 18 ...
        13 13 30 ...
        29 14 23 ...
        24 28 23 ...
        26 16 22 ...
        26 21 15 ...
        20 21];

Label = ["Normal" "Mild" "Normal" ...
         "Mild" "Mild" "MCI" ...
         "Mild" "Normal" "Normal" ...
         "Normal" "Mild" "Normal" ...
         "?" "Mild" "Mild" ...
         "Mild" "MCI" "MCI" ...
         "MCI" "Mild" "Mild" ...
         "Normal" "Mild" "Moderate" ...
         "Normal" "Normal" "MCI" ...
         "MCI" "Mild" "Mild" ...
         "Normal" "MCI" "Mild" ...
         "Mild" "Mild"];

     Pa_num = length(patients); % Total number of participants
Pa_code = "S"+[1:Pa_num]; % Participant code

Ignore = repelem("n",Pa_num);
Ignore(6) = "y"; Ignore(13) = "y";

% List of electrodes (10/20 system)
ch_list = ["Fp1","Fp2","F7","F3","Fz","F4","F8","T7","C3","Cz","C4","T8","P7","P3","Pz","P4","P8","O1","O2"];
ch_cell = {"FP1","FP2","F7","F3","FZ","F4","F8","T7","C3","CZ","C4","T8","P7","P3","PZ","P4","P8","O1","O2"};
ch_num = length(ch_list); % Number of electrodes

Fz_idx = find(ch_list == "Fz"); % Index of the Fz channel
Pz_idx = find(ch_list == "Pz"); % Index of the Pz channel

% 2D channel locations
chR = [0.511110000000000;0.511110000000000;0.511110000000000;0.333330000000000;0.255560000000000;0.333330000000000;0.511110000000000;0.511110000000000;0.255560000000000;0;0.255560000000000;0.511110000000000;0.511110000000000;0.333330000000000;0.255560000000000;0.333330000000000;0.511110000000000;0.511110000000000;0.511110000000000];
chPhi = [-18;18;-54;-39;0;39;54;-90;-90;90;90;90;-126;-141;180;141;126;-162;162]; %degree
% polarplot(deg2rad(chPhi), chR, 's')
% text(deg2rad(chPhi), chR, ch_list)



% List of sites
bch_num = 45; %14,28,4
ch1_bp = ["Fp1","F3","Fz","Fz","F4","C3","Cz","Cz","C4","P3","Pz","Pz","P4","O1","F3","F4","T7","C3","Cz","C4","T8","T7","C3","Cz","C4","T8","P3","P4"];
ch2_bp = ["Fp2","F7","F3","F4","F8","T7","C3","C4","T8","P7","P3","P4","P8","O2","Fp1","Fp2","F7","F3","Fz","F4","F8","P7","P3","Pz","P4","P8","O1","O2"];
bpSites = ch1_bp+"-"+ch2_bp;


theta_range = [4 8]; % in Hz
gamma_range = [39 41]; % in Hz
adj_40 = [38 42]; % in Hz

gr_inc = (Ignore=="n");

%% GES on channel data

W = 20; % Window length in seconds
nw = Fs*W; % Total number of samples in a window length of W sec.
f = Fs*(0:(nw/2))/nw;

freq_idx = f>=adj_40(1) & f<=adj_40(2); % Index of frequecy range aroud 40Hz
f40_idx = find(f == 40);

stim_idx = [1:2 4:5 7:8 10:11 13:14 16:17 19:20 22:23 25:26 28:29]; % Index of stimulus windows
% stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
rest_idx = [3 6 9 12 15 18 21 24 27]; % Index of rest windows
% rest_idx = [3 6 9 12 15]; % Index of rest windows

% stim_idx = [1:40 61:100 121:160 181:220 241:280 301:340]';
% rest_idx = [41:60 101:120 161:180 221:240 281:300]';

% stim_idx = [1:20; 31:50; 61:80; 91:110; 121:140; 151:170]';
% rest_idx = [21:30 51:60 81:90 111:120 141:150]';

% stim_idx = [1:20 31:50 61:80 91:110 121:140 151:170 181:200 211:230 241:260 271:290]';
% rest_idx = [21:30 51:60 81:90 111:120 141:150 171:180 201:210 231:240 261:270]';

% stim_idx = [1:4; 7:10; 13:16; 19:22; 25:28; 31:34];
% rest_idx = [5:6; 11:12; 17:18; 23:24; 29:30];

load(fullfile(main_data_path, sprintf("f_ch_w%d_Pa.mat", W)))


%% frequency spectrum
psd = fData_allWin_allPa;
psd(2:end-1,:,:,:) = psd(2:end-1,:,:,:)/2;
psd = abs(psd).^2;
psd(2:end-1,:,:,:) = psd(2:end-1,:,:,:)*2;

specs = mean(nanmean(psd(:, :, stim_idx(:), gr_inc), 3),2);
specr = mean(nanmean(psd(:, :, rest_idx(:), gr_inc), 3),2);


%% plots of psd

param = specs(1:5:end,:,:,:);
figure
colors = cbrewer('seq', 'Blues',9);
plot(f(1:5:end),pow2db(squeeze(param))', 'Color', [colors(8,:) 0.1])
hold on
param = specr(1:5:end,:,:,:);
plot(f(1:5:end),pow2db(squeeze(param))', 'Color', [colors(5,:) 0.1])

param = specs(1:5:end,:,:,:);
stdshade(pow2db(squeeze(param))', 0.3, colors(8,:), f(1:5:end), 2)
hold on
param = specr(1:5:end,:,:,:);
stdshade(pow2db(squeeze(param))', 0.3, colors(5,:), f(1:5:end), 2)

xlim([1 70])
ylim([-70 -10])
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')

% set(gcf,'Units','Inches');
% set(gcf,'PaperPosition',[0 0 3 2],'PaperUnits','Inches','PaperSize',[3, 2])
% fname = res_path+sprintf('PSD_stim');
% print(gcf, '-painters', '-dpdf', '-r300', fname+'.pdf');
% saveas(gcf, fname+'.fig')

%% GES

idx_40 = 41; % w=2 =>5, w=20 =>41
temp = (abs(psd(freq_idx, :, :, :)));
temp_z = zscore(temp,0,1);
gess = squeeze(nanmean(temp_z(idx_40,:,stim_idx(:),gr_inc),[ 3]));
gesr = squeeze(nanmean(temp_z(idx_40,:,rest_idx(:),gr_inc),[ 3]));

ges_diff = (gess - gesr);
ges_Norm = ges_diff./gesr;


%% topography plot
figure;
% subplot(1,2,1)
% plot_topography(ch_cell, mean(gess,2))
% 
% title('Stimulus','FontSize',6)
% colors = cbrewer('seq', 'Blues',20);
% colormap(colors);
% caxis([min([ mean(gess,2); mean(gesr,2)]) max([ mean(gess,2); mean(gesr,2)])])
% c = colorbar();c.FontSize=5;

% subplot(1,2,2)
plot_topography(ch_cell, mean(gesr,2))

title('Rest','FontSize',6)
colors = cbrewer('seq', 'Blues',20);
colormap(colors);
caxis([min([ mean(gess,2); mean(gesr,2)]) max([ mean(gess,2); mean(gesr,2)])])
c = colorbar();c.FontSize=5;

% set(gcf,'Units','Inches');
% set(gcf,'PaperPosition',[0 0 2 3],'PaperUnits','Inches','PaperSize',[2, 3])
% fname = res_path+sprintf('GES_topography_rest');
% print(gcf, '-painters', '-dpdf', '-r300', fname+'.pdf');
% saveas(gcf, fname+'.fig')

%% bar plot
star_loc = 7;


x = gess(Fz_idx,:);
y = gesr(Fz_idx,:);

colors = cbrewer('seq', 'Blues', 9);
% Plot error bars
figure
bar(1, mean(x), 'FaceColor',  colors(8, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on
bar(2, mean(y), 'FaceColor',  colors(5, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
% Plot samples
scatter(0.8:0.4/length(x):1.2-0.4/length(x), x, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(8, : ),...
        'LineWidth',1, 'SizeData',5)
scatter(1.8:0.4/length(y):2.2-0.4/length(y), y, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(5, : ),...
        'LineWidth',1, 'SizeData',5)
% Plot lines
fp = plot([0.8:0.4/length(x):1.2-0.4/length(x); 1.8:0.4/length(y):2.2-0.4/length(y)],...
     [x; y], 'Color', [0.5, 0.5, 0.5, 0.1]);



% Test of significance
% [~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
[~, pval, ci, stats] = ttest(x, y)
mysigstar(gca, [1 2], star_loc, pval);


set(h(1), 'marker', 'none'); 
set(gca, 'xtick', [1 1.5 2], 'xticklabel', {'Stimulus', 'Fz', 'Rest',},'xlim', [0.5 2.5],'XTickLabelRotation',45);
ylabel('Gamma entrainment score'); xlabel('Status');

% set(gcf,'Units','Inches');
% set(gcf,'PaperPosition',[0 0 2 3],'PaperUnits','Inches','PaperSize',[2, 3])
% fname = res_path+sprintf('GES_barplot');
% print(gcf, '-painters', '-dpdf', '-r300', fname+'.pdf');
% saveas(gcf, fname+'.fig')


%%
% plv

% W = 20; % Window length in seconds
% nw = Fs*W; % Total number of samples in a window length of W sec.
% f = Fs*(0:(nw/2))/nw;
% 
% freq_idx = f>=adj_40(1) & f<=adj_40(2); % Index of frequecy range aroud 40Hz
% f40_idx = find(f == 40);

% stim_idx = [1:2 4:5 7:8 10:11 13:14 16:17 19:20 22:23 25:26 28:29]; % Index of stimulus windows
% % stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
% rest_idx = [3 6 9 12 15 18 21 24 27]; % Index of rest windows
% % rest_idx = [3 6 9 12 15]; % Index of rest windows

% stim_idx = [1:40 61:100 121:160 181:220 241:280 301:340]';
% rest_idx = [41:60 101:120 161:180 221:240 281:300]';

% stim_idx = [1:20; 31:50; 61:80; 91:110; 121:140; 151:170]';
% rest_idx = [21:30 51:60 81:90 111:120 141:150]';

% stim_idx = [1:20 31:50 61:80 91:110 121:140 151:170 181:200 211:230 241:260 271:290]';
% rest_idx = [21:30 51:60 81:90 111:120 141:150 171:180 201:210 231:240 261:270]';

% stim_idx = [1:4; 7:10; 13:16; 19:22; 25:28; 31:34];
% rest_idx = [5:6; 11:12; 17:18; 23:24; 29:30];

% load(fullfile(main_data_path, sprintf("fBP45FPNorm_ch_w%d_Pa.mat", W)))
% load(fullfile(main_data_path, sprintf("tBP45FPNormAnalytic_fil39754025_ch_w%d_Pa.mat", W)))


%%
phi = angle(tData_allWin_allPa);
plv = ones(bch_num,bch_num,size(phi,3),Pa_num);

for ich=1:bch_num-1
    for jch=ich+1:bch_num
        
      deltaphi = exp( 1i.*(phi(:,ich,:,:) - phi(:,jch,:,:)) );
      value = abs(sum(deltaphi,1))./nw;
      plv(ich, jch, :, :) = value;
      plv(jch, ich, :, :) = value;
    end
end

%%
% save(fullfile(main_data_path, sprintf("PLV_FP45s1s2_w%d_Pa.mat", W)),'plv');
load(fullfile(main_data_path, sprintf("PLV_FP45s1s2_w%d_Pa.mat", W)))
%%

plvs = squeeze(nanmean(plv(:,:,stim_idx(:),gr_inc),3));
plvr = squeeze(nanmean(plv(:,:,rest_idx(:),gr_inc),3));


%%
sub = 1;

figure
heatmap(plvs(:, :, sub))
% heatmap(plvs(:, :, sub)./norm(plvs(:,:,sub),'fro'))
caxis([0 1])
figure
heatmap(plvr(:, :, sub))
% heatmap(plvr(:, :, sub)./norm(plvs(:,:,sub),'fro'))
caxis([0 1])


%% FP 45 all sites
% regenerating sites information
%long,local sites  
num = 0;
sites = [];
sitesMat = zeros(ch_num,ch_num);
loc = [1:2 4:6 14:16 18:19];
for ich=loc
    for jch=loc(find(loc>ich))
        num = num+1;
        sites(num, 1) = ich;
        sites(num, 2) = jch;
        sites(num, 3) = num;
        sitesMat(ich,jch) = 1;
        sitesMat(jch,ich) = 1;
    end
end

%% making mask for each group
mask_f = zeros(bch_num,bch_num);
mask_p = zeros(bch_num,bch_num);
mask_fp = zeros(bch_num,bch_num);

regionf = [1:2 4:6];
regionp = [14:16 18:19];
for s1=1:length(sites)
    for s2=1:length(sites)
        ch_sets = unique([sites(s1,1:2) sites(s2,1:2)]);
        if( length(ch_sets)==4 )
            if( ~isempty(find(regionf==sites(s1,1))) && ~isempty(find(regionf==sites(s1,2))) &&...
                    ~isempty(find(regionf==sites(s2,1))) && ~isempty(find(regionf==sites(s2,2))) )
                mask_f(s1,s2) = 1;
                
                
            elseif( ~isempty(find(regionp==sites(s1,1))) && ~isempty(find(regionp==sites(s1,2))) &&...
                    ~isempty(find(regionp==sites(s2,1))) && ~isempty(find(regionp==sites(s2,2))) )
                mask_p(s1,s2) = 1;
                
            elseif( ( ~isempty(find(regionf==sites(s1,1))) && ~isempty(find(regionp==sites(s1,2))) ||...
                    ~isempty(find(regionp==sites(s1,1))) && ~isempty(find(regionf==sites(s1,2))) ) &&...
                    (~isempty(find(regionf==sites(s2,1))) && ~isempty(find(regionp==sites(s2,2)))  ||...
                    ~isempty(find(regionp==sites(s2,1))) && ~isempty(find(regionf==sites(s2,2))) ) )
                mask_fp(s1,s2) = 1;
            end
        end
    end
end
   
NactiveF = sum(mask_f==1,'all');
NactiveFP = sum(mask_fp==1,'all');
NactiveP = sum(mask_p==1,'all');

mask_res = ones(bch_num,bch_num) - (mask_f + mask_fp + mask_p);
Nactiveres = sum(mask_res==1,'all');

figure
heatmap(mask_f)

figure
heatmap(mask_p)

figure
heatmap(mask_fp)

figure
heatmap(mask_res)

%% mask figures

[r,c] = find(mask_f==1);
adj_mat_f = zeros(ch_num,ch_num);
adj_mat_f(sites(r,1),sites(r,2)) = 1;
adj_mat_f(sites(c,1),sites(c,2)) = 1;

[r,c] = find(mask_fp==1);
adj_mat_fp = zeros(ch_num,ch_num);
adj_mat_fp(sites(r,1),sites(r,2)) = 1;
adj_mat_fp(sites(c,1),sites(c,2)) = 1;

[r,c] = find(mask_p==1);
adj_mat_p = zeros(ch_num,ch_num);
adj_mat_p(sites(r,1),sites(r,2)) = 1;
adj_mat_p(sites(c,1),sites(c,2)) = 1;

% save('adj_mat_f2f.edge','adj_mat_f','-ascii');
% save('adj_mat_f2p.edge','adj_mat_fp','-ascii');
% save('adj_mat_p2p.edge','adj_mat_p','-ascii');

%%
f_loc_s = squeeze(sum(plvs.*mask_f,[1 2]))/NactiveF;
f_loc_r = squeeze(sum(plvr.*mask_f,[1 2]))/NactiveF;

f_loc_diff = f_loc_s - f_loc_r;


p_loc_s = squeeze(sum(plvs.*mask_p,[1 2]))/NactiveP;
p_loc_r = squeeze(sum(plvr.*mask_p,[1 2]))/NactiveP;

p_loc_diff = p_loc_s - p_loc_r;


fp_loc_s = squeeze(sum(plvs.*mask_fp,[1 2]))/NactiveFP;
fp_loc_r = squeeze(sum(plvr.*mask_fp,[1 2]))/NactiveFP;

fp_loc_diff = fp_loc_s - fp_loc_r;


res_loc_s = squeeze(sum(plvs.*mask_res,[1 2]))/Nactiveres;
res_loc_r = squeeze(sum(plvr.*mask_res,[1 2]))/Nactiveres;

res_loc_diff = res_loc_s - res_loc_r;


%% bar plots
star_loc = 0.7;

x = f_loc_s;
y = f_loc_r;

colors = cbrewer('qual', 'Paired', 10);

% Plot error bars
figure
bar(1, mean(x), 'FaceColor',  colors(6, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on
bar(2, mean(y), 'FaceColor',  colors(5, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
% Plot samples
scatter(0.8:0.4/length(x):1.2-0.4/length(x), x, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(6, : ),...
        'LineWidth',1, 'SizeData',5)
scatter(1.8:0.4/length(y):2.2-0.4/length(y), y, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(5, : ),...
        'LineWidth',1, 'SizeData',5)

% Test of significance
% [~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
[~, pval, ci, stats] = ttest(x, y)
mysigstar(gca, [1 2], star_loc, pval);

% Plot lines
fp = plot([0.8:0.4/length(x):1.2-0.4/length(x); 1.8:0.4/length(y):2.2-0.4/length(y)],...
     [x y]', 'Color', [0.5, 0.5, 0.5, 0.1]);

set(h(1), 'marker', 'none'); 
set(gca, 'xtick', [1 1.5 2], 'xticklabel', {'Stimulus', 'F2F', 'Rest',},'xlim', [0.5 2.5],'XTickLabelRotation',45);
ylabel('PLV'); xlabel('Status');
ylim([0 1])

% set(gcf,'Units','Inches');
% set(gcf,'PaperPosition',[0 0 2 3],'PaperUnits','Inches','PaperSize',[2, 3])
% fname = res_path+sprintf('PLV_F2F_barplot');
% print(gcf, '-painters', '-dpdf', '-r300', fname+'.pdf');
% saveas(gcf, fname+'.fig')


%% correlation in same plot
figure

ger = ges_diff(5,:); %gamma entrainment respnse

colors = cbrewer('qual', 'Paired', 12);


mdl_f = fitlm(ger, f_loc_diff);
h = plot(ger, mdl_f.Fitted, 'Color', colors(10,:));
hold on
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
scatter(ger, f_loc_diff, 20, colors(10,:), 'filled')


mdl_fp = fitlm(ger, fp_loc_diff);
h = plot(ger, mdl_fp.Fitted, 'Color', colors(4,:));
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
scatter(ger, fp_loc_diff, 20, colors(4,:), 'filled')


mdl_p = fitlm(ger, p_loc_diff);
h = plot(ger, mdl_p.Fitted, 'Color', colors(8,:));
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
scatter(ger, p_loc_diff, 20, colors(8,:), 'filled')

% scatter(ges_mean, res_loc_diff, 20, 'g', 'filled')
% mdl_res = fitlm(ges_mean, res_loc_diff);
% h = plot(ges_mean, mdl_res.Fitted, 'g');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';

legend("F2F","F2P", "P2P", 'Location' ,'northwest')
xlabel(sprintf('Gamma entrainment response'))
ylabel(sprintf('PLV_{stimulus} - PLV_{rest}'))

% set(gcf,'Units','Inches');
% set(gcf,'PaperPosition',[0 0 4 3],'PaperUnits','Inches','PaperSize',[4, 3])
% fname = res_path+sprintf('PLVVsGER_all');
% print(gcf, '-painters', '-dpdf', '-r300', fname+'.pdf');
% saveas(gcf, fname+'.fig')

%% correlation on seperate plot
figure
param = fp_loc_diff;
ger = ges_diff(5,:); %gamma entrainment respnse

% scatter(ges_diff(site,:), param, 50, Label_cat, 'filled')
gscatter(ger, param, Label(gr_inc))


text(ger, param, Pa_code(gr_inc))
% text(ges_diff(site,:), param, Label)

% xlabel(sprintf('GES_{%d,stim.} - GES_{%d,rest}', site, site))
xlabel(sprintf('GES_{meanS,stim.} - GES_{meanS,rest}'))
[rho,pp] = corr(ger', param)

title(sprintf('Corr. Coeff. = %.3f', rho))
ylabel('$\bf{\Sigma PLV_{P2P,stim.} - \Sigma PLV_{P2P,rest}}$','Interpreter','latex')

% fname = out_path+'P2PPLV_Diff_45sites_corr';
% print(gcf, '-painters', '-dpng', '-r300', fname+'.png');
% saveas(gcf, fname+'.fig')



%% statistics
% pvalue of each connection
pval= [];
H = [];
pval_vec = [];
cnt=0;
for ich=1:bch_num-1
    for jch=ich+1:bch_num
      cnt = cnt + 1;  
      x = squeeze(plvs(ich, jch, :));
      y = squeeze(plvr(ich, jch, :));
      
      [H(ich,jch),pval(ich,jch)] = ttest(x,y,'Alpha',0.05);
      [H(jch,ich),pval(jch,ich)] = ttest(x,y,'Alpha',0.05);
      
      pval_vec(cnt, 1)= pval(ich,jch);
      pval_vec(cnt, 2)= ich;
      pval_vec(cnt, 3)= jch;
      
    end
end

%% Benjamini–Hochberg
[vec_sort,Idx] = sort(pval_vec(:,1));
vec_sort(:,2) = pval_vec(Idx,2);
vec_sort(:,3) = pval_vec(Idx,3);

nPair = (bch_num*(bch_num-1))/2;

figure
plot(1:nPair, vec_sort(:,1));
hold on
plot(1:nPair, 0.05*[1:nPair]/nPair);

idxs = find(vec_sort(:,1)' >= 0.05*[1:nPair]/nPair);
sig_Idx = min(idxs);

vec_sort(1:sig_Idx-1,4) = 1;

Hbenj=[];
for i=1:sig_Idx-1
    Hbenj(vec_sort(i,2),vec_sort(i,3))=1;
    Hbenj(vec_sort(i,3),vec_sort(i,2))=1;
end

%% bar plot
star_loc = 0.9;
% x = squeeze(plvs(27,28,:));
% y = squeeze(plvr(27,28,:));

x = squeeze(plvs(3,4,:));
y = squeeze(plvr(3,4,:));

colors = cbrewer('seq', 'Blues', 9);
% Plot error bars
figure
bar(1, mean(x), 'FaceColor',  colors(8, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on
bar(2, mean(y), 'FaceColor',  colors(6, : ), 'EdgeColor', 'none', 'BarWidth', 0.6);
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
% Plot samples
scatter(0.8:0.4/length(x):1.2-0.4/length(x), x, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(8, : ),...
        'LineWidth',1, 'SizeData',10)
scatter(1.8:0.4/length(y):2.2-0.4/length(y), y, 'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',colors(6, : ),...
        'LineWidth',1, 'SizeData',10)
% Plot lines
fp = plot([0.8:0.4/length(x):1.2-0.4/length(x); 1.8:0.4/length(y):2.2-0.4/length(y)],...
     [x y]', 'Color', [0.5, 0.5, 0.5, 0.5]);



% Test of significance
% [~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
[~, pval, ci, stats] = ttest(x, y)
mysigstar(gca, [1 2], star_loc, pval);


set(h(1), 'marker', 'none'); 
set(gca, 'xtick', [1 2], 'xticklabel', {'Axial sites.', 'Lateral sites',},'xlim', [0.5 2.5],'XTickLabelRotation',45);
% ylabel('Gamma entrainment score'); xlabel('Status');
ylabel('$\bf{\Sigma PLV_{stim.} - \Sigma PLV_{rest}}$','Interpreter','latex')

% 

%% EVD?!!

kappas = [];
kappar = [];
figure
for i=1:Pa_num
    [Us(:,:,i),Ss(:,:,i),Vs(:,:,i)] = svd(plvs(:, :, i)./norm(plvs(:,:,i),'fro'));
%     [Us(:,:,i),Ss(:,:,i),Vs(:,:,i)] = svd(plvs(:, :, i));
%     kappas(i) = cond(plvs(:,:,i));
    temp = diag(Ss(:,:,i));
    temp(temp==0)=[];
    kappas(i) = max(temp(:))/min(temp(:));
    plot(diag(Ss(:,:,i)))

    hold on
end

for i=1:Pa_num
    [Ur(:,:,i),Sr(:,:,i),Vr(:,:,i)] = svd(plvr(:, :, i)./norm(plvr(:,:,i),'fro'));
%     [Ur(:,:,i),Sr(:,:,i),Vr(:,:,i)] = svd(plvr(:, :, i));
%     kappar(i) = cond(plvr(:,:,i));
    temp = diag(Sr(:,:,i));
    temp(temp==0)=[];
    kappar(i) = max(temp(:))/min(temp(:));
    plot(diag(Sr(:,:,i)))

    hold on
end

%% maximum eigenvalue
S = Ss./Sr;
site = 4;

figure
scatter(ges_diff(site,:), S(1,1,:))
text(ges_diff(site,:), S(1,1,:), string(1:Pa_num))

xlabel(sprintf('GES_{Fz,stim.} - GES_{Fz,rest}'))
ylabel('$\bf{\lambda_{1,stim.} / \lambda_{1,rest}}$','Interpreter','latex')
[rho,~] = corr(ges_diff(site,:)', squeeze(S(1,1,:)));

title(sprintf('Corr. Coeff. = %.3f', rho))

%% conditon number
% 
site = 24;
kappa = kappas - kappar;

figure
% scatter(ges_diff(site,:), kappa)
gscatter(ges_diff(site,:), kappa, Label)
text(ges_diff(site,:), kappa, string(MMSE))

xlabel(sprintf('GES_{Pz-O2,stim.} - GES_{Pz-O2,rest}'))
ylabel('$\bf{\kappa_{stim.} - \kappa_{rest}}$','Interpreter','latex')
[rho,~] = corr(ges_diff(site,:)', kappa');

title(sprintf('Corr. Coeff. = %.3f', rho))

%% Mask?!!
% % fronto-parietal long distance

% frontal = [1 3 4 15 16];
% parietal = [11 12 14 27 28];
frontal = [22 26];
parietal = [23 25];
% frontal = [15 16];
% parietal = [27 28];
% frontal = [ 17:19 22:24];
% parietal = [ 19:21 24:26];

mask_long=zeros(bch_num, bch_num);
for ich=1:length(frontal)
    for jch=1:length(parietal)
        mask_long(frontal(ich), parietal(jch)) = 1;
        mask_long(parietal(jch), frontal(ich)) = 1;
    end
end
figure
heatmap(mask_long)


% % fronto-parietal local distance
mask_local=zeros(bch_num, bch_num);
for ich=1:length(frontal)
    for jch=ich+1:length(frontal)
        mask_local(frontal(ich), frontal(jch)) = 1;
        mask_local(frontal(jch), frontal(ich)) = 1;
    end
end
for ich=1:length(parietal)
    for jch=ich+1:length(parietal)
        mask_local(parietal(ich), parietal(jch)) = 1;
        mask_local(parietal(jch), parietal(ich)) = 1;
    end
end
figure
heatmap(mask_local)



%%
long_lens = squeeze(sum(plvs.*mask_long,[1 2]));
long_lenr = squeeze(sum(plvr.*mask_long,[1 2]));

long_lenDiff = long_lens - long_lenr;
long_lenNorm = long_lenDiff./long_lenr;

local_lens = squeeze(sum(plvs.*mask_local,[1 2]));
local_lenr = squeeze(sum(plvr.*mask_local,[1 2]));

local_lenDiff = local_lens - local_lenr;
local_lenNorm = local_lenDiff./local_lenr;

% laterals = squeeze(sum(plvs(1:14,1:14,:), [1 2]));
% lateralr = squeeze(sum(plvr(1:14,1:14,:), [1 2]));
% lateral_Diff = laterals - lateralr;
% 
% axials = squeeze(sum(plvs(15:end,15:end,:), [1 2]));
% axialr = squeeze(sum(plvr(15:end,15:end,:), [1 2]));
% axial_Diff = axials - axialr;

figure
plot(long_lens)
hold on
plot(long_lenr)

figure
plot(local_lens)
hold on
plot(local_lenr)

% figure
% plot(laterals)
% hold on
% plot(lateralr)
% 
% figure
% plot(axials)
% hold on
% plot(axialr)
%%
figure
tool = "long";
param = eval(tool+"_lenDiff");
site = 24;

% scatter(ges_diff(site,:), param, 50, Label_cat, 'filled')
gscatter(ges_diff(site,:), param, Label)

text(ges_diff(site,:), param, string(1:Pa_num))
% text(ges_diff(site,:), param, Label)

xlabel(sprintf('GES_{%d,stim.} - GES_{%d,rest}', site, site))
[rho,~] = corr(ges_diff(site,:)', param);

title(sprintf('%s\nCorr. Coeff. = %.3f\n %s - %s', tool, rho, strjoin(string(frontal)), strjoin(string(parietal))))
ylabel('$\bf{\Sigma PLV_{stim.} - \Sigma PLV_{rest}}$','Interpreter','latex')


% fname = out_path+sprintf('SigmaPLVDiff_%s_%s - %s_VsGES',tool,strjoin(string(frontal)),strjoin(string(parietal)))+string(site);
% print(gcf, '-painters', '-dpng', '-r300', fname+'.png');
% saveas(gcf, fname+'.fig')

%%
site = 24;
figure

scatter(ges_diff(site,:), long_lenDiff, 20, 'k', 'filled')
mdl_long = fitlm(ges_diff(site,:), long_lenDiff);
hold on
h = plot(ges_diff(site,:), mdl_long.Fitted, 'k');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


scatter(ges_diff(site,:), local_lenDiff, 20, 'b', 'filled')
mdl_local = fitlm(ges_diff(site,:), local_lenDiff);
h = plot(ges_diff(site,:), mdl_local.Fitted, 'b');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

% scatter(ges_diff(site,:), Clocal_lenDiff, 20, 'g', 'filled')
% mdl_Clocal = fitlm(ges_diff(5,:), Clocal_lenDiff);
% h = plot(ges_diff(site,:), mdl_Clocal.Fitted, 'g');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';

legend("long sites","local sites")

% fname = out_path+'Sigma_PLV_Diff_4sites';
% print(gcf, '-painters', '-dpng', '-r300', fname+'.png');
% saveas(gcf, fname+'.fig')

%% new idea, far dipoles
f_loc_s = squeeze(plvs(1,2,:));
fp_loc_s = squeeze(plvs(3,4,:));
p_loc_s = squeeze(plvs(5,6,:));

f_loc_r = squeeze(plvr(1,2,:));
fp_loc_r = squeeze(plvr(3,4,:));
p_loc_r = squeeze(plvr(5,6,:));

f_loc_diff = f_loc_s - f_loc_r;
fp_loc_diff = fp_loc_s - fp_loc_r;
p_loc_diff = p_loc_s - p_loc_r;

%%
figure
param = p_loc_diff;
ges_mean = mean(ges_diff,1);
site = 5;

% scatter(ges_diff(site,:), param, 50, Label_cat, 'filled')
gscatter(ges_mean, param, Label)

text(ges_mean, param, string(1:Pa_num))
% text(ges_diff(site,:), param, Label)

% xlabel(sprintf('GES_{%d,stim.} - GES_{%d,rest}', site, site))
xlabel(sprintf('GES_{meanS,stim.} - GES_{meanS,rest}'))
[rho,~] = corr(ges_mean', param);

title(sprintf('Corr. Coeff. = %.3f', rho))
ylabel('$\bf{\Sigma PLV_{P2P,stim.} - \Sigma PLV_{P2P,rest}}$','Interpreter','latex')

% fname = out_path+'P2PPLV_Diff_45sites_corr';
% print(gcf, '-painters', '-dpng', '-r300', fname+'.png');
% saveas(gcf, fname+'.fig')

%%
figure

scatter(ges_mean, f_loc_diff, 20, 'k', 'filled')
mdl_long = fitlm(ges_mean, f_loc_diff);
hold on
h = plot(ges_mean, mdl_long.Fitted, 'k');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

scatter(ges_mean, fp_loc_diff, 20, 'b', 'filled')
mdl_local = fitlm(ges_mean, fp_loc_diff);
h = plot(ges_mean, mdl_local.Fitted, 'b');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

scatter(ges_mean, p_loc_diff, 20, 'r', 'filled')
mdl_local = fitlm(ges_mean, p_loc_diff);
h = plot(ges_mean, mdl_local.Fitted, 'r');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


legend("F3-Fp1, F4-Fp2","F3-P3, F4-P4", "P3-O1, P4-O2")
xlabel(sprintf('GES_{meanS,stim.} - GES_{meanS,rest}'))

ylabel('$\bf{PLV_{stim.} - PLV_{rest}}$','Interpreter','latex')

% fname = out_path+'PLV_Diff_6sites';
% print(gcf, '-painters', '-dpng', '-r300', fname+'.png');
% saveas(gcf, fname+'.fig')
























%% 
figure
for site=[20 25]
% site = 24;
    phase = angle(squeeze(fData_allWin_allPa(f40_idx, site, stim_idx(:), 31)));

    
    polarhistogram(phase,15)
    hold on
end
%% delta-phi = 180 due to volume conduction?
s1_idx = 20;
s2_idx = 25;
phi = angle(fData_allWin_allPa(f40_idx, :, :, :));

delta_phis = squeeze(rad2deg(((phi(:, s1_idx, stim_idx(:), :) - phi(:, s2_idx, stim_idx(:), :)))));
delta_phir = squeeze(rad2deg(((phi(:, s1_idx, rest_idx(:), :) - phi(:, s2_idx, rest_idx(:), :)))));

% delta_phis = squeeze(rad2deg(wrapTo2Pi((phi(:, s1_idx, stim_idx(:), :) ))));
% delta_phir = squeeze(rad2deg(wrapTo2Pi((phi(:, s1_idx, rest_idx(:), :) ))));

idx_40 = 5; % w=2 =>5, w=20 =>41
site = 4;
temp = (abs(fData_allWin_allPa(freq_idx, site, :, :)));
temp_z = zscore(temp,0,1);
gess = squeeze(temp_z(idx_40,:,stim_idx(:),:));
gesr = squeeze(temp_z(idx_40,:,rest_idx(:),:));


%%
figure
p = 31;

scatter(gess(:,p), delta_phis(:,p) , 15, 'filled')
hold on
scatter(gesr(:,p), delta_phir(:,p) , 15, 'filled')
