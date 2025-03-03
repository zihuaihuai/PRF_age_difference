clear all;
%% train the scan specific function.
tic
load("data/data.mat");
num_sessions = 1;
num_feat = 14;
sess_names = {'rest1_LR', 'rest1_RL', 'rest2_LR','rest2_RL'};
parfor i = 1:length(data)
    [R_crf,R_rrf,R_prf, curve_feature,CRF,RRF] = PRF_sc(i, data);
    out(i).crf = R_crf;
    out(i).rrf = R_rrf;
    out(i).prf = R_prf;
    out(i).para = curve_feature;
    out(i).curve_c = CRF;
    out(i).curve_r = RRF;
end
toc
save('result/result.mat', 'out');
out_m = out(find(data.gender == 0));
out_f = out(find(data.gender == 1));
%% plot the crf and rrf with error bar
c_m = reshape([out_m.curve_c],[],size(out_m,2))';
r_m = reshape([out_m.curve_r],[],size(out_m,2))';
c_f = reshape([out_f.curve_c],[],size(out_f,2))';
r_f = reshape([out_f.curve_r],[],size(out_f,2))';
load curve_stand.mat;

set(gcf,'Position',get(0,'ScreenSize'));
ax1 = subplot(2,1,1);
c_y_mean = mean(c_m);
c_y_std = std(c_m)/sqrt(size(r_m,1));
shadedErrorBar(0:0.1:60,c_y_mean,c_y_std, 'lineprops','b');hold on 
c_o_mean = mean(c_f);
c_o_std = std(c_f)/sqrt(size(r_f,1));
shadedErrorBar(0:0.1:60,c_o_mean,c_o_std);hold on
shadedErrorBar(0:0.1:60,CRF_stan,zeros(601,1),'lineprops','r'); 
legend('young', 'old', 'standard');
title('cardiac response function for all');
hold off

ax2 = subplot(2,1,2);
r_y_mean = mean(r_m);
r_y_std = std(r_m)/sqrt(size(r_m,1));
shadedErrorBar(0:0.1:60,r_y_mean,r_y_std,'lineprops','b');hold on 
r_o_mean = mean(r_f);
r_o_std = std(r_f)/sqrt(size(r_f,1));
shadedErrorBar(0:0.1:60,r_o_mean,r_o_std); hold on
shadedErrorBar(0:0.1:60,RRF_stan,zeros(601,1),'lineprops','r');
legend('young', 'old', 'standard');
title('respiration response function for all');
hold off
saveas(gcf,'jpgs/MeanShadedErr.jpg')
clear gcf

c_m_all = c_m;
r_m_all = r_m;
c_f_all = c_f;
r_f_all = r_f;
% the plots for 3 sessons separately
set(gcf,'Position',get(0,'ScreenSize'));

for i = 1: num_sessions
    s = (i-1) * 12 + 1;
    e = 12 * i;
    c_m = c_m_all(s:e,:);
    r_m = r_m_all(s:e,:);
    s = (i - 1) * 9 + 1;
    e = 9 * i;
    c_f = c_f_all(s:e,:);
    r_f = r_f_all(s:e,:);
    ax1 = subplot(3,2,2*i-1);
    c_y_mean = mean(c_m);
    c_y_std = std(c_m)/sqrt(size(r_m,1));
    shadedErrorBar(0:0.1:60,c_y_mean,c_y_std, 'lineprops','b');hold on 
    c_o_mean = mean(c_f);
    c_o_std = std(c_f)/sqrt(size(r_f,1));
    shadedErrorBar(0:0.1:60,c_o_mean,c_o_std);hold on
    shadedErrorBar(0:0.1:60,CRF_stan,zeros(601,1),'lineprops','r'); 
    legend('young', 'old', 'standard');
    title(['cardiac response function',sess_names(i)]);

    ax2 = subplot(3,2,2*i);
    r_y_mean = mean(r_m);
    r_y_std = std(r_m)/sqrt(size(r_m,1));
    shadedErrorBar(0:0.1:60,r_y_mean,r_y_std,'lineprops','b');hold on 
    r_o_mean = mean(r_f);
    r_o_std = std(r_f)/sqrt(size(r_f,1));
    shadedErrorBar(0:0.1:60,r_o_mean,r_o_std); hold on
    shadedErrorBar(0:0.1:60,RRF_stan,zeros(601,1),'lineprops','r');
    legend('young', 'old', 'standard');
    title(['respiration response function',sess_names(i)]); 

end
saveas(gcf,'jpgs/ShadedError.jpg');
clear gcf
set(gcf,'Position',get(0,'ScreenSize'));
for i = 1:num_sessions
    s = (i-1) * 12 + 1;
    e = 12 * i;
    c_m = c_m_all(s:e,:);
    r_m = r_m_all(s:e,:);
    s = (i - 1) * 9 + 1;
    e = 9 * i;
    c_f = c_f_all(s:e,:);
    r_f = r_f_all(s:e,:);    
    ax1 = subplot(3,4,4*i-3)    
    for j = 1:12
        plot(c_m(j,:));hold on;
    end
    plot(CRF_stan,'LineWidth',1.5,'Color','k');
    title([sess_names(i),'CRF in young group'])

    ax2 = subplot(3,4,4*i-2)    
    for j = 1:12
        plot(r_m(j,:));hold on;
    end
    plot(RRF_stan,'LineWidth',1.5,'Color','k');
    title([sess_names(i),' RRF in young group'])


    ax3 = subplot(3,4,4*i-1)    
    for j = 1:9
        plot(c_f(j,:));hold on;
    end
    plot(CRF_stan,'LineWidth',1.5,'Color','k');
    title([sess_names(i),'CRF in old group'])

    ax4 = subplot(3,4,4*i)    
    for j = 1:9
        plot(r_f(j,:));hold on;
    end
    plot(RRF_stan,'LineWidth',1.5,'Color','k');
    title([sess_names(i),'RRF in old group'])
end
saveas(gcf,'jpgs/AllLine.jpg')
%% two pairs test of the correlation and curve features
names = {'prf_r' ,'c1_amplitude', 'c2_amplitude', 'time_to_c1',...
    'time_to_c2','c1_width', 'c2_width', 'cardic_shape', 'r1_amplitude',...
    'r2_amplitude', 'time_to_r1', 'time_to_r2','r1_width', 'r2_width', 'respiratory_shape'};

para_m = reshape([out_m.para],14,[])';
para_f = reshape([out_f.para],14,[])';
feat_m = array2table(para_m, "VariableNames",names(2:end));
feat_f = array2table(para_f, "VariableNames",names(2:end));
p(1) = ranksum([out_m.prf],[out_f.prf]);
means(1,1) = mean([out_m.prf]);
means(2,1) = mean([out_f.prf]);
for i = 1:num_feat 
    p(i+1) = ranksum(para_m(:,i),para_f(:,i));
    means(1,i+1) = mean(para_m(:,i));
    means(2,i+1) = mean(para_f(:,i));
end
means(3, :) = p;
result_table = array2table(means,'RowNames',{'male','female','p-val'},'VariableNames',names);
pair = result_table(:,find(p<0.05));
save("result.mat");

%% anova test 
%check the correlation 
nova = array2table(zeros(3,15));
nova.Properties.VariableNames = names;
nova.Properties.RowNames = {'Age','Sess','Age & Sess'};

[a, b] = anova_pre([out_m.prf,out_f.prf],num_sessions,12);
[tbl,~] = simple_mixed_anova(a,b ,{'Session'},{'Age'});
tbl = tbl([2,4,5],5); 
nova(:,1) = tbl;
clear a b tbl;
for i = 1:num_feat
    [a, b] = anova_pre([para_m(:,i);para_f(:,i)], num_sessions,12);
    [tbl,~] = simple_mixed_anova(a,b ,{'Session'},{'Age'});
    tbl = tbl([2,4,5],5); 
    nova(:,i+1) = tbl;
    clear a b tbl;
end
stat3 = [result_table; nova]
save 3sses.mat stat3
