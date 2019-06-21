clear
path = '/home/enning/Desktop/Diskstation/HCP/RawData_PC1_Michalis/';
% choose session
sess_names = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
sess = 1;
volDel=40;  
nVol=1200;
volAlreadyDel=1200-nVol;

load([path,'subject_list.mat']);
subject_list = [subject_list_R1; subject_list_R2;subject_list_R3;subject_list_R4];
HCP_table = importfile('HCP_Young_behavioral.csv', 2, 1207);

dir3 = [path,'Physio/',char(subject_list(1)),'_',char(sess_names(sess)),'/','Phys.mat'];
load(dir3, 'trig','time');
ind_BOLD=find(trig==1);
trig(ind_BOLD(1:volDel))=0; 
ind_BOLD=find(trig==1);
Fs_10=10; Ts_10=1/Fs_10; time_10 = time(1):Ts_10:time(end);
timeMR=time(find(trig==1));
ind_BOLD_10=zeros(length(timeMR),1);
for j=1:length(timeMR)   
    [~,ind]=min(abs(time_10-timeMR(j)));   %% consider number of volumes before task
    ind_BOLD_10(j)=ind;
end

for i = 1:length(subject_list)
    if mod(i,5) ==0
        i
    end
    dir1 = [path,'Atlas/',char(subject_list(i)),'_',char(sess_names(sess)),'/','TissueBasedRegressors_1199.mat'];
    dir2 = [path,'Physio/',char(subject_list(i)),'_',char(sess_names(sess)),'/','Phys_sum.mat'];
    data(i).name = char(subject_list(i));
    load(dir1,'WB','GM');
    data(i).WB = WB.MA(41:end); data(i).GM = GM.MA(41:end);
    
    load(dir2,'HRV','resp_10');
    data(i).HR = HRV; data(i).resp_10 = resp_10;
    idx = find(HCP_table.Subject == subject_list(i));
    data(i).timeMR = timeMR; data(i).time10 = time_10; data(i).ind_BOLD_10 = ind_BOLD_10;
    
    if HCP_table.Gender(idx) == 'M'
        data(i).gender = 0;
    else
        data(i).gender = 1;
    end
    
    if HCP_table.Age(idx) == '22-25'
        data(i).age = 0;
    elseif HCP_table.Age(idx) == '26-30'
        data(i).age = 1;
    else 
        data(i).age = 2;
    end
end

length(data)
save('data.mat','data')