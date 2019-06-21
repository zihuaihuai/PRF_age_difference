


baseDir='../../../RawData/';
baseDir_disk = '\\DiskStation\HCP\RawData\';

filepath.MRacq=[baseDir,'Physio/',subject,'_',task,'/phys.mat'];
filepath.phys=[baseDir,'Physio/',subject,'_',task,'/phys_sum.mat'];

load(filepath.MRacq,'Fs','trig','TR','time');
load(filepath.phys);
Fs=400; Ts=1/Fs; 

volDel=40;  

nVol=1200;
volAlreadyDel=1200-nVol;

ind_BOLD=find(trig==1);
trig(ind_BOLD(1:volDel))=0; 
ind_BOLD=find(trig==1);

Fs_10=10; Ts_10=1/Fs_10; time_10 = time(1):Ts_10:time(end);

timeMR=time(find(trig==1));
ind_BOLD_10=zeros(length(timeMR),1);
for i=1:length(timeMR)   
    [~,ind]=min(abs(time_10-timeMR(i)));   %% consider number of volumes before task
    ind_BOLD_10(i)=ind;
end
