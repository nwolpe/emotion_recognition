data = readtable('data_emotion_recognition2.csv');
load('motion_data.mat')
id=table2array(data(:,1));

PATH='aamod_smooth_00001' %folder with fMRI data
eval(['cd ' PATH])

for ii=1:length(id)
inc(ii)=length(eval(['dir(''CC' num2str(id(ii)) '/Rest/fcMRI/biHIPfinal/SPM.mat'')']));
end
incl_emo=find(inc==1);
education=data.education(incl_emo);
sex_factor = data.sex(incl_emo);
sex = strncmp(sex_factor, 'F',1)+0;
depression = data.depression(incl_emo);
depression(find(isnan(depression))) = nanmean(depression);
depression=zscore(depression);
benton = data.Benton(incl_emo);
hand = data.handedness(incl_emo);
hand(isnan(hand)) = nanmean(hand);
hand=zscore(hand);
education(isnan(education)) = nanmean(education);
benton(isnan(benton)) = nanmean(benton);
benton=zscore(benton);


age = (data.age);
factor1 = (data.factor1);

age = zscore(data.age(incl_emo));
factor1 = zscore(factor1(incl_emo));
intTerm=(factor1-mean(factor1)).*(age-mean(age));
intTermOrth=orthog(intTerm, [factor1 age]);

matchmotion=id(incl_emo);
for kk=1:length(matchmotion);
   matchmotionind(kk)=find((Final_CBU_CCID(:,1))==matchmotion(kk));
end
motion=resting_rms(matchmotionind);
motion(isnan(motion))=nanmean(motion);
motion=zscore(motion);

eval(['cd ' PATH])

ccidinc=id(incl_emo);
ROI='bihp_highpassonly'


for ii=1:length(ccidinc)
eval(['copyfile CC' num2str(ccidinc(ii)) '/Rest/fcMRI/' ROI '/beta_0002.nii ../ANALYSIS/REST/' ROI '/data/CC' num2str(ccidinc(ii)) '.nii']);
end

%%% ROI
%[b, ~, residuals] = regress(clus1, [ones(size(covar,1),1), covar]);