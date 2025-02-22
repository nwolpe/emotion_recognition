addpath /imaging/rowe/users/nw03/software/spm12_fil_r7219/
addpath(genpath('/imaging/rowe/users/nw03/software/robustcorrtool-master/'))
cd '/imaging/rowe/users/nw03/CamCAN/emotion_recognition/'
clear

%%

% create scan paths
scans = dir('/imaging/rowe/users/nw03/CamCAN/GM_release004/smoothed_8fwhm/s*.nii');
for i=1:length(scans)
    scanPath{i} = [scans(i).folder, '/', scans(i).name, ',1'];
end

% get CC IDs for scan paths
for iScanID = 1:length(scanPath)
    allScanID(iScanID) = str2double(scanPath{iScanID}(end-11:end-6));
end

% load behavioural params

data = readtable('/imaging/rowe/users/nw03/CamCAN/emotion_recognition/data_emotion_recognition2.csv');
ccid = data.ccid;

% intersect those who have scans with those that have emo task
[commmon, incl_scans, incl_emo] = intersect(allScanID', ccid);


%% prepare variables for VBM

TIV=nan(length(allScanID),1);
TGM=nan(length(allScanID),1);
for i=1:length(allScanID)
    subID=sprintf('%d.mat', allScanID(i));
    file2load = sprintf('/imaging/camcan/cc700/mri/pipeline/release004/data/aamod_structuralstats_00001/CC%d/structurals/structuralstats.mat', allScanID(i));
    load(file2load);
    TIV(i) = S.TIV.mm3_spm12;
    TGM(i) = S.parts.mm3_spm12(1);
end


education=data.education(incl_emo);
sex_factor = data.sex(incl_emo);
sex = strncmp(sex_factor, 'F',1)+0;
depression = data.depression(incl_emo);
depression(find(isnan(depression))) = nanmean(depression);
benton = data.Benton(incl_emo);

benton(isnan(benton)) = nanmean(benton);

data4pca = zscore([data.threshold_1, data.threshold_3, data.threshold_4, data.threshold_6]);
[coeff,score,~, ~, explained] = pca(data4pca);

age = data.age;

factor1 = data.factor1;
intTerm=(factor1-mean(factor1)).*(age-mean(age));

age = data.age(incl_emo);
factor1 = factor1(incl_emo);
intTerm = intTerm(incl_emo);
intTermOrth=orthog(intTerm, [factor1 age]);

TIV4analysis = TIV(incl_scans);
scans4analysis = scanPath(incl_scans)';

dir2save= {'/imaging/rowe/users/nw03/CamCAN/emotion_recognition/'};


%% load movement data

headmotion = load('motion_data_Noham.mat');

headmotion_mat = [headmotion.Final_CBU_CCID(:,1), headmotion.resting_rms',...
    headmotion.movie_rms',headmotion.avtask_rms'];

[~, incl_sofar, incl_motion] = intersect(ccid(incl_emo),headmotion_mat);

motion = nanmean(headmotion_mat(incl_motion,2:end),2);

% replace nan with mean
motion(isnan(motion)) = nanmean(motion);


%% run VBM
clearvars matlabbatch

matlabbatch{1}.spm.stats.factorial_design.dir = dir2save;
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans4analysis;

matlabbatch{1}.spm.stats.factorial_design.cov(1).c = factor1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'factor1';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(2).c = intTerm;
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'factor1_age_int';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(3).c = age;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(4).c = sex;
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(5).c = benton;
matlabbatch{1}.spm.stats.factorial_design.cov(5).cname = 'benton';
matlabbatch{1}.spm.stats.factorial_design.cov(5).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(5).iCC = 1;
%
matlabbatch{1}.spm.stats.factorial_design.cov(6).c = motion;
matlabbatch{1}.spm.stats.factorial_design.cov(6).cname = 'motion';
matlabbatch{1}.spm.stats.factorial_design.cov(6).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(6).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % multiple covariates from file
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % no threshold masking
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % implicit mask
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/rowe/users/nw03/motorLearning/VBMall/GM_0.15bin.nii,1'}; % explicit mask

% global calculation - user
matlabbatch{1}.spm.stats.factorial_design.globalc.g_user.global_uval = TIV4analysis;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1; % overall grand mean scaling 1= no
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 3; % normalisation 3=ANCOVA, 2=proportional, 1=none

% now estimation

matlabbatch{2}.spm.stats.fmri_est.spmmat = {sprintf('%s/SPM.mat', dir2save{1})};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% now run contrast

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'int_pos';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'int_neg';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'factor1_pos';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'factor1_neg';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.delete = 0;


spm('defaults', 'PET');
spm_jobman('run', matlabbatch, {});
