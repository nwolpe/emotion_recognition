%%%seed2WB_firstevel

%%% create ROI and extract subject specific time course to subj/Rest/fcMRI
%%% Then:


PATH='aamod_smooth_00001' %folder with fMRI data
eval(['cd ' PATH])
 
warning off
addpath('/usr/local/spm/spm12');
addpath('/usr/local/spm/spm2/toolbox/FieldMap/');
addpath('/usr/local/nexus-tools/matlab');
bldPrefix = 'spm_bld';


%ROIname='biHIPfinal';
ROIname='bihp'
spm_defaults

tps = 261;
RT = 1.97;
mwd = PATH;
    HPC = 1/0.01;
    LPC = 1/0.2;
   K   = spm_dctmtx(tps,tps);
nHP = fix(2*(tps*RT)/HPC + 1);
nLP = fix(2*(tps*RT)/LPC + 1);
K   = K(:,[2:nHP nLP:tps]);      % Remove initial constant
K   = K(:,[2:nHP]);      % Remove initial constant
Nk  = size(K,2);
for ii=1:Nk
    bpssname{ii}=['highpass' num2str(ii)];
end
    



sub=dir('CC*')
for ii=1:length(sub)
    
    eval(['cd ' sub(ii).name '/Rest/']);
    cd globals
    eval(['load ' sub(ii).name '.mat']); 
    cd ../motion
    eval(['load ' sub(ii).name '.mat'])
    % despike
     SpikeMovRelThr = 5;     % 5 SDs of mean?
     SpikeLag = 1;
     M=motion;
         dM  = [zeros(1,6); diff(M,1,1)];    % First-order derivatives
         dM(:,4:6) = dM(:,4:6)*50;  % Approximate way of converting rotations to translations, assuming sphere radius 50mm (and rotations in radians) from Linda Geerligs
         cdM       = sum(abs(dM),2);
         rspk = find(cdM > (mean(cdM) + SpikeMovRelThr*std(cdM)));
spk = rspk; lspk = spk;
for q = 2:SpikeLag
    lspk = [lspk; spk + (q-1)];
end
spk = unique(lspk);

if ~isempty(spk)
    RSP = zeros(tps,length(spk));
    n = 0;
    for p = 1:length(spk)
        if spk(p) <= tps
            n=n+1;
            RSP(spk(p),n) = 1;
        end
    end 
    fprintf('%d unique spikes in total\n',length(spk))
    RSP = spm_en(RSP,0);
else
    RSP = [];
end
for pp=1:size(RSP,2);
     RSPname{pp}=['spike' num2str(ii)];
end
    % end of despike
      
    cd ../fcMRI
    eval(['TC=load(' '''' sub(ii).name '_' ROIname '.txt'');'])
    mkdir([ROIname '_highpassonly'])
    eval(['cd ' ROIname '_highpassonly']);
    if ~exist([mwd '/SPM.mat']),
        SPM.nscan          = 1;            % number of scans for each of nsess sessions
        SPM.xY.RT          = RT;                   % experiment TR in seconds
        SPM.xGX.iGXcalc    = 1;
        SPM.xX.K=1;       
    % design (conditions and user specified covariates)
    %---------------------------------------------------------------------------
    
    SPM.xX.X = [ones(tps,1) TC motion [0 0 0 0 0 0; diff(motion)] motion.^2 [0 0 0 0 0 0; diff(motion)].^2 GM_WM_CSF [0 0 0; diff(GM_WM_CSF)] GM_WM_CSF.^2 [0 0 0; diff(GM_WM_CSF)].^2 K RSP] 
    SPM.xX.name = [{'intercept', 'ROI' 'mov1' 'mov2' 'mov3' 'mov4' 'mov5' 'mov6' 'mov1dt' 'mov2dt' 'mov3dt' 'mov4dt' 'mov5dt' 'mov6dt' 'mov1s' 'mov2s' 'mov3s' 'mov4s' 'mov5s' 'mov6s' 'mov1dts' 'mov2dts' 'mov3dts' 'mov4dts' 'mov5dts' 'mov6dts' 'GM' 'WM' 'CSF' 'GMdt' 'WMdt' 'CSFdt' 'GMs' 'WMs' 'CSFs' 'GMdts' 'WMdts' 'CSFdts'} bpssname RSPname];
    SPM.xX.iH = 1;
    SPM.xX.iC = 2:(48+length(RSPname));;
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Beginning estimation \n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
      SPM.xsDes.Design='Multiple regression'
      SPM.xsDes.Global_calculation='omit'
      SPM.xsDes.Grand_mean_scaling='<no grand Mean scaling>'
      SPM.xsDes.Global_normalisation='<no global normalisation>'
  
    cd ../../
    tmp=dir('*wds*nii*')
    P=[mwd '/' sub(ii).name '/Rest/' tmp.name];
    eval(['cd fcMRI/' ROIname '_highpassonly'])
    SPM.xY.P           = char(P);
    v=spm_vol(P);
    SPM.xY.VY=v;
    SPM.swd=pwd;
   
    %Estimate parameters
    %===========================================================================
    SPM = spm_spm(SPM);
    eval(['cd ' mwd]);
    clear GM_WM_CSF motion P TC
        else
            fprintf(1,'Model for current subject already defined. Only adding contrasts !\n');
    end
              
end

