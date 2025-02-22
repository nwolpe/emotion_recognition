%% script to analyse Cam-CAN Ekman Face data

%% In this code, I've used two methods (psychometric fit + Signal Detection Theory) 
% to calculate 4 summary metrics for people's responses in the emotion recognition task.
% I have not looked at RT
clear

%% read raw table for analysis

% read table of big raw data: structure is one row for one observation for
% one subject - so multiple rows per subject
data_all = readtable('~/Documents/Emotion_studies/EkmanTask/data_all.csv');

% converting RTs to ms from seconds
data_all.RTs = data_all.RTs/1000;

% find out indices for each subject 
G = findgroups(data_all.subject);


%% fit psychometric curve (to calculate threshold and width) and calculate basic SDT measures of d_prime (sensitivity) and criterion (bias)

% to fit a psychometric curve I used psignifit:
% https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/neuronale-informationsverarbeitung/research/software/psignifit/

% sigmoid_type = {'gauss', 'logistic', 'weibull'};

cd ~/Documents/Emotion_studies/EkmanTask/
% these are required parameters for the psychometric fitting using psignifit:
options             = struct;
% options.sigmoidName = 'logistic'; 
% options.sigmoidName = 'gauss'; % for normal scaling
options.sigmoidName = 'weibull'; % for log scaling
% the choice design was n Alternative Forced Choice - i.e., people chose 1 option between 6
options.expType     = 'nAFC';
% there were 6 options in the experiment on each trial
options.expN = 6;
% this handles whether you make all parameters free (fitted) or fix some
options.fixedPars = NaN(5,1); %% no fixed params - all free parameters
% fixing "guess rate" parameter gamma to zero which sets the lower
% asymptote to zero
options.fixedPars(4) = 0;

% to run psignifit you need to add it to the path:
addpath ~/Documents/MatlabScipts/psignifit-master/

% initialise variables
emotion = {'fea', 'sur', 'sad', 'hap', 'dis', 'ang'};
threshold = nan(length(unique(G)),length(emotion));
width = nan(length(unique(G)),length(emotion));
accuracy = nan(length(unique(G)),length(emotion));
rt = nan(length(unique(G)),length(emotion));
bic = nan(length(unique(G)),length(emotion));

% for each subject (unique g) fit the psychometric function 

for iSub = 1:length(unique(G))
    
    clearvars data_sub
    
    % there needs to be 25 trials per emotion
    
    data_sub = data_all(G == iSub,:);
    
    % only using 'test' phase
    data_sub = data_sub(strncmp(data_sub.trial_type,'t',1),:);
    
    % remove trials with extreme RTs
    data_sub ( data_sub.RTs < 0.2 | data_sub.RTs > 10,:) = [];
    
    % display which subject is being fit:
    text2display = ['Fitting subject ', num2str(unique(data_sub.subject))];
    disp(text2display)
    
    % emotion order: 
    % 1-fear, 2-surprise, 
    % 3-sadness, 4-happiness,
    % 5-disgust, 6-anger
%     
% colors = [
%     0.4, 0.6, 0.8;  % Subtle blue
%     0.6, 0.4, 0.4;  % Muted red
%     0.5, 0.7, 0.5;  % Soft green
%     0.7, 0.6, 0.3;  % Subtle yellow-brown
%     0.5, 0.5, 0.7;  % Soft purple
%     0.6, 0.5, 0.4   % Muted brown
% ];
%     clf
    for iEmotion = 1:length(emotion)
        
        data2fit = data2fit_gen_nozero(data_sub, emotion{iEmotion});
        
        % actual model fitting function:
        modelfit = psignifit(data2fit,options);
        
        % plot for illustration in paper
        
%         plotOptions.dataColor = colors(iEmotion, :);  % Set marker color
%         plotOptions.lineColor = colors(iEmotion, :);  % Set line color
% 
%         h= plotPsych(modelfit, plotOptions);
%         hold on
%         plotHandles(iEmotion) = h;
        
%         legendLabels{iEmotion} = emotion{iEmotion};  % Add the label for each emotion

        % Set color for each line
%         set(h, 'Color', colors(iEmotion, :));
        
        % Overlay data points with the same color as the line
%         plot(data2fit(:, 1), data2fit(:, 2) ./ data2fit(:, 3), 'o', 'MarkerFaceColor', colors(iEmotion, :));
        
        bic(iSub, iEmotion) = modelfit.deviance + 3 * log(sum(data2fit(:,3))); % 3 free params
        
        % comment and uncomment below to see individual fits
%         clf;
%         plotPsych(modelfit);
%         pause
        
        % now extract the two calculated parameters:
        threshold(iSub, iEmotion) = getThreshold(modelfit,0.5);
        width(iSub, iEmotion) = modelfit.Fit(2);
        % calculate overall accuracy for 90 and 70% stimulus
        accuracy(iSub, iEmotion) = sum(data2fit(end-1:end,2) / sum(data2fit(end-1:end,3)));
        rt(iSub, iEmotion) = compute_RT_median(data_sub,emotion{iEmotion});
        
        
        
        clearvars modelfit data2fit
    end
    


end


%% bic comparison
bic_gauss = 5.7908e+04;
bic_logistic = 5.8241e+04;
bic_weibull = 5.7073e+04;

bic_values = [bic_gauss-bic_weibull, bic_logistic-bic_weibull, bic_weibull-bic_weibull];
model_labels = {'Gaussian', 'Logistic', 'Weibull'};
figure;
set(gcf,'color','white');
bar(bic_values);

% Add labels for the x-axis
set(gca, 'XTickLabel', model_labels);

% Add axis labels and title
xlabel('Model Function', 'FontSize', 16);
ylabel('Relative BIC Value', 'FontSize', 16);
% title('BIC Comparison for Different Models');

% Adjust font size for readability
set(gca, 'FontSize', 14);
print('Fig2A', '-djpeg', '-r300');


%% plot individual's fit 


plotPsych(modelfit);
set(gcf, 'Color', 'White')
xlabel('Stimulus level (%)')

fig = gcf;
ax = gca;
ax.YLabel.String = 'Probability of choosing sad';
ax.XLabel.String = 'Percentage sad in picture';
xticks([10, 50, 100]); % Add more values if needed

xticklabels({'10', '50','100'}); % Change this to match your desired labels


% Save the figure as a high-resolution PNG
print('high_res_image', '-dpng', '-r300');


%% summarise task and demographic data for R analyses

% emotion = {'fea', 'sur', 'sad', 'hap', 'dis', 'ang'};
data_all = readtable('~/Documents/Emotion_studies/EkmanTask/data_all.csv');

% load required demographic data
load('~/Documents/Emotion_studies/demog_cc700.mat')
load('~/Documents/Emotion_studies/EkmanTask/Participation.mat') % dates of sessions
load('~/Documents/Emotion_studies/EkmanTask/hi_dataset.mat')
load('~/Documents/Emotion_studies/covariates_additional.mat')
additional = load('~/Documents/CamCAN//hi_additional.mat');

sum(cellfun(@strcmp,hi.CCID, demog.ccid)); %sanity check for same ccid order
demog.interview_date = hi.date_iv;

% ccid is participant's id
ccid_Ekman = unique(data_all.subject);

% find out ccid in demographic data
ccid_mat = cell2mat(demog.ccid); ccid_700 = str2num(ccid_mat(:,3:end));

% now intersect them
[~, i700, iEkman] = intersect(ccid_700, ccid_Ekman);
% to calculate people's age when doing the task, extract the data for
% performing the task
EkmanDate = datetime(Participation.EkmanEmotionRecognition(i700), 'Format', "dd/MM/yyyy");
% calculate age at time of testing
age = between(datetime(demog.dob(i700), 'Format', 'MM/yyyy'), datetime(EkmanDate(iEkman), 'Format', 'dd/MM/yyyy'), 'years');
age_num = split(age, 'years');

% save Benton faces data - this is a basic task for identifying faces (no
% emotion) so can serve as a control
BentonFacesScore = covariates_additional.BentonFaces_totalScore(i700);
BentonFacesScore = cellfun(@str2double, BentonFacesScore);

% Cattell is a fluid intelligence test
CattellScore = covariates_additional.Cattell_totalScore(i700);
CattellScore = cellfun(@str2double, CattellScore);

handedness = demog.handedness(i700);

% save sex and education
sex = demog.sex(i700); sex_num = strncmp(sex, 'F',1);
education = demog.qualifications(i700);
education_num(strncmp(education, '0',1))=0; % 0 none of the above
education_num(strncmp(education, '4',1))=1; % 4 cse/gcse
education_num(strncmp(education, '3',1))=1; % 3 o levels

education_num(strncmp(education, '2',1))=2; % 2 a levels
education_num(strncmp(education, '5',1))=3; % 5 nvq / national vocational qualification, NVQ or HND or HNC or equivalent 
education_num(strncmp(education, '6',1))=3; % 6 other professional, Other professional qualifications e.g.: nursing, teaching
education_num(strncmp(education, '1',1))=3; % 1 college or university
education_num(strncmp(education, '8',1))=nan; % 8 no answer

acer_orientation = additional.hi.attention_orientation(i700);
acer_memory = additional.hi.memory(i700);
acer_fluency = additional.hi.fluencies(i700);
acer_language = additional.hi.language(i700);
acer_visuospatial = additional.hi.visuospatial(i700);
acer_total = additional.hi.acer(i700);

hads_depression = demog.hads_depression(i700);
hads_depression(hads_depression>100) = NaN; % means they did not respond in the interview

hads_anxiety = demog.hads_anxiety(i700);
hads_anxiety(hads_anxiety>100) = NaN;
hads_total = demog.hads_anxiety(i700) + demog.hads_depression(i700);

% save ACER (cognitive/dementia test) score
acer =  demog.acer(i700);


% emotion = {'fea', 'sur', 'sad', 'hap', 'dis', 'ang'};
glm_table = table(hads_depression, threshold(:,1), threshold(:,2),threshold(:,3), ...
    threshold(:,4), threshold(:,5), threshold(:,6), age_num, age_num.^2,CattellScore, BentonFacesScore,...
    rt(:,1),rt(:,2),rt(:,3),rt(:,4),rt(:,5),rt(:,6), acer_fluency, acer_language, acer_memory,...
    acer_orientation, acer_total, acer_visuospatial,score(:,1),apathy,sex_num, education_num',...
    hads_anxiety,handedness,score(:,2),...
    'VariableNames', {'depression', 'threshold_fear', 'threshold_surprise', 'threshold_sadness', ...
    'threshold_happiness', 'threshold_disgust', 'threshold_anger', 'age', 'age_squared', 'cattell','Benton',...
    'rt1','rt2','rt3','rt4','rt5','rt6','acer_fluency', 'acer_language', 'acer_memory',...
    'acer_orientation', 'acer_total', 'acer_visuospatial','factor1','enthusiasm', 'sex', 'education',...
    'anxiety','handedness', 'factor2'});


writetable(glm_table, '/Users/nwolpe/Documents/Emotion_studies/EkmanTask/glm_table.csv')

% readtable('/Users/nwolpe/Documents/Emotion_studies/EkmanTask/glm_table.csv')
