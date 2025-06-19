%% Targeted memory reactivation elicits temporally compressed reactivation linked to spindles
% Mahmoud E. A. Abdellahi, Martyna Rakowska, Matthias S. Treder & Penelope A. Lewis

% Code by Mahmoud Abdellahi, 2024
% emails: abdellahime@cardiff.ac.uk   m.eid@fci-cu.edu.eg

% Feel free to use the code to reproduce the results or for reference in
% your own work, but please cite the work: Targeted memory reactivation elicits temporally compressed reactivation linked to spindles
% (link to the article with details will be provided here upon publication).

%% general parameters 
% pathAppend is the path of main folder that should contain the following folders:
% 1. code
% 2. data
% 3. fieldtrip-20190419
% (please change pathAppend accordingly)
pathAppend = 'D:/work/jittered exp/Data and code for_Targeted memory reactivation elicits temporally compressed reactivation linked to spindles';
cd([pathAppend '/code']);
analysis = 'behavioural_analyses'; % please change according to the analysis you would like to run: erp_tf, classification, high_spindle_pw, low_spindle_pw, temporal_compression, behavioural_analyses

%% configuring plots, and participants ids, importing fieldtrip default functions
addpath([pathAppend '\fieldtrip-20190419'])
ft_defaults
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontSize',16)
set(0,'DefaultAxesTitleFontWeight','normal');
sbj = [3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,20,21,22,23,24,26,2,3,5,6,7,8,9,11,13,14,15,16,17,19,20,21,22,23,24,26,27,28,29,30,31,32,33];
%% classification and all analyses
% this is the main loop body that loops on different participants and
% executes the required analysis, by the end of that loop we save the
% output and then each of following blocks represent a visualisation and
% statistical analysis for a specific analysis

% general parameters settings to load the correct files6
trn_ds = 'img';
tst_ds = 'sleep';
if strcmp(trn_ds,'img')==1, sleep_stage_trn=0; else trn_ds='sleep'; sleep_stage_trn=3; end
if strcmp(tst_ds,'img')==1, sleep_stage_tst=0; else tst_ds='sleep'; sleep_stage_tst=3; end
sequence = {'lv_feature_extractor','erp 100'}; % to run many sequences back to back put new seq. in rows so it will be run(s) x sequence
if strcmp(analysis,'temporal_compression')
    sequence = {'lv_feature_extractor','erp 5'}; 
end
[trn,tst,classification_result,fidelity_pt,fidelity_pt2,lens,lens_sleep,recurrence_distrib]=deal([]);
temp=[]; erp=[]; erp_trn=[];  spindle_likelihood=[];
vschance = 1; % vs chance or do we have control condition
if vschance==1, conditions=1; else, conditions=2; end
cleaned_path = [pathAppend '/data/part']; 
if ~exist('sequence','var') % if the sequence not provided then load from excel sheet
    table_specs = readtable([pwd '\lv_automation_runs.xlsx']);
    titles = table_specs.Properties.VariableNames;
    variables = table_specs.Variables;
    run_id=[];
    for i=1:size(variables,1), if strcmp(char(variables(i,end)),'done')==0, run_id=[run_id ; i]; end, end
    sequence = variables(run_id,2:end);
    fprintf(['\n lv: Applying pipeline sequence(s): \n']); disp(sequence);
    record=repmat({' '},size(variables,1),1); record(run_id,1)={'done'};
    table_specs = [table_specs record];
    writetable(table_specs, [pwd '\lv_automation_runs.xlsx']); 
end
axtimeAcc=[];

% main loop on participants 
for run=1:size(sequence,1)
    for cond=1:conditions 
        for nn=1:numel(sbj)
            if strcmp(analysis,'behavioural_analyses')
                break;
            end
            if nn>=22 
                cleaned_path = [pathAppend '/data/MRI_part']; 
            end
            fprintf(['\n lv: Working on pipeline for sbj: ' num2str(sbj(nn)) ' \n']);
            trn.data = lv_load([cleaned_path num2str(sbj(nn)) '_' trn_ds '_manual_cleaned_N' num2str(sleep_stage_trn)],'trial');
            tst.data = lv_load([cleaned_path num2str(sbj(nn)) '_' tst_ds '_manual_cleaned_N' num2str(sleep_stage_tst)],'trial');
            
            % cutting the trials with respect to their variable lengths
            temp = tst.data.trial_lens;
            [~,~,common_all] = intersect([tst.data.trialinfo(:,2) tst.data.trialinfo(:,5)],[temp(:,2) temp(:,5)], 'rows','stable');
            temp = temp(common_all,:);
            if ~isequal(tst.data.trialinfo(:,1:5),temp(:,1:5)), error('lv_ mismatch between records'); end
            tst.data.trialinfo = [tst.data.trialinfo temp(:,6)]; 
            % rejecting bad trials and fixing long ones
            temp = tst.data.trialinfo(:,6);
            bad_id = find(temp<2.5);
            if length(bad_id)>0, bad_id, warning('trial <2.5'); end % not real error, but just a check
            tst.data.trial(bad_id,:,:)=[]; tst.data.trialinfo(bad_id,:)=[]; tst.data.sampleinfo(bad_id,:)=[];
            temp = tst.data.trialinfo(:,6);
            long_id = find(temp>3.5);
            tst.data.trialinfo(long_id,6)=3.5; clear temp;
               

            if strcmp(analysis,'erp_tf')
                % time frequency analysis
                hold_data = tst.data;
                cfg=[]; cfg.latency = [-0.6 3.8]; tst.data = ft_selectdata(cfg, tst.data);
                TF_struct= do_tf(tst.data, [-0.3 -0.1], [1 30]);
                TF_struct.trial = squeeze(mean(TF_struct.trial, 1)); % average of channels
                tf_pw(nn,:,:) = TF_struct.trial; % participant x frequency x time
    
                % ERP analysis
                tst.data = hold_data;
                baseline_id = nearest(tst.data.time,-0.2):nearest(tst.data.time,0); % time
    
                baseline = (mean(tst.data.trial(:,:,baseline_id),3));
                baseline = repmat(baseline, 1,1,size(tst.data.trial,3));
                erp(nn,:) = mean(mean(tst.data.trial - baseline, 2),1);
    
                continue;
            end    

            % extracting features
            cfg=[];
            cfg.data=trn.data;
            cfg.sequence = sequence;
            result_trn = lv_build_pipeline(cfg); % for trn

            cfg.data=tst.data;
            result_tst = lv_build_pipeline(cfg); % for tst
        
            % checking class labels 
            classes=unique(result_trn.trialinfo(:,1)), if length(classes)~=4, error('lv: classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([1])) ,1)  = 1;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([2])) ,1)  = 2;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([3])) ,1)  = 3;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([4])) ,1)  = 4;
            classes=unique(result_tst.trialinfo(:,1)), if length(classes)~=4, error('lv: classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([1])) ,1)  = 1;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([2])) ,1)  = 2;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([3])) ,1)  = 3;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([4])) ,1)  = 4;
            
            % cutting trials to 2.5s if performing temporal compression
            % analysis (as this is the minimum length and more than that is not consistent among trials
            % as this is jittered)
            if strcmp(analysis,'temporal_compression')==1
                cfg=[];
                cfg.latency=[0 2.5];
            else
                cfg=[];
                cfg.latency=[0 3.5];
            end
            result_tst=ft_selectdata(cfg,result_tst);
            % cutting trials for wake motor imagery
            if strcmp(trn_ds,'img')==1, cfg.latency=[0 1.15]; result_trn=ft_selectdata(cfg,result_trn);  end
            if strcmp(tst_ds,'img')==1, cfg.latency=[0 1.15]; result_tst=ft_selectdata(cfg,result_tst); end
            sz = size(result_tst.trial);

            % variable length nans, to append trials' lengths so that they
            % match
            idx = cell2mat(arrayfun(@(x) (nearest(result_tst.time,x)), result_tst.trialinfo(:,6),'Un',0));
            for i=1:size(result_tst.trial,1), result_tst.trial(i,:,idx(i)+1:end) = nan; end
            hold_trialinfo = result_tst.trialinfo;  

            if (strcmp(analysis,'high_spindle_pw')==1) || (strcmp(analysis,'low_spindle_pw')==1)
                % filtering in spindle band and performing hilbert
                % transform to get power 
                cfg_preprocessing                 = [];
                cfg_preprocessing.bpfilter        = 'yes';
                cfg_preprocessing.bpfreq          = [11 16];
                data_bp= ft_preprocessing(cfg_preprocessing, tst.data  );
                data_bp_bl = data_bp;
                cfg=[]; cfg.latency=[0 2.5]; cfg.channel = lower({'cz'});
                data_bp=ft_selectdata(cfg,data_bp);
    
                hilbert_dat=[];
                for h=1:size(data_bp.trial,1) %trls
                    for j=1:size(data_bp.trial,2) %ch
                        hilbert_dat(h,j,:) = hilbert( squeeze(data_bp.trial(h,j,:)) ); % time should be in the first dimension.
                    end
                end 
                hilbert_dat = squeeze(hilbert_dat);
                hilbert_dat = abs(hilbert_dat).^2; % abs to get the magnitude and squaring that to get power 
                D1 = squeeze(median(hilbert_dat,2));
                % median split based on power
                if strcmp(analysis,'high_spindle_pw')
                    keep_id = find(D1 > median(D1,1)); 
                end
                if strcmp(analysis,'low_spindle_pw')
                    keep_id = find(D1 < median(D1,1)); 
                end
                result_tst.trial = result_tst.trial(keep_id,:,:);
                result_tst.trialinfo = result_tst.trialinfo(keep_id,:);
            end    
 
            result_trn.trial = single(result_trn.trial); result_tst.trial = single(result_tst.trial);
            
            % calculate pca using sleep data
            cfg=[];
            cfg.data = result_tst;
            cfg.method = 'pca';
            cfg.step = 'calculate';
            cfg.centered = 1; % data is centered with demean in preprocessing
            comp = lv_component_analysis(cfg); %cfg=rmfield(cfg,'data2');
            % transform both sleep and wake data using the PCs calculated
            % using sleep
            cfg.eigVects = comp.eigVects;
            cfg.chosen = zeros(length(comp.eigVals),1);
            id = find(cumsum(comp.eigVals)> 0.95); % 95% of variance
            cfg.chosen(1:id(1))=1;
            cfg.step = 'transform'; 
            result_tst.trial = lv_component_analysis(cfg);
            cfg.data = result_trn;
            result_trn.trial = lv_component_analysis(cfg);

            % temporal compression with classification
            if strcmp(analysis,'temporal_compression')==1
                cfg=[]; cfg.data_long=result_tst; cfg.data_short=result_trn; cfg.method='compression_wilson_jittered';
                cfg.centered = 0; % for PCA because it will perform centering when calculating pca inside lv_post_classification
                warning('for compression analysis no smoothing should be done and the length of trials should be the same');
                temp = lv_post_classification(cfg);
                compressions{1,nn} = temp.compressionRatio(temp.score==max(temp.score));
                compression_scores(:,nn) = temp.score;
                continue;
    
            end
            % organising time points
            temp2 = repmat(result_tst.trialinfo(:,1),1,size(result_tst.trial,3)); result_tst.trialinfo=[];
            result_tst.trialinfo(:,1)=reshape(temp2,[],1);
            temp = permute(result_tst.trial, [2 1 3]);
            result_tst.trial = (reshape(temp, size(temp,1),[]))';
            id = find(isnan(sum(result_tst.trial,2)));  keep_id = find(~isnan(sum(result_tst.trial,2))); temp=nan(size(result_tst.trial,1),1);
            result_tst.trial(id,:,:)=[]; result_tst.trialinfo(id,:)=[]; 
 
            % sleep to wake classification (sleep for training and wake for testing)
            result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            cfg=[]; cfg.method = 'timextime';
            cfg.classifier_type = {'lda'}; % linear discriminant analysis classifier
            cfg.perf_measure = 'acc'; % evaluate with accuracy
            cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1); cfg.tst.trial = single(cfg.tst.trial);
            cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1);
            cfg.folds=nan; cfg.do_parallel=1;
            cfg.trn.trial = single(cfg.trn.trial);
            cfg.weights = ones(size(cfg.trn.trial,1),1);
            fidelity_pt(nn,:) = lv_classify(cfg); 
            continue; 
        end
    end
    if (strcmp(analysis,'high_spindle_pw')==1)
        fidelity_pt_sp = fidelity_pt;
        save fidelity_pt_sp fidelity_pt_sp;
    end
    if (strcmp(analysis,'low_spindle_pw')==1)
        fidelity_pt_sp_low = fidelity_pt;
        save fidelity_pt_sp_low fidelity_pt_sp_low;
    end
    if strcmp(analysis,'temporal_compression')==1
        new_compression_scores = compression_scores;
        new_compressionRatio = temp.compressionRatio;
        save new_compression_scores new_compression_scores;
        save new_compressionRatio new_compressionRatio;
    end
    if strcmp(analysis,'classification')==1
        save fidelity_pt fidelity_pt
        % figure 2 b
        load fidelity_pt fidelity_pt
        res = lv_pretty_errorbar(result_trn.time, fidelity_pt, (fidelity_pt*0)+0.25, 1);
        fig = gcf;
        set(fig, 'NumberTitle', 'off', 'Name', 'Figure 2b');
    end
end
%% ERP and TF analysis, figure 2a
if strcmp(analysis,'erp_tf')
    save erp erp
    save tf_pw tf_pw
    save TF_struct TF_struct
end
if strcmp(analysis,'erp_tf')
    figure, 
    load tf_pw tf_pw
    load TF_struct TF_struct
    id1 = nearest(TF_struct.time,0):nearest(TF_struct.time,2.5); % not after 2.5 because of the jittering
    TF_struct.time = TF_struct.time(id1);
    id2 = nearest(TF_struct.freq,5):nearest(TF_struct.freq,30); 
    TF_struct.freq = TF_struct.freq(id2);
    temp = tf_pw(:,id2,id1);
    b = imagesc(TF_struct.time, TF_struct.freq, squeeze(mean(temp,1)).*100 ); set(gca,'YDir','normal') 
    xlabel('Time (sec.)', 'Interpreter','none');
    ylabel('Frequency Hz', 'Interpreter','none');
    h = colorbar; title(h,'Power');
    caxis([-25 25])
    % ERP analysis
    load erp erp
    hold on,
    yyaxis right
    id1 = nearest(tst.data.time,0):nearest(tst.data.time,2.5);
    plot(tst.data.time(id1), mean(erp(:,id1),1), ...
        'black-')  
    set(gca,'YColor','k');
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 2b');
end

%% high and low sigma classification .. requires running both high_spindle_pw and fidelity_pt_sp_low conditions 
% as this will compare the results 
if strcmp(analysis,'low_spindle_pw') || strcmp(analysis,'high_spindle_pw')
    load fidelity_pt_sp_low fidelity_pt_sp_low;
    load fidelity_pt_sp fidelity_pt_sp; 
    % figure 4a
    res = lv_pretty_errorbar(result_trn.time, fidelity_pt_sp, (fidelity_pt_sp*0)+0.25, 1);
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 4a');
    % figure 4b
    figure, 
    res = lv_pretty_errorbar(result_trn.time, fidelity_pt_sp, fidelity_pt_sp_low, 1);
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 4b');
end
%% behavioural analyses
%%  cued vs. un-cued which is the aggregation from different sessions
if strcmp(analysis,'behavioural_analyses')
    pwd
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,3); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,3);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    s3r = isnan(nanmean(re,2)) | isnan(nanmean(nre,2)) | isnan(nanmean(re_random,2)) | isnan(nanmean(nre_random,2));
    
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,4); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,4);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    rppnt = find(isnan(nanmean(re,2)) | isnan(nanmean(nre,2)) | isnan(nanmean(re_random,2)) | isnan(nanmean(nre_random,2)) | s3r);
    
    % preMRI (more participants)
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,2:4); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,2:4);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    
    re22 = re(:,1:24); nre22 = nre(:,1:24); re_random22 = re_random(:,1:2); nre_random22 = nre_random(:,1:2); 
    for i=1:size(re22,1), q1=prctile(re22(i,:),5); re22(i, re22(i,:)>q1)=nan; end
    for i=1:size(nre22,1), q1=prctile(nre22(i,:),5); nre22(i, nre22(i,:)>q1)=nan; end
    
    re33 = re(:,25:48); nre33 = nre(:,25:48); re_random33 = re_random(:,3:4); nre_random33 = nre_random(:,3:4);
    for i=1:size(re33,1), q1=prctile(re33(i,:),5); re33(i, re33(i,:)>q1)=nan; end
    for i=1:size(nre33,1), q1=prctile(nre33(i,:),5); nre33(i, nre33(i,:)>q1)=nan; end
    
    re44 = re(:,49:end); nre44 = nre(:,49:end); re_random44 = re_random(:,5:6); nre_random44 = nre_random(:,5:6);
    for i=1:size(re44,1), q1=prctile(re44(i,:),5); re44(i, re44(i,:)>q1)=nan; end
    for i=1:size(nre44,1), q1=prctile(nre44(i,:),5); nre44(i, nre44(i,:)>q1)=nan; end
    
    re22 = nanmedian(re22,2); re33 = nanmedian(re33,2); re44 = nanmedian(re44,2); nre22 = nanmedian(nre22,2); nre33 = nanmedian(nre33,2); nre44 = nanmedian(nre44,2);
    
    re1 = [re22  re33  re44]; nre1 = [nre22  nre33  nre44];
    re1(rppnt,:)=[]; nre1(rppnt,:)=[];
    
    % pre-sleep session1
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,1); 
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,1); 
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    for i=1:size(re,1), q1=prctile(re(i,:),5); re(i, re(i,:)>q1)=nan; end
    for i=1:size(nre,1), q1=prctile(nre(i,:),5); nre(i, nre(i,:)>q1)=nan; end 
    
    
    re = nanmedian(re,2); nre = nanmedian(nre,2);

    % % correcting with random for pre session
    % re = re - nanmedian(nanmin(re_random,[],2),2); % effect - random(higher rt) so lower/negative is better 
    % nre = nre - nanmedian(nanmin(nre_random,[],2),2);
    % % for post sessions 2 3 4
    % nre_rand_post = [nanmedian(nanmin(nre_random22,[],2),2),nanmedian(nanmin(nre_random33,[],2),2),nanmedian(nanmin(nre_random44,[],2),2)]; nre_rand_post(rppnt,:)=[];
    % re_rand_post = [nanmedian(nanmin(re_random22,[],2),2),nanmedian(nanmin(re_random33,[],2),2),nanmedian(nanmin(re_random44,[],2),2)]; re_rand_post(rppnt,:)=[];
    % re1 = nanmedian(re1,2) - nanmedian(re_rand_post,2);
    % nre1 = nanmedian(nre1,2) - nanmedian(nre_rand_post,2);

    re(rppnt,:)=[]; nre(rppnt,:)=[]; re_random(rppnt,:)=[]; nre_random(rppnt,:)=[];
    re1 = bsxfun(@minus,nanmedian(re,2), re1); % pre - post
    nre1 = bsxfun(@minus,nanmedian(nre,2), nre1);
    % figure 2c
    all_stats = lv_pretty_errorbar(nanmedian(re1 ,2),nanmedian(nre1 ,2),' ',' ');
    hold_post = nanmedian(re1,2); 
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 2c');
    
    
    pre = nanmedian(re,2); 
    winlen=14; preTime=[];
    r = 1:length(pre)-winlen; 
    range = [r' r'+winlen];
    R = nanmedian(re1,2); NR = nanmedian(nre1,2);
 
    [val, id] = sort(nanmedian(re,2)); Rdat=[]; NRdat=[];
    % axtime with improvement
    for i=1:size(range,1)
        Rdat = [Rdat R( id(range(i,1):range(i,2)) )]; NRdat = [NRdat NR( id(range(i,1):range(i,2)) )];
        preTime(i,1) = mean( val( range(i,1):range(i,2) ) );
    end
    % figure, lv_pretty_errorbar(round(preTime'), Rdat-NRdat, (Rdat-NRdat).*0, 0); % difference
    % h=gca; h.XTickLabelRotation = 90;
    % xlabel('Encoding reaction time (ms)');
    % ylabel('Improvement (cued - uncued)');
    % fig = gcf;
    % set(fig, 'NumberTitle', 'off', 'Name', 'Figure 4a');
    %% partial correlation between improvement after sleep and classification figure 2d
    pwd
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,3); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,3);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    s3r = isnan(nanmean(re,2)) | isnan(nanmean(nre,2)) | isnan(nanmean(re_random,2)) | isnan(nanmean(nre_random,2));
    
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,4); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,4);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    rppnt = find(isnan(nanmean(re,2)) | isnan(nanmean(nre,2)) | isnan(nanmean(re_random,2)) | isnan(nanmean(nre_random,2)) | s3r);
    
    % preMRI (more participants)
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,2:4); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,2:4);  % MRI
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    
    re22 = re(:,1:24); nre22 = nre(:,1:24); re_random22 = re_random(:,1:2); nre_random22 = nre_random(:,1:2); 
    for i=1:size(re22,1), q1=prctile(re22(i,:),5); re22(i, re22(i,:)>q1)=nan; end
    for i=1:size(nre22,1), q1=prctile(nre22(i,:),5); nre22(i, nre22(i,:)>q1)=nan; end
    
    re33 = re(:,25:48); nre33 = nre(:,25:48); re_random33 = re_random(:,3:4); nre_random33 = nre_random(:,3:4);
    for i=1:size(re33,1), q1=prctile(re33(i,:),5); re33(i, re33(i,:)>q1)=nan; end
    for i=1:size(nre33,1), q1=prctile(nre33(i,:),5); nre33(i, nre33(i,:)>q1)=nan; end
    
    re44 = re(:,49:end); nre44 = nre(:,49:end); re_random44 = re_random(:,5:6); nre_random44 = nre_random(:,5:6);
    for i=1:size(re44,1), q1=prctile(re44(i,:),5); re44(i, re44(i,:)>q1)=nan; end
    for i=1:size(nre44,1), q1=prctile(nre44(i,:),5); nre44(i, nre44(i,:)>q1)=nan; end
    re22 = nanmedian(re22,2); re33 = nanmedian(re33,2); re44 = nanmedian(re44,2); nre22 = nanmedian(nre22,2); nre33 = nanmedian(nre33,2); nre44 = nanmedian(nre44,2);
    re1 = [re22]; nre1 = [nre22];
    
    % pre-sleep session1
    [re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,1); 
    [re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,1); 
    re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
    for i=1:size(re,1), q1=prctile(re(i,:),5); re(i, re(i,:)>q1)=nan; end
    for i=1:size(nre,1), q1=prctile(nre(i,:),5); nre(i, nre(i,:)>q1)=nan; end 
    
    re = nanmedian(re,2); nre = nanmedian(nre,2);
    
    % figure 2d
    % correlating with classification .. for session 2 because after sleep
    load fidelity_pt fidelity_pt;
    % the correlation without presleep correction
    stats = lv_vec2corr(max(fidelity_pt,[],2), (nanmedian(nre1 ,2)-nanmedian(re1 ,2)),'Classification performance (CCR)','Improvement after sleep (ms)'); % correlating with classification
    % partial correlation
    [rho,pval] = partialcorr(max(fidelity_pt,[],2), nanmedian(nre1,2) - nanmedian(re1,2) ,nanmedian(nre,2) - nanmedian(re,2), 'type', 'Spearman')
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 2d');
end

%% temporal compression
% significant scaling factors
if strcmp(analysis,'temporal_compression')==1
    load new_compression_scores new_compression_scores
    load new_compressionRatio new_compressionRatio
    lv_pretty_errorbar(1:size(new_compression_scores,1), new_compression_scores', (new_compression_scores'.*0)+0.25,0);
    new_compressionRatio = round(new_compressionRatio, 2);
    new_compressionRatio(new_compressionRatio<1) = 1./new_compressionRatio(new_compressionRatio<1); new_compressionRatio = round(new_compressionRatio, 1);
    xticklabels(new_compressionRatio);
    xtickangle(90)
    fig = gcf;
    set(fig, 'NumberTitle', 'off', 'Name', 'Figure 3');
    textHandles = findall(fig, 'Type', 'text'); 
    delete(textHandles);
    xlabel('compression/dilation');
    ylabel('CCR');
end
%% helping functions
%% TF analysis
function TFdat = do_tf(dat, baseline, frequencies) % takes 3d in .trial (trls_ch_time) and returns (ch_freq_time)

cfg              = [];
cfg.output       = 'pow';

cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:0.5:30; %linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1)));
cfg.t_ftimwin    = 5./cfg.foi;  % 5 cycles as a minimum to describe the frequency well
cfg.toi          = dat.time; % .time for max resolution 
cfg.pad          ='nextpow2'; % rounds the maximum trial length up to the next power of 2
cfg.keeptrials = 'no';
TFdat = ft_freqanalysis(cfg, dat);


if ~isempty(baseline) && baseline(1)~=0
    cfg              = [];
    cfg.baseline     = [baseline(1) baseline(2)];
    cfg.baselinetype = 'relchange';
    [TFdat] = ft_freqbaseline(cfg, TFdat); % ch x freq x time
end

TFdat.trial = TFdat.powspctrm;

end
function [re,nre,re_random,nre_random] = extract_blocks(ds_name,sbj,sessions) 
% takes name of the dataset and returns the blocks of all sessions
[re,nre,re_random,nre_random] = deal([]); 
for i=1:24, R_pre{i,1} = ['R_pre_' num2str(i)]; NR_pre{i,1} = ['NR_pre_' num2str(i)]; end % blocks names
for i=1:2, R_random{i,1} = ['R_random_pre_' num2str(i)]; NR_random{i,1} = ['NR_random_pre_' num2str(i)]; end
if strcmp(ds_name,'myDat_s')==1, sbj=sbj(1:21); else, sbj=sbj(22:end); end
for session=sessions
    load behav_lbls behav_lbls;
    var = [ds_name num2str(session)];
    load (var); myDat = eval(var);
    id = find(ismember(myDat(:,1), sbj)); % ids of sbj in excel
    dat = myDat(id,3:end); behav_lbls = behav_lbls(3:end);

    re = [re dat(:, ismember(behav_lbls,R_pre))];% reactivated seq. blocks
    nre = [nre dat(:, ismember(behav_lbls,NR_pre))];
    re_random = [re_random dat(:, ismember(behav_lbls,R_random))];
    nre_random = [nre_random dat(:, ismember(behav_lbls,NR_random))];
end
end


