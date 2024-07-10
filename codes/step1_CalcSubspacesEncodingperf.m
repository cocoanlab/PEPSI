% All Dependencies (in README.md) should be installed and added in the path
clear;clc 
% basedir   = '/path/to/github/PEPSI';
figoutdir = fullfile(basedir, 'figures/step1_CalcSubspacesEncodingperf');
if ~exist(figoutdir, 'dir'), mkdir(figoutdir); end
NTWRK   = {'Cerebrum', 'CerebrumStriatumCerebellum'};
parcelnames = {'Visual', 'Somatomotor', 'dAttention', 'vAttention', ...
    'Limbic', 'Frontoparietal', 'Default'};
condName = {'highLV2', 'highLV3', 'highLV4', 'highLV5', ...
            'lowLV1' , 'lowLV2', 'lowLV3', 'lowLV4', ...
            'noLV2'  , 'noLV3',  'noLV4'};
Cues     = {'high', 'low', 'no'};
Stims    = {'LV1', 'LV2', 'LV3', 'LV4', 'LV5'};

% load data
load(fullfile(basedir, 'data', 'data_for_replication.mat'));

% key variables are ...
%   neurAvg.FIR 
%       : subject averaged FIR response. comprising whole voxel.
%   neurAvg.CuetBeta, neurAvg.StimtBeta 
%       : temporal encoding weights using neurAvg.FIR as Y, and "CueStimX" as X
%   CueStimX
%       : [Intercept, CueInfo, StimInfo]. Normalized. Condition ordered as in
%       variable "condName"
%   parcelIndx
%       : spatial index of each network
%   behvout
%       : mat size of #cond X #sub. pain reports. condition is ordered as
%       in variable "condName"
% Brain mask
%   obj = fmri_data(fullfile(basedir, 'data', 'template.nii'))

%% 1. Selection of PC.
pccrit      = 70;
cmaps       = cbrewer('qual', 'Dark2', 7);

for iN = 1:2
    % iN = 2. we used 7-network parcellation that encompasses
    % Cerebrum, Striatum, and also Cerebellum.
    % iN = 1. This one includes only Cerebrum. 
    for iiP = 1:numel(parcelnames)
        int_ntw = NTWRK{iN};
        int_par = parcelnames{iiP};
        int_idx = parcelIndx.(int_ntw).(int_par);
        cue_tBetas  = neurAvg.CuetBeta(int_idx, :)';
        stim_tBetas = neurAvg.StimtBeta(int_idx, :)';
        [CuePC, ~, ~, ~, CueExpl]  = pca(cue_tBetas); 
        [StimPC, ~, ~, ~, StimExpl] = pca(stim_tBetas);
        CueExpls_ntw(iN, iiP, :) = CueExpl;
        StimExpls_ntw(iN, iiP, :) = StimExpl;
    end
end

for iN = 1:2
    figure('Name', NTWRK{iN});
    cuePCexpls  = squeeze(CueExpls_ntw(iN, :, :));
    stimPCexpls = squeeze(StimExpls_ntw(iN, :, :));
    for i = 1:7
        subplot 211;
        p1(i) = plot(cuePCexpls(i, :), 'Color', cmaps(i, :), 'LineWidth', 1); hold on;
        numPC = min(find(cumsum(cuePCexpls(i, :)) >= pccrit));
        line([numPC numPC], [0, 80], 'Color', cmaps(i, :), 'LineStyle', '--', 'LineWidth', 1);
        xlabel('Number of PC');ylabel('Percentage Explained'); title('In Cue Subspace');
        xlim([1 20]);ylim([0 60])
        set(gca, 'TickDir', 'out', 'FontSize', 10, 'FontName', 'Helvetica', 'LineWidth', 1.5); box off;
        subplot 212;
        p2(i) = plot(stimPCexpls(i, :), 'Color', cmaps(i, :), 'LineWidth', 1); hold on;
        numPC = min(find(cumsum(stimPCexpls(i, :)) >= pccrit));
        line([numPC numPC], [0, 80], 'Color', cmaps(i, :), 'LineStyle', '--', 'LineWidth', 1);
        xlabel('Number of PC');ylabel('Percentage Explained'); title('In Stim Subspace');
        xlim([1 20]);ylim([0 60])
        set(gca, 'TickDir', 'out', 'FontSize', 10, 'FontName', 'Helvetica', 'LineWidth', 1.5); box off;
    end
end

% Dashed line is the number of PC explaining 70% PCs.
% Going to select FIXED number of PC. Since after PC exceeding number of 10,
% Percentage Explained is almost flat ... 
% Important thing is, all the PC selected doesn't make a significant
% differences in the results.

%% 2. Calculate Encoding performances in each subspace
mean_behvout = mean(behvout, 2);
pccrits = 20; % [10 15 20 25]; 
for iP = 1:numel(pccrits)
    pccrit = pccrits(iP);
    cue_effects = {};   stim_effects = {};
    for iN = 1:2
        for iiP = 1:numel(parcelnames)
            int_ntw = NTWRK{iN};
            int_par = parcelnames{iiP};
            int_idx = parcelIndx.(int_ntw).(int_par);
            cue_tBetas  = neurAvg.CuetBeta(int_idx, :)';
            stim_tBetas = neurAvg.StimtBeta(int_idx, :)';
            
            % cuePC and stimPC are Cue Subspace and Stim Subspace respectively.
            [CuePC, ~, ~, ~, CueExpls]  = pca(cue_tBetas); 
            [StimPC, ~, ~, ~, StimExpls] = pca(stim_tBetas);
            CuePC  = CuePC(:, 1:pccrit);
            StimPC = StimPC(:, 1:pccrit);
            pcstr  = sprintf('NumPC%d', pccrit);
            
            % Reshaping for projection to each subspace and then reshape to
            % original space.
            FIRs        = neurAvg.FIR(int_idx, :, :);
            [n_vox, n_time, n_cond] = size(FIRs);
            r_FIRs      = reshape(FIRs, n_vox, []);  
            r_FIRs0     = r_FIRs - mean(r_FIRs, 2);
            r_cuePrjs   = CuePC' * r_FIRs0; % Projected trajectory in cue subspace.
            r_stimPrjs  = StimPC'* r_FIRs0; % Projected trajectory in stim subspace.
            FIRs0       = reshape(r_FIRs0, n_vox, n_time, n_cond);
            cuePrj      = reshape(r_cuePrjs, [], n_time, n_cond);
            stimPrj     = reshape(r_stimPrjs, [], n_time, n_cond);
    
            % ------ Start Calculating Encoding performances ------
            for iT = 1:n_time
                FIRs0_t   = squeeze(FIRs0(:, iT, :));
                cuePrj_t  = squeeze(cuePrj(:, iT, :))';
                stimPrj_t = squeeze(stimPrj(:, iT, :))';
                
                % Single Vector based on maximal distance.
                [cuemaxidx, ~]  = find(squareform(pdist(cuePrj_t)) == max(pdist(cuePrj_t)));
                [stimmaxidx, ~] = find(squareform(pdist(stimPrj_t)) == max(pdist(stimPrj_t)));
                cue_diffs       = (cuePrj_t(cuemaxidx(1), :) - cuePrj_t(cuemaxidx(2), :));
                stim_diffs      = (stimPrj_t(stimmaxidx(1), :) - stimPrj_t(stimmaxidx(2), :));
                cue_normvec     = cue_diffs ./ norm(cue_diffs);
                stim_normvec    = stim_diffs ./ norm(stim_diffs);
                cue1dim         = cuePrj_t*cue_normvec';
                stim1dim        = stimPrj_t*stim_normvec';
                    % Now the variable "cue1dim" and "stim1dim" has
                    % distance information between conditions.
                        
                mdl_cue  = CueStimX(:, 2); % normalized [-1 0 1]
                mdl_stim = CueStimX(:, 3); % normalized [1 2 3 4 5]
                CueCue   = fitlm(mdl_cue, cue1dim); % function "fitlm" adds Intcercept as default.
                                                    % direction of the axis doesn't matter here
                StimCue  = fitlm(mdl_stim, cue1dim);
                StimStim = fitlm(mdl_stim, stim1dim);
                CueStim  = fitlm(mdl_cue, stim1dim);
                
                RSQs.(int_ntw).(int_par).XY_CueCue(iT, :)  = CueCue.Rsquared.Ordinary;
                RSQs.(int_ntw).(int_par).XY_StimCue(iT, :) = StimCue.Rsquared.Ordinary;
                RSQs.(int_ntw).(int_par).XY_CueStim(iT, :) = CueStim.Rsquared.Ordinary;
                RSQs.(int_ntw).(int_par).XY_StimStim(iT, :) = StimStim.Rsquared.Ordinary;
            end
            % ------ End Calculating Null Encoding performances ------
            fprintf('DONE calculating     Enc. perf. in the %s %d / 2 \n', int_par, iN)
            
    
            % ------ Start Calculating Null Encoding performances ------
            nullCuetBeta  = NaN(size(FIRs0, [1 2]));
            nullStimtBeta = NaN(size(FIRs0, [1 2]));
            n_perm        = 1; % This will take a while.
            rng('default');
            for ii = 1:n_perm
                nullidx = randperm(n_cond);
                for iT = 1:n_time
                    y_vox    = squeeze(FIRs0(:, iT, :))';
                    NullX    = [CueStimX(:, 1) CueStimX(nullidx, 2:3)];
                    nulltBeta    = pinv(NullX) * y_vox;
                    nullCuetBeta(:, iT)  = nulltBeta(2, :);
                    nullStimtBeta(:, iT) = nulltBeta(3, :);
                end
                [nullCuePC, ~]  = pca(nullCuetBeta', 'NumComponents', pccrit);
                [nullStimPC, ~] = pca(nullStimtBeta', 'NumComponents', pccrit);
    
                for iT = 1:n_time
                    y_vox    = squeeze(FIRs0(:, iT, :))';
                    nullCue_prj  = y_vox * nullCuePC;
                    nullStim_prj = y_vox * nullStimPC;
    
                    [~, cue1dim] = calc_normprjs(nullCue_prj);
                    [~, stim1dim] = calc_normprjs(nullStim_prj);
    
                    CueCue   = fitlm(mdl_cue, cue1dim); 
                    StimCue  = fitlm(mdl_stim, cue1dim);
                    StimStim = fitlm(mdl_stim, stim1dim);
                    CueStim  = fitlm(mdl_cue, stim1dim);
    
                    nullRSQs.(int_ntw).(int_par).XY_CueCue(iT, ii) = CueCue.Rsquared.Ordinary;
                    nullRSQs.(int_ntw).(int_par).XY_StimCue(iT, ii) = StimCue.Rsquared.Ordinary;
                    nullRSQs.(int_ntw).(int_par).XY_CueStim(iT, ii) = CueStim.Rsquared.Ordinary;
                    nullRSQs.(int_ntw).(int_par).XY_StimStim(iT, ii) = StimStim.Rsquared.Ordinary;
                end
            end
            % ------ End Calculating Null Encoding performances ------
            fprintf('DONE calculating Null Enc. perf. in the %s %d / 2 \n', int_par, iN)
            
            Cue2Cue      = RSQs.(int_ntw).(int_par).XY_CueCue;
            Cue2CueNull  = mean(nullRSQs.(int_ntw).(int_par).XY_CueCue, 2);
            Stim2Stim    = RSQs.(int_ntw).(int_par).XY_StimStim;
            Stim2StimNull= mean(nullRSQs.(int_ntw).(int_par).XY_StimStim, 2);
            cue_effects{iN, iiP} = Cue2Cue - Cue2CueNull;
            stim_effects{iN, iiP} = Stim2Stim - Stim2StimNull;
            
        end
    end
    save(fullfile(basedir, 'results', sprintf('DecodingInSubspaces_Result_%s.mat', pcstr)), ...
       'RSQs', 'nullRSQs', 'cue_effects', 'stim_effects')
end

%% 3. Draw Encoding Results
pccrit = 20; pcstr = sprintf('NumPC%d', pccrit);
load(fullfile(basedir, 'results', sprintf('DecodingInSubspaces_Result_%s.mat', pcstr)), ...
   'RSQs', 'nullRSQs');
figsize     = [120 135]; % orig [150 150]. Cereb Supple. [120 135]

savefig     = false;

cmaps       = cbrewer('qual', 'Set2', 3);
NTWRK       = fields(RSQs);
parcelnames = fields(RSQs.CerebrumStriatumCerebellum);  

kk = 1;
for iN = 1:numel(NTWRK)
    for iP = 1:numel(parcelnames)
        int_ntw = NTWRK{iN};
        int_par = parcelnames{iP};
        % figure('Name', sprintf('%s %s', int_ntw, int_par)); 
        
        Cue2Cue      = RSQs.(int_ntw).(int_par).XY_CueCue;
        Cue2Stim     = RSQs.(int_ntw).(int_par).XY_CueStim;
        Cue2CueNull  = mean(nullRSQs.(int_ntw).(int_par).XY_CueCue, 2);
        CueEff   = Cue2Cue - Cue2CueNull;
        [h1, p1, ~, stats1] = ttest(CueEff);
        
        Stim2Stim    = RSQs.(int_ntw).(int_par).XY_StimStim;
        Stim2Cue     = RSQs.(int_ntw).(int_par).XY_StimCue;
        Stim2StimNull= mean(nullRSQs.(int_ntw).(int_par).XY_StimStim, 2);
        StimEff  = Stim2Stim - Stim2StimNull;
        [h2, p2, ~, stats2] = ttest(StimEff);
        [h3, p3, ~, stats3] = ttest(Cue2Cue - Stim2Stim);

        pt_stats.CueEffs.(int_ntw).(int_par)  = [p1 stats1.tstat];
        pt_stats.StimEffs.(int_ntw).(int_par) = [p2 stats2.tstat];
        pt_stats.CueStim.(int_ntw).(int_par)  = [p3 stats3.tstat];

        plot_specificity_box_ycgosu(Cue2Cue, Stim2Stim, 'colors', ...
            cmaps([1 2], :));
        box off;set(gca, 'TickLength', [.02 .02], 'FontSize', 9, 'FontName', 'Helvetica'); 
        ylim([-0.3 1.4]);yticks([0 0.5 1])
        set(gcf, 'Position', [150*(iP-1) 300 figsize]);

        if savefig
            exportgraphics(gcf, fullfile(basedir, 'figures', 'step1_CalcSubspacesEncodingperf', ...
                sprintf('%s_%s_numPC%d_CueStimEncPerf.eps', int_ntw, int_par, pccrit)), 'ContentType','vector');
        end
        
        plot_specificity_box_ycgosu(Cue2Cue, Cue2CueNull, 'colors', ...
            cmaps([1 3], :));
        box off;set(gca, 'TickLength', [.02 .02], 'FontSize', 9, 'FontName', 'Helvetica'); 
        ylim([-0.3 1.4]);yticks([0 0.5 1])
        set(gcf, 'Position', [150*(iP-1) 300 figsize]);

        if savefig
            exportgraphics(gcf, fullfile(basedir, 'figures', 'step1_CalcSubspacesEncodingperf', ...
                sprintf('%s_%s_numPC%d_CueEncPerf.eps', int_ntw, int_par, pccrit)), 'ContentType','vector');
        end

        plot_specificity_box_ycgosu(Stim2Stim, Stim2StimNull, 'colors', ...
            cmaps([2 3], :));hold(gca, 'on')
        box off;set(gca, 'TickLength', [.02 .02], 'FontSize', 9, 'FontName', 'Helvetica'); 
        ylim([-0.3 1.4]);yticks([0 0.5 1])
        set(gcf, 'Position', [150*(iP-1) 150 figsize])
        
        if savefig
            exportgraphics(gcf, fullfile(basedir, 'figures', 'step1_CalcSubspacesEncodingperf', ...
                sprintf('%s_%s_numPC%d_StimEncPerf.eps', int_ntw, int_par, pccrit)), 'ContentType','vector');
        end
    end
end

%% 4. Draw Domain Effects
pccrit = 20; pcstr = sprintf('NumPC%d', pccrit);
load(fullfile(basedir, 'results', sprintf('DecodingInSubspaces_Result_%s.mat', pcstr)), ...
   'cue_effects', 'stim_effects')
cmaps = [102 195 165;
    252 141 98]./255;

savefig     = false;
% int_ntw = 'Cerebrum';

for iN = 1:numel(NTWRK)
    int_ntw = NTWRK{iN};
    for iP = 1:numel(parcelnames)
        int_par = parcelnames{iP};
        Effects = [cue_effects(iN, iP), stim_effects(iN, iP)];
        [~, ps_CueStim(iN, iP, :)] = cellfun(@ttest, Effects);
        violinplot(Effects, 'medc', 'k', 'mc', 'none', 'nopoints', ...
            'facecolor', cmaps, 'edgecolor', 'none', 'bw', 0.08, 'facealpha', 1, ...
            'x', [1 1.8]);box off;legend off;
        yticks([0 0.5 1]); ylim([-.2 1]);
        set(gcf, 'Position', [100 100 95 150])
        % set(gcf, 'Position', [100 100 110 120])
        yline(0, 'LineStyle','--', 'Color',[.3 .3 .3], 'LineWidth', 1.3);
        set(gca, 'LineWidth', 1.1, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', ...
            'TickLength', [.02 .02]);
        if savefig
            exportgraphics(gcf, fullfile(figoutdir, sprintf('effects_%s_%s_numPC%d.eps', ...
                int_ntw, int_par, pccrit)), 'ContentType','vector');
        end
        close all;
    end
end

%% 4-1. Draw cue and stim effects seperately.
pccrit = 20; pcstr = sprintf('NumPC%d', pccrit);
load(fullfile(basedir, 'results', sprintf('DecodingInSubspaces_Result_%s.mat', pcstr)), ...
   'cue_effects', 'stim_effects')
cmaps = [102 195 165;
    252 141 98]./255;

savefig     = false;
% int_ntw = 'Cerebrum';
xtickvals  = 1:2:14;
intervals  = 0.3;
xticks_are = [];
for ixt = 1:numel(xtickvals)
    xticks_are = [xticks_are, ...
        xtickvals(ixt) - intervals xtickvals(ixt) + intervals];
end

for iN = 1:numel(NTWRK)
    int_ntw = NTWRK{iN};
    cuestim_effect   = [cue_effects(iN, :);stim_effects(iN, :)];
    r_cuestim_effect = reshape(cuestim_effect, [], 1);
    
    violinplot(r_cuestim_effect', 'medc', 'k', 'mc', 'none', 'nopoints', ...
            'facecolor', repmat(cmaps, 7, 1), 'edgecolor', 'none', ...
            'bw', 0.08, 'facealpha', 1, ...
            'x', xticks_are)
            box off;legend off;
    xticks(xtickvals); yticks([0 0.5 1]); ylim([-.2 1]);
    set(gcf, 'Position', [100 100 300 120])
    yline(0, 'LineStyle','--', 'Color',[.3 .3 .3], 'LineWidth', 1.3);
    set(gca, 'LineWidth', 1.1, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', ...
        'TickLength', [.02 .02]);
    if savefig
        exportgraphics(gcf, fullfile(figoutdir, sprintf('effectsCueStim7_%s_%s_numPC%d.eps', ...
            int_ntw, int_par, pccrit)), 'ContentType','vector');
    end
    close all;

    transIdx = [5 6 7];
    uniIdx   = [1 2];

    transCueEffs  = mean(cat(2, cue_effects{iN, transIdx}), 2);
    transStimEffs = mean(cat(2, stim_effects{iN, transIdx}), 2);

    uniCueEffs    = mean(cat(2, cue_effects{iN, uniIdx}), 2);
    uniStimEffs   = mean(cat(2, stim_effects{iN, uniIdx}), 2);

    violinplot([uniCueEffs  transCueEffs], 'medc', 'k', 'mc', 'none', 'nopoints', ...
            'facecolor', repmat(cmaps(1, :), 2, 1), 'edgecolor', 'none', ...
            'bw', 0.08, 'facealpha', 1)
            box off;legend off;
    set(gcf, 'Position', [100 100 100 120])
    yline(0, 'LineStyle','--', 'Color',[.3 .3 .3], 'LineWidth', 1.3);
    set(gca, 'LineWidth', 1.1, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', ...
        'TickLength', [.02 .02]);
    if savefig
        exportgraphics(gcf, fullfile(figoutdir, sprintf('effectsCueUniTrans_%s_%s_numPC%d.eps', ...
            int_ntw, int_par, pccrit)), 'ContentType','vector');
    end
    close all;

    violinplot([uniStimEffs  transStimEffs], 'medc', 'k', 'mc', 'none', 'nopoints', ...
            'facecolor', repmat(cmaps(2, :), 2, 1), 'edgecolor', 'none', ...
            'bw', 0.08, 'facealpha', 1)
            box off;legend off;
    set(gcf, 'Position', [100 100 100 120])
    yline(0, 'LineStyle','--', 'Color',[.3 .3 .3], 'LineWidth', 1.3);
    set(gca, 'LineWidth', 1.1, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', ...
        'TickLength', [.02 .02]);
    if savefig
        exportgraphics(gcf, fullfile(figoutdir, sprintf('effectsStimUniTrans_%s_%s_numPC%d.eps', ...
            int_ntw, int_par, pccrit)), 'ContentType','vector');
    end
    close all;


end

%% SUBFUCTIONS
function [normvec, prj1dim] = calc_normprjs(inputPrj_t)
[maxidx, ~]   = find(squareform(pdist(inputPrj_t)) == max(pdist(inputPrj_t)));
cond_diffs    = (inputPrj_t(maxidx(1), :) - inputPrj_t(maxidx(2), :));
normvec       = cond_diffs ./ norm(cond_diffs);
prj1dim       = inputPrj_t*normvec';
end
