% All Dependencies (in README.md) should be installed and added in the path
clear;clc % 
% basedir   = '/path/to/github/PEPSI';
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


%% Calculating behavioral fits
mean_behvout = mean(behvout, 2);
pccrit = 20;
for iN = 1:numel(NTWRK)
    for iiP = 1:numel(parcelnames)
        int_ntw = NTWRK{iN};
        int_par = parcelnames{iiP};
        int_idx = parcelIndx.(int_ntw).(int_par);
        cue_tBetas  = neurAvg.CuetBeta(int_idx, :)';
        stim_tBetas = neurAvg.StimtBeta(int_idx, :)';
        [cuePC, ~, ~, ~, CueExpls]  = pca(cue_tBetas);
        [stimPC, ~, ~, ~, StimExpls] = pca(stim_tBetas);
        
        if isnumeric(pccrit)
            cuePC = cuePC(:, 1:pccrit);
            stimPC = stimPC(:, 1:pccrit);
        elseif ischar(pccrit)
            nPC = min(find(cumsum(CueExpls) > str2num(pccrit)));
            cuePC = cuePC(:, 1:nPC);
            nPC = min(find(cumsum(StimExpls) > str2num(pccrit)));
            stimPC = stimPC(:, 1:nPC);
        end

        FIRs        = neurAvg.FIR(int_idx, :, :);
        [n_vox, n_time, n_cond] = size(FIRs);
        timegroups = divide_ycgosu(n_time, n_time);

        r_FIRs = reshape(FIRs, n_vox, []);
        r_cuePrjs  = ((r_FIRs - mean(r_FIRs, 2))' * cuePC)';
        r_stimPrjs = ((r_FIRs - mean(r_FIRs, 2))' * stimPC)';
        cuePrj  = reshape(r_cuePrjs, [], n_time, n_cond);
        stimPrj = reshape(r_stimPrjs, [], n_time, n_cond);

        for iT = 1:numel(timegroups)
            cuePrj_t  = squeeze(mean(cuePrj(:, timegroups{iT}, :), 2))';
            stimPrj_t = squeeze(mean(stimPrj(:, timegroups{iT}, :), 2))';
            
    %       % Vector based on maximal distance.
            [R1, ~] = find(squareform(pdist(cuePrj_t)) == max(pdist(cuePrj_t)));
            [R2, ~] = find(squareform(pdist(stimPrj_t)) == max(pdist(stimPrj_t)));
            cue_diffs    = (cuePrj_t(R1(1), :) - cuePrj_t(R1(2), :));
            stim_diffs   = (stimPrj_t(R2(1), :) - stimPrj_t(R2(2), :));
            cue_normvec  = cue_diffs ./ norm(cue_diffs);
            stim_normvec = stim_diffs ./ norm(stim_diffs);
            cue1dim      = cuePrj_t*cue_normvec';
            stim1dim     = stimPrj_t*stim_normvec';

            Cue_HLN = [];
            for iC = 1:numel(Cues)
                Cue_HLN(iC, :) = mean(cue1dim(contains(condName, Cues{iC}), :));
            end
            if (Cue_HLN(1) - Cue_HLN(2) < 0) 
                Cue_HLN = -Cue_HLN;
                cue1dim = -cue1dim;
            end
            
            Stim_LVs = [];
            for iS = 1:numel(Stims)
                Stim_LVs(iS, :) = mean(stim1dim(contains(condName, Stims{iS}), :));
            end
            
            if (Stim_LVs(end) - Stim_LVs(1) < 0) 
                Stim_LVs = -Stim_LVs;
                stim1dim = -stim1dim;
            end
            
            pred_report = NaN(numel(condName), 1);
            pred_report(1:4) = Cue_HLN(1) + stim1dim(1:4);
            pred_report(5:8) = Cue_HLN(2) + stim1dim(5:8);
            pred_report(9:end) = Cue_HLN(3) + stim1dim(9:end);

            y = mean_behvout; % centering to No Cue LV3.
            y0 = [y - y(end-1)] ./ std(y);
            yfit = pred_report;
            yfit0 = [yfit - yfit(end-1)] ./ std(yfit);
            
            predout_bhvs.(int_ntw).(int_par)(iT, :) = 1 - (sum(abs((y0-yfit0)))./sum(abs((y0-mean(y0)))));
            predout_bhvs_rsq.(int_ntw).(int_par)(iT, :) = 1 - (sum((y0-yfit0).^2)./sum((y0-mean(y0)).^2));
            mseout_bhvs.(int_ntw).(int_par)(iT, :) = mean(abs(y0-yfit0));
            yfits.(int_ntw).(int_par)(iT, :) = yfit0;
        end
    end
end

%% Draw Real Behv & Recon Behv.
close all
idces = {1:4, 5:8, 9:11};
idloc = {2:5, 1:4, 2:4};
cmaps = [ 0.9882    0.6824    0.5686
        0.9843    0.4157    0.2902
        0.8706    0.1765    0.1490
        0.6471    0.0588    0.0824
        0.7412    0.8431    0.9059
        0.4196    0.6824    0.8392
        0.1922    0.5098    0.7412
        0.0314    0.3176    0.6118
        0.8000    0.8000    0.8000
        0.5882    0.5882    0.5882
        0.3882    0.3882    0.3882];

for iN = 1:numel(NTWRK)
    for iP = 1:numel(parcelnames)
        int_ntw = NTWRK{iN};
        int_par = parcelnames{iP}; 
        filestr = sprintf('%s_%s_numPC%d_BhvInEveryTime', ...
            int_ntw, int_par, pccrit);

        kk = 1;
        figure('Name', int_par)
        for i = 1:numel(idces)
            plot(idloc{i}, y0(idces{i}), 'LineWidth', 1, 'Color', [.1 .1 .1], ...
                'LineStyle','-');hold on;
            xlim([1 5]); ylim([-1.7 2.2])
            line([1 5], [0 0], 'Color', [.8 .8 .8], 'LineStyle', '--', 'LineWidth', 0.3);
            line([3 3], [2.5 -2], 'Color', [.8 .8 .8], 'LineStyle', '--', 'LineWidth', 0.3);
            meanvals = mean(yfits.(int_ntw).(int_par));
            stevals = ste(yfits.(int_ntw).(int_par));
            for iC = 1:11
                ww = cat(2, idloc{:});
                scatter(ww(iC), y0(iC), 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
            end
            errorbar(idloc{i}, meanvals(idces{i}), stevals(idces{i}), 'Color', [mean(cmaps(idces{i}, :)), 0.1], ...
                'LineWidth', 1.5, 'CapSize', 0.5)

            for ii = 1:numel(idloc{i})
                scatter(idloc{i}(ii), meanvals(idces{i}(ii)), 9, 'MarkerFaceColor', 'white', ...
                    'MarkerEdgeColor', cmaps(kk, :));
                kk = kk + 1;
            end
            box off;xticks(1:5)
            set(gca, 'LineWidth', 1.2, 'FontSize', 9, 'TickDir', 'out');
            set(gcf, 'Position', [100 100 170 160])
        end
        saveas(gcf, fullfile(basedir, 'figures', ...
            'step3_Traj2Behv', filestr), 'epsc');
        close all;
    end
end

%% Across 7 networks
close all;
% int_ntw = 'CerebrumStriatumCerebellum';
int_ntw = 'Cerebrum';

% Based on L1 error.
Drawmat1 = struct2array(predout_bhvs.(int_ntw), 2);
[~, srtidx] = sort(mean(Drawmat1), 'descend');
parcelnames(logical(ttest(Drawmat1)))
pDrawmat = Drawmat1(:, srtidx);
ylims = [-3 2]; ytick_lab = [-1:1];savestr = 'abs';

% Based on L2 error.
Drawmat2 = struct2array(predout_bhvs_rsq.(int_ntw), 2);
[~, srtidx] = sort(median(Drawmat2), 'descend');
for iP = 1:size(Drawmat2, 2)
    [l2_p(iP), l2_h(iP), l2_stats(iP)] = signrank(Drawmat2(:, iP), 0, 'tail', 'right');
    if l2_h(iP) ~= 0
        parcelnames(iP)
    end
end
pDrawmat = Drawmat2(:, srtidx);
ylims = [-4 1.6]; ytick_lab = [-4:1];savestr = 'rsq';


%
close all;
violinplot(pDrawmat, 'medc', 'black', 'mc', [0 249 0]./255, 'nopoints', ...
        'facecolor', cmaps7(srtidx, :), 'edgecolor', [.8 .8 .8], 'bw', 0.2, 'facealpha', 1);box off;
yline(0, 'LineStyle','--','LineWidth',1.2, 'Color',[.5 .5 .5])
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02], 'LineWidth', 1.5, 'FontSize', 9);
set(gcf, 'Position', [100 100 210 150]);legend off;
ylim(ylims);
yticks(ytick_lab)
saveas(gcf, fullfile(basedir, 'figures', 'step3_Traj2Behv',...
    sprintf('%s_numPC%d_ReconFits_%s', int_ntw, pccrit, savestr)), 'epsc');
% close all;
