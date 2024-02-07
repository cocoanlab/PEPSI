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

%%
close all;
% mkdir(outdir)
for iN = 2
    for iP = 1:numel(parcelnames)
        int_ntw  = NTWRK{iN};
        int_par  = parcelnames{iP};
        int_idx  = parcelIndx.(int_ntw).(int_par);
        outname  = sprintf('3dTraj_%s_%s.png', int_ntw, int_par);
        FIRs     = neurAvg.FIR(int_idx, :, :);
        CuetBeta = neurAvg.CuetBeta(int_idx, :);
        StimtBeta= neurAvg.StimtBeta(int_idx, :);
        CuePC  = pca(CuetBeta');
        StimPC = pca(StimtBeta');
        FIRs0    = FIRs - mean(FIRs, [2 3]); 
        [n_vox, n_time, n_cond] = size(FIRs0);

        figure('Name', int_par);subplot 211;
        for iC = 1:n_cond
            CueDraw  = smoothdata(squeeze(FIRs0(:, :, iC))' * CuePC, 'gaussian', 9);
            plot3(CueDraw(:, 1), CueDraw(:, 2), CueDraw(:, 3), 'LineWidth', 1.4, 'Color', cmaps(iC, :));hold on;
            scatter3(CueDraw(end, 1), CueDraw(end, 2), CueDraw(end, 3), 30, [.2 .2 .2], 'filled');
            set(gca, 'LineWidth', 1.5, 'fontsize', 9, 'GridAlpha', 0.05, ...
                'TickDir' ,'out', 'TickLength', [0.04 0.04]); 
            grid on;
        end
        % xlabel('PC1');ylabel('PC2');zlabel('PC3');
        
        subplot 212;
        for iC = 1:n_cond
            StimDraw = smoothdata(squeeze(FIRs0(:, :, iC))' * StimPC, 'gaussian', 9);
            plot3(StimDraw(:, 1), StimDraw(:, 2), StimDraw(:, 3), 'LineWidth', 1.4, 'Color', cmaps(iC, :));hold on;
            scatter3(StimDraw(end, 1), StimDraw(end, 2), StimDraw(end, 3), 30, [.2 .2 .2], 'filled');
            set(gca, 'LineWidth', 1.5, 'fontsize', 9, 'GridAlpha', 0.05, 'TickDir' ,'out'); 
            grid on;
        end
        % xlabel('PC1');ylabel('PC2');zlabel('PC3');
        set(gcf, 'Position', [200*(iP-1) 300 200 300])
    end
end
