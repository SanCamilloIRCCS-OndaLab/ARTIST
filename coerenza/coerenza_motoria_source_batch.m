function coerenza_motoria_source_batch(outputMatlabRoot, outRoot, overwriteOutput, sourceRoot)
%COERENZA_MOTORIA_SOURCE_BATCH
% Coerenza cortico-sinergica con ROI motorie EEG in source-space (NET output: sources_eeg.mat).
% Non modifica NET: legge solo file gia' prodotti.
% Se il source non e' disponibile per un pair, fa fallback scalp-proxy e lo segnala nel summary.
%
% English overview:
% - consumes aligned MFF outputs produced by allineamento_eeg_emg_sync;
% - builds trials from reconstructed logical events (beep->start window);
% - computes cortico-synergy coherence with source ROI when available,
%   otherwise uses a scalp-proxy fallback.
%
% Input default:
%   outputMatlabRoot = ../batch_otb_mff/output_matlab
%   outRoot          = ./   (cartella coerenza)
%
% Requisiti:
%   - report.json prodotto da allineamento_eeg_emg_sync.m
%   - plugin EEGLAB MFFMatlabIO disponibile (pop_mffimport)

if nargin < 1 || isempty(outputMatlabRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    outputMatlabRoot = fullfile(fileparts(scriptDir), 'batch_otb_mff', 'output_matlab_sync');
end
if nargin < 2 || isempty(outRoot)
    outRoot = fullfile(fileparts(mfilename('fullpath')), 'coerenza_motoria_source_sync');
end
if nargin < 3 || isempty(overwriteOutput)
    overwriteOutput = true;
end
if nargin < 4
    sourceRoot = '';
end

if ~isfolder(outputMatlabRoot)
    error('Output matlab root non trovato: %s', outputMatlabRoot);
end
if ~isfolder(outRoot)
    mkdir(outRoot);
end

initEeglab();
roiCfg = defaultMotorRoiConfig();
sourceDb = indexSourceFiles(sourceRoot, outputMatlabRoot);

reports = dir(fullfile(outputMatlabRoot, '**', 'report.json'));
if isempty(reports)
    error('Nessun report.json trovato in %s', outputMatlabRoot);
end

summaryHeader = {'pair','phase','status','message','n_trials','n_eeg_ch','n_emg_ch', ...
    'n_synergies','vaf_global','vaf_min_channel','c95', ...
    'beta_mean_all','beta_peak_all','beta_peak_hz_all','beta_auc_all', ...
    'zcoh_beta_mean_all','zcoh_beta_peak_all','zcoh_beta_peak_hz_all','zcoh_threshold_beta_all','zcoh_sigfrac_beta_all', ...
    'zcoh_surrogate_p_min','zcoh_surrogate_q_min','zcoh_surrogate_sig_count', ...
    'beta_mean_lagopt_all','lag_opt_ms_all','psi_beta_all', ...
    'im_beta_mean_all','beta_ci95_low_all','beta_ci95_high_all','beta_stability_cv_all', ...
    'gesture_delta_late_minus_early_all','lateralization_li_all', ...
    'm1_l_beta_all','m1_r_beta_all','m1_delta_l_minus_r_all', ...
    'm1_delta_l_minus_r_early','m1_delta_l_minus_r_middle','m1_delta_l_minus_r_late', ...
    'm1_delta_l_minus_r_13_20','m1_delta_l_minus_r_20_30', ...
    'surrogate_p_min','surrogate_q_min','surrogate_sig_count', ...
    'multivar_r1','multivar_r1_perm_p','covar_corr_emgdrive_beta','covar_corr_duration_beta', ...
    'nmf_stability', ...
    'eeg_mode','source_file','source_status', ...
    'delta_beta_mean_source_minus_scalp','delta_beta_peak_source_minus_scalp','delta_beta_auc_source_minus_scalp', ...
    'roi_used','roi_channels_used', ...
    'out_dir'};
summaryRows = {};

for i = 1:numel(reports)
    reportPath = fullfile(reports(i).folder, reports(i).name);
    try
        R = jsondecode(fileread(reportPath));
        pairName = getFieldDef(R, 'pair', sprintf('pair_%03d', i));
        phase = getFieldDef(R, 'phase', 'session');
        outDir = fullfile(outRoot, pairName, phase);

        if isfolder(outDir) && overwriteOutput
            rmdir(outDir, 's');
        end
        if ~isfolder(outDir)
            mkdir(outDir);
        end

        if ~isfield(R,'output_mff') || ~isfolder(R.output_mff)
            error('output_mff mancante o non valido in report');
        end
        if ~isfield(R,'source_otb') || ~isfile(R.source_otb)
            error('source_otb mancante o non valido in report');
        end

        fprintf('\n[COERENZA] %s/%s\n', pairName, phase);

        % Load aligned EEG and derive trial windows from logical beep->start events.
        [eegData, eegFs, eegLabels] = loadMffEeg(R.output_mff);
        trials = buildTrialsFromMffEvents(R.output_mff, eegFs, size(eegData,2));
        if isempty(trials)
            error('Nessun trial beep/start valido trovato in MFF');
        end
        trialsSec = [(trials(:,1)-1)/eegFs, (trials(:,2)-1)/eegFs];

        [emgRaw, emgFs, emgLabels] = readOtbBipolar10(R.source_otb);
        trA = 1;
        trB = 0;
        if isfield(R,'transform') && isstruct(R.transform)
            if isfield(R.transform,'a'), trA = double(R.transform.a); end
            if isfield(R.transform,'b'), trB = double(R.transform.b); end
        end
        eegMode = 'scalp_proxy';
        sourceFile = '';
        sourceStatus = 'not_found';
        deltaBetaMean = NaN;
        deltaBetaPeak = NaN;
        deltaBetaAuc  = NaN;

        % Prefer source-space ROI extraction; fallback to scalp-proxy if source is unavailable.
        [sourceFile, sourceStatus] = resolveSourceFile(R, reportPath, pairName, phase, sourceDb);
        if ~isempty(sourceFile)
            try
                [eegSourceCont, eegFsUse, roiMapSource] = loadMotorRoiFromSource(sourceFile, roiCfg);
                emgAligned = preprocessAndAlignEmg(emgRaw, emgFs, trA, trB, eegFsUse, size(eegSourceCont,2));
                trialsUse = secondsToTrials(trialsSec, eegFsUse, size(eegSourceCont,2));
                [eegTrials, emgTrials, trialInfo] = sliceTrials(eegSourceCont, emgAligned, trialsUse);
                if isempty(eegTrials)
                    error('Nessun trial valido dopo slicing in source-space');
                end
                eegMode = 'source_net';
                eegSelIdx = 1:size(eegTrials{1},1);
                eegSelLabels = {roiMapSource.name};
                roiMap = roiMapSource;
            catch MEsrc
                sourceStatus = ['failed: ' MEsrc.message];
            end
        end

        if ~strcmp(eegMode, 'source_net')
            emgAligned = preprocessAndAlignEmg(emgRaw, emgFs, trA, trB, eegFs, size(eegData,2));
            [eegTrialsRaw, emgTrials, trialInfo] = sliceTrials(eegData, emgAligned, trials);
            if isempty(eegTrialsRaw)
                error('Dopo slicing non ci sono trial validi');
            end
            [eegTrials, eegSelIdx, eegSelLabels, roiMap] = selectMotorRoiChannels(eegTrialsRaw, eegLabels, roiCfg);
            eegFsUse = eegFs;
        end

        [W, Htrials, nSyn, vafGlobal, vafMinChannel, nmfStability] = extractSynergiesNmf(emgTrials, 6);
        [cohMean, f, metrics] = computeCorticoSynergyCoherence(eegTrials, Htrials, eegFsUse, eegSelLabels, emgTrials);
        metrics = addZcohAndSynergyClusters(metrics, cohMean, f);

        % Confronto source vs scalp (stesso fs/segmentazione del run source)
        if strcmp(eegMode, 'source_net')
            [scalpCont, roiMapScalp] = buildScalpRoiContinuous(eegData, eegFs, eegLabels, roiCfg, eegFsUse, size(emgAligned,2)); %#ok<ASGLU>
            nCmp = min(size(scalpCont,2), size(emgAligned,2));
            trialsUse = secondsToTrials(trialsSec, eegFsUse, nCmp);
            scalpCont = scalpCont(:,1:nCmp);
            emgCmp = emgAligned(:,1:nCmp);
            [eegTrialsScalp, ~, ~] = sliceTrials(scalpCont, emgCmp, trialsUse);
            if ~isempty(eegTrialsScalp)
                scalpLabels = {roiMapScalp.name};
                [~, ~, metricsScalp] = computeCorticoSynergyCoherence(eegTrialsScalp, Htrials, eegFsUse, scalpLabels, [], struct('basicOnly', true));
                deltaBetaMean = metrics.betaMeanAll - metricsScalp.betaMeanAll;
                deltaBetaPeak = metrics.betaPeakAll - metricsScalp.betaPeakAll;
                deltaBetaAuc  = metrics.betaAucAll  - metricsScalp.betaAucAll;
            end
        end

        save(fullfile(outDir, 'risultati_coerenza.mat'), ...
            'W','Htrials','nSyn','vafGlobal','vafMinChannel', ...
            'cohMean','f','metrics','trialInfo','eegSelIdx','eegSelLabels','roiMap','roiCfg','emgLabels', ...
            'nmfStability', ...
            'eegMode','sourceFile','sourceStatus','deltaBetaMean','deltaBetaPeak','deltaBetaAuc', ...
            'trA','trB','reportPath','-v7.3');

        % Skip figure generation to reduce runtime and avoid PNG outputs.
        writeSynergyMetricsCsv(fullfile(outDir, 'metriche_sinergie.csv'), metrics);
        writeZcoherenceMetricsCsv(fullfile(outDir, 'metriche_zcoherence.csv'), metrics);
        writeSynergyClusterCsv(fullfile(outDir, 'cluster_sinergie.csv'), metrics);
        writeRoiMetricsCsv(fullfile(outDir, 'metriche_roi.csv'), metrics);
        writeM1LateralizationCsv(fullfile(outDir, 'metriche_m1_dx_sx.csv'), metrics);
        writeRoiStateDependentCsv(fullfile(outDir, 'metriche_csc_state_dependent.csv'), metrics);
        writeTrialInfoCsv(fullfile(outDir, 'trial_info.csv'), trialInfo, eegFsUse);

        fprintf('[OK] %s/%s | mode=%s | trials=%d | chEEG=%d | chEMG=%d | Nsyn=%d | VAF=%.3f\n', ...
            pairName, phase, eegMode, numel(eegTrials), size(eegTrials{1},1), size(emgTrials{1},1), nSyn, vafGlobal);

        summaryRows(end+1,:) = {pairName,phase,'OK','', ...
            num2str(numel(eegTrials)), num2str(size(eegTrials{1},1)), num2str(size(emgTrials{1},1)), ...
            num2str(nSyn), sprintf('%.6f',vafGlobal), sprintf('%.6f',vafMinChannel), ...
            sprintf('%.6f',metrics.c95), ...
            sprintf('%.6f',metrics.betaMeanAll), sprintf('%.6f',metrics.betaPeakAll), ...
            sprintf('%.6f',metrics.betaPeakHzAll), sprintf('%.6f',metrics.betaAucAll), ...
            fmt6(metrics.zCohBetaMeanAll), fmt6(metrics.zCohBetaPeakAll), fmt6(metrics.zCohBetaPeakHzAll), ...
            fmt6(metrics.zCohThresholdBetaAll), fmt6(metrics.zCohSigFracBetaAll), ...
            fmt6(metrics.zCohSurrogatePMin), fmt6(metrics.zCohSurrogateQMin), num2str(metrics.zCohSurrogateSigCount), ...
            sprintf('%.6f',metrics.betaMeanLagOptAll), sprintf('%.6f',metrics.lagOptMsAll), sprintf('%.6f',metrics.psiBetaAll), ...
            sprintf('%.6f',metrics.imBetaMeanAll), sprintf('%.6f',metrics.betaCi95LowAll), ...
            sprintf('%.6f',metrics.betaCi95HighAll), sprintf('%.6f',metrics.betaStabilityCvAll), ...
            sprintf('%.6f',metrics.gestureDeltaLateEarlyAll), sprintf('%.6f',metrics.lateralizationLiAll), ...
            fmt6(metrics.m1LeftBetaAll), fmt6(metrics.m1RightBetaAll), fmt6(metrics.m1DeltaLminusRAll), ...
            fmt6(metrics.m1DeltaLminusREarlyAll), fmt6(metrics.m1DeltaLminusRMiddleAll), fmt6(metrics.m1DeltaLminusRLateAll), ...
            fmt6(metrics.m1DeltaLminusRLowBetaAll), fmt6(metrics.m1DeltaLminusRHighBetaAll), ...
            sprintf('%.6f',metrics.surrogatePMin), sprintf('%.6f',metrics.surrogateQMin), ...
            num2str(metrics.surrogateSigCount), sprintf('%.6f',metrics.multivarR1), ...
            sprintf('%.6f',metrics.multivarR1PermP), sprintf('%.6f',metrics.covarCorrEmgDriveBeta), ...
            sprintf('%.6f',metrics.covarCorrDurationBeta), sprintf('%.6f',nmfStability), ...
            eegMode, sourceFile, sourceStatus, fmt6(deltaBetaMean), fmt6(deltaBetaPeak), fmt6(deltaBetaAuc), ...
            strjoin({roiMap.name}, '|'), summarizeRoiMap(roiMap), ...
            outDir}; %#ok<AGROW>

    catch ME
        [pairName, phase] = safePairPhaseFromReportPath(reportPath);
        fprintf('[ERR] %s/%s | %s\n', pairName, phase, ME.message);
        errRow = repmat({''}, 1, numel(summaryHeader));
        errRow{1} = pairName;
        errRow{2} = phase;
        errRow{3} = 'ERR';
        errRow{4} = ME.message;
        summaryRows(end+1,:) = errRow; %#ok<AGROW>
    end
end

summaryCsv = fullfile(outRoot, 'summary_coerenza_motoria_source.csv');
writeCsv(summaryCsv, summaryHeader, summaryRows);
writePrePostBetaSummary(fullfile(outRoot, 'summary_beta_pre_post_motoria_source.csv'), summaryRows, summaryHeader);
fprintf('\nSummary: %s\n', summaryCsv);
fprintf('Nota: nessuna modifica a NET. Modalita'' source usa solo sources_eeg.mat se disponibile.\n');
end

function initEeglab()
persistent done
if ~isempty(done) && done
    return;
end
eeglabRoot = '/Users/martinaregazzetti/Desktop/KU Leuven/Artist project/matlab /eeglab2025.1.0';
if ~isfolder(eeglabRoot)
    error('EEGLAB non trovato: %s', eeglabRoot);
end
addpath(eeglabRoot);
eeglab nogui;
if isempty(which('pop_mffimport'))
    error('pop_mffimport non disponibile. Verifica plugin MFFMatlabIO.');
end
done = true;
end

function [data, fs, labels] = loadMffEeg(mffPath)
EEG = pop_mffimport(mffPath);
data = double(EEG.data);
fs = double(EEG.srate);
labels = cell(1, EEG.nbchan);
for i = 1:EEG.nbchan
    if isfield(EEG,'chanlocs') && numel(EEG.chanlocs) >= i && isfield(EEG.chanlocs(i),'labels')
        labels{i} = char(string(EEG.chanlocs(i).labels));
    else
        labels{i} = sprintf('EEG_%03d', i);
    end
end
end

function sourceDb = indexSourceFiles(sourceRoot, outputMatlabRoot)
sourceDb = struct('paths',{{}}, 'names',{{}});
roots = {};
if nargin >= 1 && ~isempty(sourceRoot) && isfolder(sourceRoot)
    roots{end+1} = sourceRoot; %#ok<AGROW>
end
if nargin >= 2 && ~isempty(outputMatlabRoot) && isfolder(outputMatlabRoot)
    roots{end+1} = outputMatlabRoot; %#ok<AGROW>
    p = fileparts(outputMatlabRoot);
    if isfolder(p), roots{end+1} = p; end %#ok<AGROW>
end
if isempty(roots)
    return;
end
allp = {};
for i = 1:numel(roots)
    d = dir(fullfile(roots{i}, '**', 'sources_eeg.mat'));
    for k = 1:numel(d)
        allp{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end
allp = unique(allp, 'stable');
sourceDb.paths = allp;
sourceDb.names = lower(string(allp));
end

function [sourceFile, status] = resolveSourceFile(R, reportPath, pairName, phase, sourceDb)
sourceFile = '';
status = 'not_found';
side = inferSideFromPair(pairName);

% 1) esplicito nel report
candidateFields = {'source_eeg','source_source','source_sources_eeg'};
for i = 1:numel(candidateFields)
    f = candidateFields{i};
    if isfield(R,f) && ischar(R.(f)) && isfile(R.(f))
        sourceFile = R.(f);
        status = 'from_report';
        return;
    end
end

% 2) vicino al report
localCandidates = { ...
    fullfile(fileparts(reportPath), 'sources_eeg.mat'), ...
    fullfile(fileparts(reportPath), 'eeg_source', 'sources_eeg.mat')};
for i = 1:numel(localCandidates)
    if isfile(localCandidates{i})
        sourceFile = localCandidates{i};
        status = 'local_match';
        return;
    end
end

% 3) indice globale (match per token)
if ~isempty(sourceDb.paths)
    toks = tokenizePair(pairName, phase);
    if ~isempty(toks)
        scores = zeros(numel(sourceDb.paths),1);
        for i = 1:numel(sourceDb.paths)
            nm = lower(sourceDb.paths{i});
            if ~isSideCompatible(side, nm)
                continue;
            end
            for t = 1:numel(toks)
                if contains(nm, toks{t})
                    scores(i) = scores(i) + 1;
                end
            end
        end
        [mx, ix] = max(scores);
        if mx >= 2
            sourceFile = sourceDb.paths{ix};
            status = sprintf('indexed_match_%d', mx);
            return;
        end
    end
    status = 'indexed_not_matched';
end

function tf = isSideCompatible(side, pathLower)
if isempty(side)
    tf = true;
    return;
end
switch side
    case 'dx'
        tf = contains(pathLower, 'dx') || contains(pathLower, 'destra') || contains(pathLower, 'right');
    case 'sx'
        tf = contains(pathLower, 'sx') || contains(pathLower, 'sin') || contains(pathLower, 'left');
    otherwise
        tf = true;
end
end

function side = inferSideFromPair(pairName)
s = lower(char(string(pairName)));
if contains(s, '_dx_') || endsWith(s,'_dx') || contains(s,'destra') || contains(s,'right')
    side = 'dx';
elseif contains(s, '_sx_') || endsWith(s,'_sx') || contains(s,'sin') || contains(s,'left')
    side = 'sx';
else
    side = '';
end
end
end

function toks = tokenizePair(pairName, phase)
s = lower(char(string(pairName)));
s = regexprep(s, '[^a-z0-9]+', '_');
parts = regexp(s, '_+', 'split');
parts = parts(~cellfun(@isempty, parts));
stop = {'bids','ses','session','pair','pre','post'};
toks = {};
for i = 1:numel(parts)
    p = parts{i};
    if numel(p) < 2 || any(strcmp(p, stop))
        continue;
    end
    toks{end+1} = p; %#ok<AGROW>
end
if nargin >= 2
    ph = lower(char(string(phase)));
    if ~isempty(ph) && ~any(strcmp(toks, ph))
        toks{end+1} = ph; %#ok<AGROW>
    end
end
toks = unique(toks, 'stable');
end

function [eegCont, fs, roiMap] = loadMotorRoiFromSource(sourceFile, roiCfg)
S = load(sourceFile, 'source');
if ~isfield(S,'source')
    error('Campo source mancante in %s', sourceFile);
end
src = S.source;

req = {'sensor_data','imagingkernel','pca_projection','pos_mni','inside_mni','spatial_filter_mni','time'};
for i = 1:numel(req)
    if ~isfield(src, req{i})
        error('Campo source.%s mancante', req{i});
    end
end

sensor = double(src.sensor_data);
K = double(src.imagingkernel);
P = double(src.pca_projection);

if size(P,2) == size(K,1)
    Kp = P * K; % [Ninside_native x Nsensor]
elseif mod(size(K,1),3) == 0
    % fallback conservativo se pca_projection non compatibile
    Kp = (K(1:3:end,:) + K(2:3:end,:) + K(3:3:end,:)) / 3;
else
    error('Mappatura kernel non compatibile (pca_projection/imagingkernel).');
end

if size(Kp,2) ~= size(sensor,1)
    if size(Kp,2) == size(sensor,2)
        sensor = sensor';
    else
        error('Dimensioni non compatibili tra kernel e sensor_data.');
    end
end

tt = double(src.time(:)');
if numel(tt) < 2
    error('source.time non valido');
end
fs = 1 / median(diff(tt));

posMni = double(src.pos_mni);
insideMni = logical(src.inside_mni(:));
if size(posMni,1) ~= numel(insideMni)
    error('pos_mni/inside_mni non allineati');
end
posInsideMni = posMni(insideMni,:);
Wmni = double(src.spatial_filter_mni); % [Ninside_mni x Ninside_native]
if size(Wmni,1) ~= size(posInsideMni,1) || size(Wmni,2) ~= size(Kp,1)
    error('spatial_filter_mni non compatibile con pos_mni/inside_mni o kernel proiettato.');
end

nRoi = numel(roiCfg);
eegCont = zeros(nRoi, size(sensor,2));
roiMap = struct('name',{},'mni_xyz',{},'mni_nearest_xyz',{},'mni_nearest_dist_mm',{},'mni_n_voxels',{},'idx',{},'labels',{});

for r = 1:nRoi
    xyz = double(roiCfg(r).mni_xyz(:)');
    d2 = sum((posInsideMni - xyz).^2, 2);
    [sd, ord] = sort(d2, 'ascend');
    nvox = min(25, numel(ord));
    nn = ord(1:nvox);
    w = mean(Wmni(nn,:), 1);
    if ~any(abs(w) > 0)
        nn = ord(1);
        w = Wmni(nn,:);
    end
    filt = w * Kp;          % [1 x Nsensor]
    sig = filt * sensor;    % [1 x Nt]
    sig = sig - mean(sig);
    ssd = std(sig);
    if ssd > eps, sig = sig ./ ssd; end
    eegCont(r,:) = sig;

    roiMap(r).name = roiCfg(r).name; %#ok<AGROW>
    roiMap(r).mni_xyz = xyz; %#ok<AGROW>
    roiMap(r).mni_nearest_xyz = posInsideMni(ord(1),:); %#ok<AGROW>
    roiMap(r).mni_nearest_dist_mm = sqrt(sd(1)); %#ok<AGROW>
    roiMap(r).mni_n_voxels = nvox; %#ok<AGROW>
    roiMap(r).idx = NaN; %#ok<AGROW>
    roiMap(r).labels = {sprintf('source_mni_%s', roiCfg(r).name)}; %#ok<AGROW>
end
end

function trials = secondsToTrials(trialsSec, fs, nPts)
trials = zeros(0,2);
minTrialSamples = minReachTrialSamples(fs);
for i = 1:size(trialsSec,1)
    i1 = max(1, floor(trialsSec(i,1)*fs) + 1);
    i2 = min(nPts, ceil(trialsSec(i,2)*fs));
    if i2 <= i1
        continue;
    end
    % Keep only trials long enough to represent a real reaching segment.
    if (i2 - i1 + 1) < minTrialSamples
        continue;
    end
    trials(end+1,:) = [i1 i2]; %#ok<AGROW>
end
end

function [roiCont, roiMap] = buildScalpRoiContinuous(eegData, eegFs, eegLabels, roiCfg, fsTarget, nTarget)
[tmpTrials, ~, ~, roiMap] = selectMotorRoiChannels({eegData}, eegLabels, roiCfg);
roiCont = tmpTrials{1};
if abs(fsTarget - eegFs) < 1e-12 && size(roiCont,2) == nTarget
    return;
end
tOrig = (0:size(roiCont,2)-1) / eegFs;
tNew = (0:nTarget-1) / fsTarget;
Y = zeros(size(roiCont,1), numel(tNew));
for r = 1:size(roiCont,1)
    Y(r,:) = interp1(tOrig, roiCont(r,:), tNew, 'linear', 'extrap');
end
roiCont = Y;
end

function makeSourceVsScalpComparePlot(metricsSource, metricsScalp, outDir)
fig = figure('Visible','off','Color','w');
vals = [metricsSource.betaMeanAll metricsScalp.betaMeanAll; ...
        metricsSource.betaPeakAll metricsScalp.betaPeakAll; ...
        metricsSource.betaAucAll  metricsScalp.betaAucAll];
bar(vals);
set(gca, 'XTickLabel', {'beta mean','beta peak','beta AUC'});
ylabel('Valore');
legend({'source','scalp'}, 'Location', 'best');
title('Confronto Source vs Scalp (ROI motorie)');
grid on;
saveas(fig, fullfile(outDir, 'confronto_source_vs_scalp_beta.png'));
close(fig);
end

function trials = buildTrialsFromMffEvents(mffPath, fs, nPts)
[times, codes] = readMffEvents(mffPath);
% Trials are defined on the movement preparation window:
% beep -> first start before next beep.
beeps = sort(times(strcmp(codes,'beep')));
starts = sort(times(strcmp(codes,'start')));
trials = zeros(0,2);
minTrialSamples = minReachTrialSamples(fs);
if isempty(beeps) || isempty(starts)
    return;
end

for i = 1:numel(beeps)
    s = beeps(i);
    if i < numel(beeps)
        sNext = beeps(i+1);
    else
        sNext = inf;
    end
    candStart = starts(starts > s & starts < sNext);
    if isempty(candStart)
        continue;
    end
    e = candStart(1);
    i1 = max(1, floor(s*fs) + 1);
    i2 = min(nPts, ceil(e*fs));
    if i2 <= i1
        continue;
    end
    % Reject very short beep->start windows because they destabilize
    % coherence estimates and do not reflect the full reaching movement.
    if (i2 - i1 + 1) < minTrialSamples
        continue;
    end
    trials(end+1,:) = [i1 i2]; %#ok<AGROW>
end
end

function n = minReachTrialSamples(fs)
% Enforce a minimum reaching duration of 0.50 s in every sampling space.
n = max(50, round(0.50 * fs));
end

function [times, codes] = readMffEvents(mffPath)
evPath = fullfile(mffPath, 'Events_8 DINs.xml');
infoPath = fullfile(mffPath, 'info.xml');
if ~isfile(evPath)
    error('Events_8 DINs.xml non trovato in %s', mffPath);
end
if ~isfile(infoPath)
    error('info.xml non trovato in %s', mffPath);
end

infoTxt = fileread(infoPath);
recordIso = oneToken(infoTxt, '<recordTime>([^<]+)</recordTime>');
recordEpoch = isoToEpoch(recordIso);

txt = fileread(evPath);
allEv = regexp(txt,'<event>\s*<beginTime>([^<]+)</beginTime>[\s\S]*?<code>([^<]+)</code>[\s\S]*?</event>','tokens');
if isempty(allEv)
    error('Nessun evento in %s', evPath);
end

times = zeros(1, numel(allEv));
codes = cell(1, numel(allEv));
% Convert absolute event timestamps to seconds from recording start.
for i = 1:numel(allEv)
    times(i) = isoToEpoch(allEv{i}{1}) - recordEpoch;
    codes{i} = char(string(allEv{i}{2}));
end
end

function [emgRaw, fs, labels] = readOtbBipolar10(otbPath)
tmpDir = tempname;
mkdir(tmpDir);
c = onCleanup(@() cleanupTemp(tmpDir));

untar(otbPath, tmpDir);
xmlList = dir(fullfile(tmpDir, '*.xml'));
sigList = dir(fullfile(tmpDir, '*.sig'));
if isempty(xmlList) || isempty(sigList)
    error('OTB non valido: %s', otbPath);
end
sigPath = fullfile(sigList(1).folder, sigList(1).name);

DOM = [];
dev = [];
fs = NaN;
nChannels = NaN;
for ix = 1:numel(xmlList)
    xmlPath = fullfile(xmlList(ix).folder, xmlList(ix).name);
    try
        D = xmlread(xmlPath);
    catch
        continue;
    end
    r = D.getDocumentElement;
    fsTry = str2double(char(r.getAttribute('SampleFrequency')));
    chTry = str2double(char(r.getAttribute('DeviceTotalChannels')));
    if ~isnan(fsTry) && ~isnan(chTry) && chTry > 0
        DOM = D;
        dev = r;
        fs = fsTry;
        nChannels = chTry;
        break;
    end
end
if isempty(DOM)
    error('Metadata OTB non validi nel file XML');
end

adapters = dev.getElementsByTagName('Adapter');
selIdx = [];
labels = {};
for i = 0:(adapters.getLength()-1)
    ad = adapters.item(i);
    desc = lower(char(ad.getAttribute('Description')));
    idv = lower(char(ad.getAttribute('ID')));
    if ~(contains(desc,'bipolar') || contains(idv,'ad8x2'))
        continue;
    end
    st = str2double(char(ad.getAttribute('ChannelStartIndex')));
    chNodes = ad.getElementsByTagName('Channel');
    for j = 0:(chNodes.getLength()-1)
        ch = chNodes.item(j);
        idx = str2double(char(ch.getAttribute('Index')));
        mus = char(ch.getAttribute('Muscle'));
        if isnan(idx), continue; end
        selIdx(end+1) = st + idx + 1; %#ok<AGROW>
        if isempty(strtrim(mus))
            labels{end+1} = sprintf('BIP_%02d', numel(selIdx)); %#ok<AGROW>
        else
            labels{end+1} = sprintf('BIP_%02d_%s', numel(selIdx), regexprep(strtrim(mus), '\s+', '_')); %#ok<AGROW>
        end
    end
end

if isempty(selIdx)
    selIdx = 1:min(10, nChannels);
    labels = arrayfun(@(k) sprintf('BIP_%02d',k), 1:numel(selIdx), 'UniformOutput', false);
else
    selIdx = selIdx(1:min(10,numel(selIdx)));
    labels = labels(1:numel(selIdx));
end

fid = fopen(sigPath,'r');
raw = fread(fid, inf, 'int16=>double', 0, 'l');
fclose(fid);
N = floor(numel(raw)/nChannels);
raw = raw(1:N*nChannels);
data = reshape(raw, [nChannels,N]);
emgRaw = data(selIdx, :);
end

function emgAligned = preprocessAndAlignEmg(emgRaw, fsEmg, a, b, fsEeg, nEegPts)
tOtb = (0:size(emgRaw,2)-1) / fsEmg;
tEegFromOtb = a * tOtb + b;
tGrid = (0:nEegPts-1) / fsEeg;

emgAligned = zeros(size(emgRaw,1), numel(tGrid));
for ch = 1:size(emgRaw,1)
    s = double(emgRaw(ch,:));
    s = preprocessEmgSignal(s, fsEmg);
    v = interp1(tEegFromOtb, s, tGrid, 'linear', NaN);
    bad = isnan(v);
    if any(bad)
        good = ~bad;
        if any(good)
            v(bad) = interp1(tGrid(good), v(good), tGrid(bad), 'nearest', 'extrap');
        else
            v(:) = 0;
        end
    end
    emgAligned(ch,:) = v;
end
end

function s = preprocessEmgSignal(s, fs)
s = s(:)';
s = applyBandpass(s, fs, 20, min(450, fs/2 - 1));
for f0 = [50 100 150 200 250 300 350 400]
    if f0 + 1 < fs/2
        s = applyBandstop(s, fs, f0-1, f0+1);
    end
end
s = abs(s);
s = applyLowpass(s, fs, min(10, fs/4));
end

function y = applyBandpass(x, fs, f1, f2)
if f2 <= f1
    y = x;
    return;
end
try
    y = bandpass(x, [f1 f2], fs);
catch
    [b,a] = butter(4, [f1 f2]/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
end
end

function y = applyBandstop(x, fs, f1, f2)
if f2 <= f1 || f2 >= fs/2
    y = x;
    return;
end
try
    y = bandstop(x, [f1 f2], fs);
catch
    [b,a] = butter(2, [f1 f2]/(fs/2), 'stop');
    y = filtfilt(b,a,x);
end
end

function y = applyLowpass(x, fs, fc)
if fc >= fs/2
    y = x;
    return;
end
try
    y = lowpass(x, fc, fs);
catch
    [b,a] = butter(4, fc/(fs/2), 'low');
    y = filtfilt(b,a,x);
end
end

function [eegTrials, emgTrials, trialInfo] = sliceTrials(eegData, emgAligned, trials)
eegTrials = {};
emgTrials = {};
trialInfo = zeros(0,3); % idx, start, end

for i = 1:size(trials,1)
    i1 = trials(i,1);
    i2 = trials(i,2);
    if i2 <= i1 || i1 < 1 || i2 > size(eegData,2)
        continue;
    end
    e = eegData(:, i1:i2);
    m = emgAligned(:, i1:i2);
    if size(e,2) < 50
        continue;
    end
    if ~all(isfinite(m(:)))
        continue;
    end
    eegTrials{end+1} = e; %#ok<AGROW>
    emgTrials{end+1} = max(m, 0); %#ok<AGROW>
    trialInfo(end+1,:) = [i, i1, i2]; %#ok<AGROW>
end
end

function [eegTrialsRoi, selIdx, selLabels, roiMap] = selectMotorRoiChannels(eegTrials, labels, roiCfg)
% Costruisce time series ROI motorie mediando i canali scalp definiti.
roiMap = struct('name', {}, 'mni_xyz', {}, 'idx', {}, 'labels', {});
for r = 1:numel(roiCfg)
    idx = find(ismember(lower(string(labels)), lower(string(roiCfg(r).labels))));
    if isempty(idx)
        warning('ROI %s senza canali trovati: verrà esclusa.', roiCfg(r).name);
        continue;
    end
    roiMap(end+1).name = roiCfg(r).name; %#ok<AGROW>
    roiMap(end).mni_xyz = roiCfg(r).mni_xyz; %#ok<AGROW>
    roiMap(end).idx = idx(:)'; %#ok<AGROW>
    roiMap(end).labels = labels(idx); %#ok<AGROW>
end

if numel(roiMap) < 2
    error('ROI motorie insufficienti trovate nel file EEG (trovate=%d, richieste>=2).', numel(roiMap));
end

selIdx = cell2mat({roiMap.idx});
selIdx = unique(selIdx, 'stable');
selLabels = {roiMap.name};
eegTrialsRoi = cell(size(eegTrials));
for t = 1:numel(eegTrials)
    X = eegTrials{t};
    Y = zeros(numel(roiMap), size(X,2));
    for r = 1:numel(roiMap)
        Y(r,:) = mean(X(roiMap(r).idx,:), 1);
    end
    eegTrialsRoi{t} = Y;
end
end

function cfg = defaultMotorRoiConfig()
% ROI motorie reaching/hand/upper-limb (MNI, impostazione meta-analitica
% tipo Neurosynth). Manteniamo mot_lpcg/mot_rpcg per compatibilita'
% con la lateralizzazione gia' implementata nello script.
%
% Nota: per scalp-proxy le labels sono una approssimazione EGI256.
cfg = struct('name', {}, 'labels', {}, 'mni_xyz', {});
cfg(1).name = 'mot_lpcg';
cfg(1).mni_xyz = [-38 -24 56];
cfg(1).labels = {'E111','E112','E120','E133','E145'};

cfg(2).name = 'mot_rpcg';
cfg(2).mni_xyz = [38 -24 56];
cfg(2).labels = {'E229','E230','E231','E234','E235'};

cfg(3).name = 'pmd_l';
cfg(3).mni_xyz = [-26 -6 58];
cfg(3).labels = {'E90','E91','E92','E103','E146','E156'};

cfg(4).name = 'pmd_r';
cfg(4).mni_xyz = [26 -6 58];
cfg(4).labels = {'E239','E240','E242','E243','E245','E247','E248','E250','E251','E256'};

cfg(5).name = 'sma_l';
cfg(5).mni_xyz = [-6 -6 58];
cfg(5).labels = {'E173','E174','E187','E188','E197','E198'};

cfg(6).name = 'sma_r';
cfg(6).mni_xyz = [6 -6 58];
cfg(6).labels = {'E208','E216','E217','E226','E227','E236'};

cfg(7).name = 'spl_l';
cfg(7).mni_xyz = [-24 -58 62];
cfg(7).labels = {'E76','E77','E78','E86','E87','E88'};

cfg(8).name = 'spl_r';
cfg(8).mni_xyz = [24 -56 62];
cfg(8).labels = {'E154','E155','E156','E164','E165','E166'};

cfg(9).name = 'ips_l';
cfg(9).mni_xyz = [-34 -46 52];
cfg(9).labels = {'E83','E84','E85','E94','E95','E96'};

cfg(10).name = 'ips_r';
cfg(10).mni_xyz = [34 -44 52];
cfg(10).labels = {'E145','E146','E147','E156','E157','E158'};
end

function s = summarizeRoiMap(roiMap)
if isempty(roiMap)
    s = '';
    return;
end
parts = cell(1, numel(roiMap));
for i = 1:numel(roiMap)
    parts{i} = sprintf('%s[%d]:%s', roiMap(i).name, numel(roiMap(i).idx), strjoin(string(roiMap(i).labels), ';'));
end
s = strjoin(parts, ' | ');
end

function [W, Htrials, nSyn, vafGlobal, vafMinChannel, nmfStability] = extractSynergiesNmf(emgTrials, maxK)
X = [];
trialLen = zeros(1,numel(emgTrials));
for i = 1:numel(emgTrials)
    D = double(emgTrials{i});
    D(D < 0) = 0;
    X = [X D]; %#ok<AGROW>
    trialLen(i) = size(D,2);
end
if isempty(X)
    error('Matrice EMG vuota per NMF');
end

scale = prct(X, 95, 2);
scale(scale <= eps) = 1;
Xn = X ./ scale;
Xn(~isfinite(Xn)) = 0;

best = struct('k',1,'W',[],'H',[],'vafG',-Inf,'vafMin',-Inf);
targetG = 0.90;
targetMin = 0.75;
maxK = min(maxK, size(Xn,1));

for k = 1:maxK
    try
        opts = statset('MaxIter',150,'Display','off');
        [Wk,Hk] = nnmf(Xn, k, 'algorithm','mult', 'replicates',1, 'options',opts);
    catch
        [Wk,Hk] = nnmf(Xn, k, 'algorithm','mult', 'replicates',1);
    end
    Xh = Wk * Hk;
    [vG, vMin] = computeVaf(Xn, Xh);
    if vG > best.vafG
        best = struct('k',k,'W',Wk,'H',Hk,'vafG',vG,'vafMin',vMin);
    end
    if vG >= targetG && vMin >= targetMin
        best = struct('k',k,'W',Wk,'H',Hk,'vafG',vG,'vafMin',vMin);
        break;
    end
end

W = best.W;
H = best.H;
nSyn = best.k;
vafGlobal = best.vafG;
vafMinChannel = best.vafMin;
nmfStability = estimateNmfStability(Xn, nSyn, W, 8);

Htrials = cell(1,numel(trialLen));
p = 1;
for i = 1:numel(trialLen)
    L = trialLen(i);
    Htrials{i} = H(:, p:(p+L-1));
    p = p + L;
end
end

function stab = estimateNmfStability(Xn, k, Wref, nRep)
if isempty(Wref) || k < 1
    stab = NaN;
    return;
end
scores = nan(nRep,1);
W0 = normalizeColumns(Wref);
for r = 1:nRep
    try
        opts = statset('MaxIter',120,'Display','off');
        [Wr,~] = nnmf(Xn, k, 'algorithm','mult', 'replicates',1, 'options',opts);
    catch
        [Wr,~] = nnmf(Xn, k, 'algorithm','mult', 'replicates',1);
    end
    Wr = normalizeColumns(Wr);
    S = abs(W0' * Wr);
    scores(r) = greedyMatchMean(S);
end
stab = meanNoNan(scores);
end

function Wn = normalizeColumns(W)
Wn = W;
for i = 1:size(W,2)
    n = norm(W(:,i));
    if n > eps
        Wn(:,i) = W(:,i) / n;
    end
end
end

function s = greedyMatchMean(S)
Swork = S;
vals = nan(min(size(S)),1);
for i = 1:numel(vals)
    [mx, idx] = max(Swork(:));
    if isempty(mx) || ~isfinite(mx)
        break;
    end
    [r,c] = ind2sub(size(Swork), idx);
    vals(i) = mx;
    Swork(r,:) = -Inf;
    Swork(:,c) = -Inf;
end
s = meanNoNan(vals);
end

function [vG, vMin] = computeVaf(X, Xh)
err = X - Xh;
vG = 1 - sum(err(:).^2) / max(sum(X(:).^2), eps);
vC = zeros(size(X,1),1);
for c = 1:size(X,1)
    xc = X(c,:);
    ec = err(c,:);
    vC(c) = 1 - sum(ec.^2) / max(sum(xc.^2), eps);
end
vMin = min(vC);
end

function [cohMean, f, metrics] = computeCorticoSynergyCoherence(eegTrials, Htrials, fs, channelLabels, emgTrials, opts)
if nargin < 5
    emgTrials = [];
end
if nargin < 6 || isempty(opts)
    opts = struct();
end
basicOnly = isfield(opts,'basicOnly') && opts.basicOnly;

nSyn = size(Htrials{1},1);
nCh = size(eegTrials{1},1);
nTrials = numel(eegTrials);

downStep = max(1, round(fs / 250)); % accelera il calcolo mantenendo banda utile
fsC = fs / downStep;

win = max(128, round(fsC*0.8));
ov = round(0.5*win);
nfft = max(512, 2^nextpow2(win));

cohMean = [];
imCohMean = [];
f = [];
segTotal = 0;
for t = 1:nTrials
    L = size(eegTrials{t},2);
    segTotal = segTotal + max(0, floor((L - ov) / max(win - ov,1)));
end
if segTotal > 1
    c95 = 1 - 0.05^(1/(segTotal - 1));
else
    c95 = NaN;
end

betaMeanTrials = nan(nSyn, nTrials);
imBetaMeanTrials = nan(nSyn, nTrials);
phaseBetaTrials = nan(nSyn, nTrials);
psiTrials = nan(nSyn, nTrials);
lagOptTrialsMs = nan(nSyn, nTrials);
betaLagOptTrials = nan(nSyn, nTrials);
betaWinTrials = nan(nSyn, nTrials, 3); % early/mid/late
liTrials = nan(nSyn, nTrials);         % lateralization index per trial
roiBetaTrials = nan(nCh, nSyn, nTrials);
roiActBetaTrials = nan(nCh, nTrials);  % stato di attivazione ROI (potenza beta EEG)
m1LeftTrials = nan(nSyn, nTrials);
m1RightTrials = nan(nSyn, nTrials);
m1DeltaTrials = nan(nSyn, nTrials);
m1DeltaWinTrials = nan(nSyn, nTrials, 3); % early/mid/late
m1DeltaLowBetaTrials = nan(nSyn, nTrials);  % 13-20 Hz
m1DeltaHighBetaTrials = nan(nSyn, nTrials); % 20-30 Hz

% ROI indices for lateralization
leftIdx = [];
rightIdx = [];
if nargin >= 4 && ~isempty(channelLabels)
    labs = lower(string(channelLabels));
    leftIdx = find(contains(labs,'mot_lpcg'));
    rightIdx = find(contains(labs,'mot_rpcg'));
end

for k = 1:nSyn
    trialCurves = {};
    trialImCurves = {};
    for t = 1:nTrials
        perChC = [];
        perChI = [];
        psiBetaCh = nan(1,nCh);
        phaseBetaCh = nan(1,nCh);
        winBetaCh = nan(nCh,3);
        lowBetaCh = nan(1,nCh);
        highBetaCh = nan(1,nCh);
        liLeft = [];
        liRight = [];
        eAvgForLag = [];
        hForLag = [];
        for ch = 1:nCh
            e = double(eegTrials{t}(ch,1:downStep:end));
            h = double(Htrials{t}(k,1:downStep:end));
            L = min(numel(e), numel(h));
            if L < win
                continue;
            end
            e = e(1:L);
            h = h(1:L);
            [cxy, icxy, ph, ff, cohComp] = computeCoherencePair(e, h, win, ov, nfft, fsC);
            if k == 1 && ~isfinite(roiActBetaTrials(ch,t))
                [Pxx, ffP] = safePwelch(e, win, ov, nfft, fsC);
                roiActBetaTrials(ch,t) = bandMeanFromCurve(Pxx(:), ffP(:), [13 30]);
            end
            if isempty(f)
                f = ff(:);
            end
            perChC = [perChC cxy(:)]; %#ok<AGROW>
            perChI = [perChI icxy(:)]; %#ok<AGROW>
            psiBetaCh(ch) = computePsi(cohComp, ff(:), [13 30]);
            [bm, ~, ~] = betaStatsFromCurve(cxy(:), ff(:), [13 30]);
            lowBetaCh(ch) = bandMeanFromCurve(cxy(:), ff(:), [13 20]);
            highBetaCh(ch) = bandMeanFromCurve(cxy(:), ff(:), [20 30]);
            if isfinite(bm)
                roiBetaTrials(ch,k,t) = bm;
            end
            if ismember(ch, leftIdx), liLeft(end+1) = bm; end %#ok<AGROW>
            if ismember(ch, rightIdx), liRight(end+1) = bm; end %#ok<AGROW>
            phaseBetaCh(ch) = circularMean(ph(ff>=13 & ff<=30));
            eAvgForLag = [eAvgForLag; e(:)']; %#ok<AGROW>
            hForLag = h;

            % Gesture windows: split trial in 3 thirds
            [w1, w2, w3] = splitThirds(numel(e));
            winBetaCh(ch,1) = computeWindowBeta(e(w1), h(w1), win, ov, nfft, fsC);
            winBetaCh(ch,2) = computeWindowBeta(e(w2), h(w2), win, ov, nfft, fsC);
            winBetaCh(ch,3) = computeWindowBeta(e(w3), h(w3), win, ov, nfft, fsC);
        end

        if ~isempty(perChC)
            cTrial = mean(perChC,2);
            iTrial = mean(perChI,2);
            trialCurves{end+1} = cTrial; %#ok<AGROW>
            trialImCurves{end+1} = iTrial; %#ok<AGROW>
            [betaMeanTrials(k,t), ~, ~] = betaStatsFromCurve(cTrial, f, [13 30]);
            [imBetaMeanTrials(k,t), ~, ~] = betaStatsFromCurve(iTrial, f, [13 30]);
            phaseBetaTrials(k,t) = circularMean(phaseBetaCh(isfinite(phaseBetaCh)));
            psiTrials(k,t) = meanNoNan(psiBetaCh);
            betaWinTrials(k,t,1) = meanNoNan(winBetaCh(:,1));
            betaWinTrials(k,t,2) = meanNoNan(winBetaCh(:,2));
            betaWinTrials(k,t,3) = meanNoNan(winBetaCh(:,3));
            if ~isempty(eAvgForLag) && ~isempty(hForLag)
                eLag = mean(eAvgForLag, 1);
                [lagMs, betaLag] = computeLagOptimizedBeta(eLag, hForLag, fsC, win, ov, nfft, 50);
                lagOptTrialsMs(k,t) = lagMs;
                betaLagOptTrials(k,t) = betaLag;
            end
            if ~isempty(liLeft) && ~isempty(liRight)
                Lm = mean(liLeft);
                Rm = mean(liRight);
                liTrials(k,t) = (Rm - Lm) / max(Rm + Lm, eps);
            end
            if ~isempty(leftIdx) && ~isempty(rightIdx)
                m1L = meanNoNan(roiBetaTrials(leftIdx,k,t));
                m1R = meanNoNan(roiBetaTrials(rightIdx,k,t));
                m1LeftTrials(k,t) = m1L;
                m1RightTrials(k,t) = m1R;
                if isfinite(m1L) && isfinite(m1R)
                    m1DeltaTrials(k,t) = m1L - m1R;
                end
                m1eL = meanNoNan(winBetaCh(leftIdx,1));
                m1eR = meanNoNan(winBetaCh(rightIdx,1));
                m1mL = meanNoNan(winBetaCh(leftIdx,2));
                m1mR = meanNoNan(winBetaCh(rightIdx,2));
                m1lL = meanNoNan(winBetaCh(leftIdx,3));
                m1lR = meanNoNan(winBetaCh(rightIdx,3));
                m1DeltaWinTrials(k,t,1) = m1eL - m1eR;
                m1DeltaWinTrials(k,t,2) = m1mL - m1mR;
                m1DeltaWinTrials(k,t,3) = m1lL - m1lR;
                m1DeltaLowBetaTrials(k,t) = meanNoNan(lowBetaCh(leftIdx)) - meanNoNan(lowBetaCh(rightIdx));
                m1DeltaHighBetaTrials(k,t) = meanNoNan(highBetaCh(leftIdx)) - meanNoNan(highBetaCh(rightIdx));
            end
        end
    end

    if isempty(trialCurves)
        cohMean(k,:) = zeros(1, numel(f)); %#ok<AGROW>
        imCohMean(k,:) = zeros(1, numel(f)); %#ok<AGROW>
    else
        cohMean(k,:) = mean(cell2mat(trialCurves),2)'; %#ok<AGROW>
        imCohMean(k,:) = mean(cell2mat(trialImCurves),2)'; %#ok<AGROW>
    end
end

metrics = struct();
metrics.c95 = c95;
metrics.peakCoh = zeros(nSyn,1);
metrics.peakHz = zeros(nSyn,1);
metrics.meanBand = zeros(nSyn,1);
metrics.betaPeakCoh = zeros(nSyn,1);
metrics.betaPeakHz = zeros(nSyn,1);
metrics.betaMean = zeros(nSyn,1);
metrics.betaAuc = zeros(nSyn,1);
metrics.betaMeanLagOpt = zeros(nSyn,1);
metrics.lagOptMs = zeros(nSyn,1);
metrics.psiBeta = zeros(nSyn,1);
metrics.imBetaMean = zeros(nSyn,1);
metrics.phaseBetaMeanRad = zeros(nSyn,1);
metrics.betaMeanCiLow = zeros(nSyn,1);
metrics.betaMeanCiHigh = zeros(nSyn,1);
metrics.betaStabilityStd = zeros(nSyn,1);
metrics.betaStabilityCv = zeros(nSyn,1);
metrics.betaEarly = zeros(nSyn,1);
metrics.betaMiddle = zeros(nSyn,1);
metrics.betaLate = zeros(nSyn,1);
metrics.gestureDeltaLateEarly = zeros(nSyn,1);
metrics.lateralizationLi = zeros(nSyn,1);
metrics.m1LeftBeta = zeros(nSyn,1);
metrics.m1RightBeta = zeros(nSyn,1);
metrics.m1DeltaLminusR = zeros(nSyn,1);
metrics.m1DeltaLminusREarly = zeros(nSyn,1);
metrics.m1DeltaLminusRMiddle = zeros(nSyn,1);
metrics.m1DeltaLminusRLate = zeros(nSyn,1);
metrics.m1DeltaLminusRLowBeta = zeros(nSyn,1);
metrics.m1DeltaLminusRHighBeta = zeros(nSyn,1);
metrics.surrogateP = nan(nSyn,1);
metrics.surrogateQ = nan(nSyn,1);
metrics.surrogateSig = false(nSyn,1);
metrics.zCohSurrogateP = nan(nSyn,1);
metrics.zCohSurrogateQ = nan(nSyn,1);
metrics.zCohSurrogateSig = false(nSyn,1);

band = (f >= 1 & f <= 45);
beta = (f >= 13 & f <= 30);
for k = 1:nSyn
    ck = cohMean(k,:);
    if any(band)
        [pk, idx] = max(ck(band));
        fb = f(band);
        metrics.peakCoh(k) = pk;
        metrics.peakHz(k) = fb(idx);
        metrics.meanBand(k) = mean(ck(band));
    else
        metrics.peakCoh(k) = NaN;
        metrics.peakHz(k) = NaN;
        metrics.meanBand(k) = NaN;
    end
    if any(beta)
        cb = ck(beta);
        fb = f(beta);
        [pkb, idxb] = max(cb);
        metrics.betaPeakCoh(k) = pkb;
        metrics.betaPeakHz(k) = fb(idxb);
        metrics.betaMean(k) = mean(cb);
        metrics.betaAuc(k) = trapz(fb, cb);
        metrics.imBetaMean(k) = mean(imCohMean(k,beta));
    else
        metrics.betaPeakCoh(k) = NaN;
        metrics.betaPeakHz(k) = NaN;
        metrics.betaMean(k) = NaN;
        metrics.betaAuc(k) = NaN;
        metrics.imBetaMean(k) = NaN;
    end

    bt = betaMeanTrials(k,:);
    bt = bt(isfinite(bt));
    [lo, hi] = bootstrapMeanCI(bt, 500, 0.05);
    metrics.betaMeanCiLow(k) = lo;
    metrics.betaMeanCiHigh(k) = hi;
    metrics.betaStabilityStd(k) = std(bt, 0, 'omitnan');
    metrics.betaStabilityCv(k) = metrics.betaStabilityStd(k) / max(mean(bt,'omitnan'), eps);
    metrics.phaseBetaMeanRad(k) = circularMean(phaseBetaTrials(k,isfinite(phaseBetaTrials(k,:))));
    metrics.psiBeta(k) = meanNoNan(psiTrials(k,:));
    metrics.lagOptMs(k) = meanNoNan(lagOptTrialsMs(k,:));
    metrics.betaMeanLagOpt(k) = meanNoNan(betaLagOptTrials(k,:));
    metrics.betaEarly(k) = meanNoNan(squeeze(betaWinTrials(k,:,1)));
    metrics.betaMiddle(k) = meanNoNan(squeeze(betaWinTrials(k,:,2)));
    metrics.betaLate(k) = meanNoNan(squeeze(betaWinTrials(k,:,3)));
    metrics.gestureDeltaLateEarly(k) = metrics.betaLate(k) - metrics.betaEarly(k);
    metrics.lateralizationLi(k) = meanNoNan(liTrials(k,:));
    metrics.m1LeftBeta(k) = meanNoNan(m1LeftTrials(k,:));
    metrics.m1RightBeta(k) = meanNoNan(m1RightTrials(k,:));
    metrics.m1DeltaLminusR(k) = meanNoNan(m1DeltaTrials(k,:));
    metrics.m1DeltaLminusREarly(k) = meanNoNan(squeeze(m1DeltaWinTrials(k,:,1)));
    metrics.m1DeltaLminusRMiddle(k) = meanNoNan(squeeze(m1DeltaWinTrials(k,:,2)));
    metrics.m1DeltaLminusRLate(k) = meanNoNan(squeeze(m1DeltaWinTrials(k,:,3)));
    metrics.m1DeltaLminusRLowBeta(k) = meanNoNan(m1DeltaLowBetaTrials(k,:));
    metrics.m1DeltaLminusRHighBeta(k) = meanNoNan(m1DeltaHighBetaTrials(k,:));
end

% Surrogati trial-shuffle + FDR sulle sinergie
if ~basicOnly
    surrN = 40;
    pvals = nan(nSyn,1);
    pvalsZ = nan(nSyn,1);
    for k = 1:nSyn
        obs = metrics.betaMean(k);
        obsZ = zcohBetaMeanFromCurve(cohMean(k,:), f, [13 30]);
        if ~isfinite(obs) && ~isfinite(obsZ)
            continue;
        end
        [nullBeta, nullZ] = surrogateBetaAndZcohNull(eegTrials, Htrials, k, downStep, win, ov, nfft, fsC, surrN);
        if isfinite(obs)
            nb = nullBeta(isfinite(nullBeta));
            if ~isempty(nb)
                pvals(k) = (1 + sum(nb >= obs)) / (numel(nb) + 1);
            end
        end
        if isfinite(obsZ)
            nz = nullZ(isfinite(nullZ));
            if ~isempty(nz)
                pvalsZ(k) = (1 + sum(nz >= obsZ)) / (numel(nz) + 1);
            end
        end
    end
    qvals = bhFdr(pvals);
    qvalsZ = bhFdr(pvalsZ);
    metrics.surrogateP = pvals;
    metrics.surrogateQ = qvals;
    metrics.surrogateSig = qvals < 0.05;
    metrics.zCohSurrogateP = pvalsZ;
    metrics.zCohSurrogateQ = qvalsZ;
    metrics.zCohSurrogateSig = qvalsZ < 0.05;
end

% Accoppiamento multivariato ROI<->sinergie (CCA)
if ~basicOnly
    [metrics.multivarR1, metrics.multivarR1PermP] = computeMultivarCoupling(eegTrials, Htrials, downStep, 100);
else
    metrics.multivarR1 = NaN;
    metrics.multivarR1PermP = NaN;
end

% Covariate trial-level (proxy comportamentali)
if ~basicOnly
    [metrics.covarCorrEmgDriveBeta, metrics.covarCorrDurationBeta] = ...
        computeTrialCovariateCorrelations(betaMeanTrials, Htrials, fs);
else
    metrics.covarCorrEmgDriveBeta = NaN;
    metrics.covarCorrDurationBeta = NaN;
end

metrics.betaMeanAll = meanNoNan(metrics.betaMean);
metrics.betaPeakAll = maxNoNan(metrics.betaPeakCoh);
metrics.betaAucAll = meanNoNan(metrics.betaAuc);
metrics.betaMeanLagOptAll = meanNoNan(metrics.betaMeanLagOpt);
metrics.lagOptMsAll = meanNoNan(metrics.lagOptMs);
metrics.psiBetaAll = meanNoNan(metrics.psiBeta);
metrics.imBetaMeanAll = meanNoNan(metrics.imBetaMean);
metrics.betaCi95LowAll = meanNoNan(metrics.betaMeanCiLow);
metrics.betaCi95HighAll = meanNoNan(metrics.betaMeanCiHigh);
metrics.betaStabilityCvAll = meanNoNan(metrics.betaStabilityCv);
metrics.gestureDeltaLateEarlyAll = meanNoNan(metrics.gestureDeltaLateEarly);
metrics.lateralizationLiAll = meanNoNan(metrics.lateralizationLi);
metrics.m1LeftBetaAll = meanNoNan(metrics.m1LeftBeta);
metrics.m1RightBetaAll = meanNoNan(metrics.m1RightBeta);
metrics.m1DeltaLminusRAll = meanNoNan(metrics.m1DeltaLminusR);
metrics.m1DeltaLminusREarlyAll = meanNoNan(metrics.m1DeltaLminusREarly);
metrics.m1DeltaLminusRMiddleAll = meanNoNan(metrics.m1DeltaLminusRMiddle);
metrics.m1DeltaLminusRLateAll = meanNoNan(metrics.m1DeltaLminusRLate);
metrics.m1DeltaLminusRLowBetaAll = meanNoNan(metrics.m1DeltaLminusRLowBeta);
metrics.m1DeltaLminusRHighBetaAll = meanNoNan(metrics.m1DeltaLminusRHighBeta);
metrics.surrogatePMin = minNoNan(metrics.surrogateP);
metrics.surrogateQMin = minNoNan(metrics.surrogateQ);
metrics.surrogateSigCount = sum(metrics.surrogateSig);
metrics.zCohSurrogatePMin = minNoNan(metrics.zCohSurrogateP);
metrics.zCohSurrogateQMin = minNoNan(metrics.zCohSurrogateQ);
metrics.zCohSurrogateSigCount = sum(metrics.zCohSurrogateSig);
if all(~isfinite(metrics.betaPeakCoh))
    metrics.betaPeakHzAll = NaN;
else
    [~, ii] = max(metrics.betaPeakCoh);
    metrics.betaPeakHzAll = metrics.betaPeakHz(ii);
end

% Metriche per-ROI (ranking): media beta per ROI su sinergie/trial
if nargin >= 4 && ~isempty(channelLabels) && numel(channelLabels) == nCh
    metrics.roiLabels = cellstr(string(channelLabels(:)));
else
    metrics.roiLabels = arrayfun(@(k) sprintf('ROI_%02d',k), 1:nCh, 'UniformOutput', false)';
end
metrics.roiBetaMeanBySynergy = nan(nCh, nSyn);
for ch = 1:nCh
    for k = 1:nSyn
        metrics.roiBetaMeanBySynergy(ch,k) = meanNoNan(squeeze(roiBetaTrials(ch,k,:)));
    end
end
metrics.roiBetaMeanAll = nan(nCh,1);
for ch = 1:nCh
    metrics.roiBetaMeanAll(ch) = meanNoNan(metrics.roiBetaMeanBySynergy(ch,:));
end
[mx, ix] = max(metrics.roiBetaMeanAll);
if isempty(ix) || ~isfinite(mx)
    metrics.roiTopLabel = '';
    metrics.roiTopBetaMean = NaN;
else
    metrics.roiTopLabel = metrics.roiLabels{ix};
    metrics.roiTopBetaMean = mx;
end

% State-dependent CSC: confronto ROI attiva vs inattiva (soglia percentile su potenza beta ROI)
metrics.roiActivityPercentile = 70;
metrics.roiActivityThreshold = nan(nCh,1);
metrics.roiStateNActive = zeros(nCh,1);
metrics.roiStateNInactive = zeros(nCh,1);
metrics.roiStateCscActive = nan(nCh, nSyn);
metrics.roiStateCscInactive = nan(nCh, nSyn);
metrics.roiStateDelta = nan(nCh, nSyn);
metrics.roiStateRatio = nan(nCh, nSyn);
metrics.roiStateDeltaAll = nan(nCh,1);
metrics.roiStateRatioAll = nan(nCh,1);
metrics.roiStateTopLabel = '';
metrics.roiStateTopDelta = NaN;

for ch = 1:nCh
    a = roiActBetaTrials(ch,:);
    valid = isfinite(a);
    if sum(valid) < 6
        continue;
    end
    thr = prct(a(valid), metrics.roiActivityPercentile, 2);
    metrics.roiActivityThreshold(ch) = thr;
    maskAct = a >= thr;
    maskIn = ~maskAct & valid;
    metrics.roiStateNActive(ch) = sum(maskAct);
    metrics.roiStateNInactive(ch) = sum(maskIn);
    if metrics.roiStateNActive(ch) < 3 || metrics.roiStateNInactive(ch) < 3
        continue;
    end
    for k = 1:nSyn
        v = squeeze(roiBetaTrials(ch,k,:))';
        vA = v(maskAct);
        vI = v(maskIn);
        cA = meanNoNan(vA);
        cI = meanNoNan(vI);
        metrics.roiStateCscActive(ch,k) = cA;
        metrics.roiStateCscInactive(ch,k) = cI;
        metrics.roiStateDelta(ch,k) = cA - cI;
        if isfinite(cA) && isfinite(cI) && abs(cI) > eps
            metrics.roiStateRatio(ch,k) = cA / cI;
        end
    end
    metrics.roiStateDeltaAll(ch) = meanNoNan(metrics.roiStateDelta(ch,:));
    metrics.roiStateRatioAll(ch) = meanNoNan(metrics.roiStateRatio(ch,:));
end
[mxDelta, iTop] = max(metrics.roiStateDeltaAll);
if ~isempty(iTop) && isfinite(mxDelta)
    metrics.roiStateTopLabel = metrics.roiLabels{iTop};
    metrics.roiStateTopDelta = mxDelta;
end
end

function [cxy, icxy, phxy, ff, cohComp] = computeCoherencePair(e, h, win, ov, nfft, fsC)
[cxy, ff] = safeMscohere(e, h, win, ov, nfft, fsC);
[Pxy, ~] = safeCpsd(e, h, win, ov, nfft, fsC);
[Pxx, ~] = safePwelch(e, win, ov, nfft, fsC);
[Pyy, ~] = safePwelch(h, win, ov, nfft, fsC);
icxy = abs(imag(Pxy)) ./ sqrt(max(Pxx.*Pyy, eps));
phxy = angle(Pxy);
cohComp = Pxy ./ sqrt(max(Pxx.*Pyy, eps));
end

function [m, pk, fpk] = betaStatsFromCurve(curve, ff, brange)
beta = ff >= brange(1) & ff <= brange(2);
if ~any(beta)
    m = NaN; pk = NaN; fpk = NaN; return;
end
cb = curve(beta);
fb = ff(beta);
[pk, ii] = max(cb);
fpk = fb(ii);
m = mean(cb);
end

function zbm = zcohBetaMeanFromCurve(curve, ff, brange)
beta = ff >= brange(1) & ff <= brange(2);
if ~any(beta)
    zbm = NaN;
    return;
end
zc = fisherZcoh(curve(beta));
zbm = meanNoNan(zc);
end

function m = bandMeanFromCurve(curve, ff, brange)
idx = ff >= brange(1) & ff <= brange(2);
if ~any(idx)
    m = NaN;
else
    m = mean(curve(idx));
end
end

function v = computeWindowBeta(e, h, win, ov, nfft, fsC)
L = min(numel(e), numel(h));
if L < win
    v = NaN; return;
end
[cxy, ff] = safeMscohere(e(1:L), h(1:L), win, ov, nfft, fsC);
[v, ~, ~] = betaStatsFromCurve(cxy(:), ff(:), [13 30]);
end

function [w1, w2, w3] = splitThirds(L)
e1 = floor(L/3);
e2 = floor(2*L/3);
w1 = 1:max(e1,1);
w2 = max(e1+1,1):max(e2,1);
w3 = max(e2+1,1):L;
if isempty(w2), w2 = w1; end
if isempty(w3), w3 = w2; end
end

function mu = circularMean(x)
x = x(isfinite(x));
if isempty(x)
    mu = NaN;
else
    mu = angle(mean(exp(1i*x)));
end
end

function [lo, hi] = bootstrapMeanCI(x, nboot, alpha)
x = x(isfinite(x));
if isempty(x)
    lo = NaN; hi = NaN; return;
end
if numel(x) == 1
    lo = x; hi = x; return;
end
nb = max(100, nboot);
m = zeros(nb,1);
n = numel(x);
for b = 1:nb
    idx = randi(n, [n 1]);
    m(b) = mean(x(idx));
end
m = sort(m);
ilo = max(1, floor((alpha/2)*nb));
ihi = min(nb, ceil((1-alpha/2)*nb));
lo = m(ilo);
hi = m(ihi);
end

function psi = computePsi(cohComp, ff, brange)
beta = ff >= brange(1) & ff <= brange(2);
idx = find(beta);
if numel(idx) < 2
    psi = NaN;
    return;
end
C = cohComp(idx);
psi = sum(imag(conj(C(1:end-1)) .* C(2:end)));
end

function [lagMs, bestBeta] = computeLagOptimizedBeta(e, h, fsC, win, ov, nfft, maxLagMs)
if isempty(e) || isempty(h)
    lagMs = NaN; bestBeta = NaN; return;
end
maxLag = max(1, round((maxLagMs/1000) * fsC));
lags = -maxLag:maxLag;
vals = nan(size(lags));
for i = 1:numel(lags)
    lg = lags(i);
    if lg >= 0
        ee = e(1:end-lg);
        hh = h(1+lg:end);
    else
        ee = e(1-lg:end);
        hh = h(1:end+lg);
    end
    L = min(numel(ee), numel(hh));
    if L < win
        continue;
    end
    [cxy, ff] = safeMscohere(ee(1:L), hh(1:L), win, ov, nfft, fsC);
    [bm,~,~] = betaStatsFromCurve(cxy(:), ff(:), [13 30]);
    vals(i) = bm;
end
[bestBeta, ii] = max(vals);
if isempty(ii) || ~isfinite(bestBeta)
    lagMs = NaN;
else
    lagMs = (lags(ii) / fsC) * 1000;
end
end

function [nullBeta, nullZ] = surrogateBetaAndZcohNull(eegTrials, Htrials, k, downStep, win, ov, nfft, fsC, surrN)
nTrials = numel(eegTrials);
nullBeta = nan(surrN,1);
nullZ = nan(surrN,1);
eTrials = cell(1,nTrials);
hTrials = cell(1,nTrials);
for t = 1:nTrials
    e = mean(double(eegTrials{t}(:,1:downStep:end)),1);
    h = double(Htrials{t}(k,1:downStep:end));
    L = min(numel(e), numel(h));
    eTrials{t} = e(1:L);
    hTrials{t} = h(1:L);
end
for s = 1:surrN
    bt = nan(1,nTrials);
    zt = nan(1,nTrials);
    for t = 1:nTrials
        e = eTrials{t};
        h = hTrials{t};
        L = min(numel(e), numel(h));
        if L < win
            continue;
        end
        sh = randi([max(1,round(0.1*L)) max(1,round(0.9*L))],1,1);
        h = circshift(h, [0 sh]);
        [cxy, ff] = safeMscohere(e(1:L), h(1:L), win, ov, nfft, fsC);
        [bt(t),~,~] = betaStatsFromCurve(cxy(:), ff(:), [13 30]);
        zt(t) = zcohBetaMeanFromCurve(cxy(:), ff(:), [13 30]);
    end
    nullBeta(s) = meanNoNan(bt);
    nullZ(s) = meanNoNan(zt);
end
end

function [Pxx, ff] = safePwelch(x, win, ov, nfft, fsC)
try
    [Pxx, ff] = pwelch(x, win, ov, nfft, fsC);
catch ME
    if ~isOverlapFractionError(ME)
        rethrow(ME);
    end
    [Pxx, ff] = pwelch(x, win, overlapFraction(ov, win), nfft, fsC);
end
end

function [Pxy, ff] = safeCpsd(x, y, win, ov, nfft, fsC)
try
    [Pxy, ff] = cpsd(x, y, win, ov, nfft, fsC);
catch ME
    if ~isOverlapFractionError(ME)
        rethrow(ME);
    end
    [Pxy, ff] = cpsd(x, y, win, overlapFraction(ov, win), nfft, fsC);
end
end

function [cxy, ff] = safeMscohere(x, y, win, ov, nfft, fsC)
try
    [cxy, ff] = mscohere(x, y, win, ov, nfft, fsC);
catch ME
    if ~isOverlapFractionError(ME)
        rethrow(ME);
    end
    [cxy, ff] = mscohere(x, y, win, overlapFraction(ov, win), nfft, fsC);
end
end

function tf = isOverlapFractionError(ME)
msg = lower(ME.message);
tf = contains(msg, 'overlap') && ...
    (contains(msg, '0.950000') || contains(msg, '0 to 0.95') || contains(msg, '0 to 0.950000'));
end

function ovFrac = overlapFraction(ov, win)
if numel(win) > 1
    winLen = numel(win);
else
    winLen = double(win);
end
if ~isfinite(winLen) || winLen <= 0
    ovFrac = 0;
else
    ovFrac = ov / winLen;
end
ovFrac = min(0.95, max(0, ovFrac));
end

function q = bhFdr(p)
p = p(:);
q = nan(size(p));
ok = isfinite(p);
if ~any(ok)
    return;
end
pv = p(ok);
[ps, ord] = sort(pv);
m = numel(ps);
qs = ps .* m ./ (1:m)';
for i = m-1:-1:1
    qs(i) = min(qs(i), qs(i+1));
end
qs = min(qs, 1);
tmp = nan(size(pv));
tmp(ord) = qs;
q(ok) = tmp;
end

function [r1, pperm] = computeMultivarCoupling(eegTrials, Htrials, downStep, nPerm)
X = [];
Y = [];
for t = 1:numel(eegTrials)
    E = double(eegTrials{t}(:,1:downStep:end))';
    H = double(Htrials{t}(:,1:downStep:end))';
    L = min(size(E,1), size(H,1));
    if L < 20
        continue;
    end
    X = [X; E(1:L,:)]; %#ok<AGROW>
    Y = [Y; H(1:L,:)]; %#ok<AGROW>
end
if size(X,1) < 40 || size(X,2) < 2 || size(Y,2) < 2
    r1 = NaN; pperm = NaN; return;
end
X = zscore(X);
Y = zscore(Y);
try
    [~,~,r] = canoncorr(X,Y);
    r1 = r(1);
catch
    C = corrcoef([X Y]);
    r1 = max(abs(C(1:size(X,2), size(X,2)+1:end)), [], 'all');
end
if ~isfinite(r1)
    pperm = NaN;
    return;
end
rp = nan(nPerm,1);
N = size(X,1);
for i = 1:nPerm
    idx = randperm(N);
    Yp = Y(idx,:);
    try
        [~,~,r] = canoncorr(X,Yp);
        rp(i) = r(1);
    catch
        C = corrcoef([X Yp]);
        rp(i) = max(abs(C(1:size(X,2), size(X,2)+1:end)), [], 'all');
    end
end
pperm = (1 + sum(rp >= r1)) / (nPerm + 1);
end

function [rEmg, rDur] = computeTrialCovariateCorrelations(betaMeanTrials, Htrials, fs)
bt = mean(betaMeanTrials,1,'omitnan')';
nT = numel(Htrials);
emgDrive = nan(nT,1);
dur = nan(nT,1);
for t = 1:nT
    H = Htrials{t};
    emgDrive(t) = mean(H(:), 'omitnan');
    dur(t) = size(H,2) / fs;
end
[rEmg,~] = corrSafe(bt, emgDrive);
[rDur,~] = corrSafe(bt, dur);
end

function [r,p] = corrSafe(x,y)
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 6
    r = NaN; p = NaN; return;
end
try
    [r,p] = corr(x(ok), y(ok), 'Type', 'Spearman');
catch
    C = corrcoef(x(ok), y(ok));
    r = C(1,2);
    p = NaN;
end
end


function metrics = addZcohAndSynergyClusters(metrics, cohMean, f)
nSyn = size(cohMean,1);
metrics.zCohCurve = nan(size(cohMean));
metrics.zCohThreshold = nan(nSyn,1);
metrics.zCohBetaMean = nan(nSyn,1);
metrics.zCohBetaPeak = nan(nSyn,1);
metrics.zCohBetaPeakHz = nan(nSyn,1);
metrics.zCohSigFracBeta = nan(nSyn,1);
metrics.zCohBetaMeanAll = NaN;
metrics.zCohBetaPeakAll = NaN;
metrics.zCohBetaPeakHzAll = NaN;
metrics.zCohThresholdBetaAll = NaN;
metrics.zCohSigFracBetaAll = NaN;
metrics.synergyClusterId = nan(nSyn,1);
metrics.synergyClusterK = NaN;
metrics.synergyClusterQuality = NaN;
metrics.synergyClusterCentroids = [];

if isempty(cohMean) || isempty(f)
    return;
end

band = (f >= 1 & f <= 45);
beta = (f >= 13 & f <= 30);
for k = 1:nSyn
    zc = fisherZcoh(cohMean(k,:));
    metrics.zCohCurve(k,:) = zc;
    if any(band)
        metrics.zCohThreshold(k) = prct(zc(band), 95, 2);
    end
    if any(beta)
        zb = zc(beta);
        fb = f(beta);
        metrics.zCohBetaMean(k) = meanNoNan(zb);
        [pk, ii] = max(zb);
        if ~isempty(ii) && isfinite(pk)
            metrics.zCohBetaPeak(k) = pk;
            metrics.zCohBetaPeakHz(k) = fb(ii);
        end
        thr = metrics.zCohThreshold(k);
        ok = isfinite(zb);
        if any(ok) && isfinite(thr)
            metrics.zCohSigFracBeta(k) = mean(zb(ok) > thr);
        end
    end
end

metrics.zCohBetaMeanAll = meanNoNan(metrics.zCohBetaMean);
metrics.zCohBetaPeakAll = maxNoNan(metrics.zCohBetaPeak);
if all(~isfinite(metrics.zCohBetaPeak))
    metrics.zCohBetaPeakHzAll = NaN;
else
    [~, ii] = max(metrics.zCohBetaPeak);
    metrics.zCohBetaPeakHzAll = metrics.zCohBetaPeakHz(ii);
end
metrics.zCohThresholdBetaAll = meanNoNan(metrics.zCohThreshold);
metrics.zCohSigFracBetaAll = meanNoNan(metrics.zCohSigFracBeta);

[clId, kBest, qBest, cBest] = clusterSynergiesByProfile(metrics.zCohCurve, f, 5);
metrics.synergyClusterId = clId;
metrics.synergyClusterK = kBest;
metrics.synergyClusterQuality = qBest;
metrics.synergyClusterCentroids = cBest;
end

function zc = fisherZcoh(c)
c = max(min(c, 1 - 1e-9), 0);
zc = atanh(sqrt(c));
end

function [labels, kBest, qBest, centroids] = clusterSynergiesByProfile(curves, f, maxK)
n = size(curves,1);
labels = ones(n,1);
kBest = 1;
qBest = NaN;
centroids = [];
if n < 2 || isempty(curves) || isempty(f)
    return;
end

idx = (f >= 10 & f <= 28);
if nnz(idx) < 4
    idx = (f >= 13 & f <= 30);
end
if nnz(idx) < 2
    return;
end

X = curves(:, idx);
X = fillmissing(X, 'constant', 0);
X = bsxfun(@minus, X, mean(X,2));
X = bsxfun(@rdivide, X, sqrt(sum(X.^2,2)) + eps);

kMax = min(maxK, n);
if kMax < 2
    return;
end

D = corrDistanceMatrix(X);
bestLoss = Inf;
qBest = -Inf;
bestL = labels;
bestC = mean(X,1);
for k = 2:kMax
    [lab, C, loss] = kmeansCorrSimple(X, k, 25, 80);
    q = silhouetteFromDistance(D, lab);
    if isfinite(q) && (q > qBest || (abs(q-qBest) < 1e-9 && loss < bestLoss))
        qBest = q;
        bestLoss = loss;
        bestL = lab;
        bestC = C;
        kBest = k;
    end
end
labels = bestL;
centroids = bestC;
if ~isfinite(qBest)
    qBest = NaN;
end
end

function D = corrDistanceMatrix(X)
S = X * X';
D = 1 - S;
D(1:size(D,1)+1:end) = 0;
D = max(D, 0);
end

function q = silhouetteFromDistance(D, labels)
n = numel(labels);
if n < 3
    q = NaN;
    return;
end
s = nan(n,1);
for i = 1:n
    ci = labels(i);
    same = find(labels == ci);
    same(same == i) = [];
    if isempty(same)
        continue;
    end
    a = mean(D(i, same));
    b = Inf;
    cls = unique(labels(:)');
    for c = cls
        if c == ci
            continue;
        end
        oth = find(labels == c);
        if isempty(oth)
            continue;
        end
        b = min(b, mean(D(i, oth)));
    end
    if isfinite(a) && isfinite(b) && max(a,b) > eps
        s(i) = (b - a) / max(a, b);
    end
end
q = meanNoNan(s);
end

function [labels, C, bestLoss] = kmeansCorrSimple(X, k, nRep, nIter)
n = size(X,1);
bestLoss = Inf;
labels = ones(n,1);
C = X(1:min(k,n),:);
if n < k
    return;
end
for r = 1:nRep
    idx = randperm(n, k);
    Ccur = X(idx,:);
    lcur = ones(n,1);
    for it = 1:nIter
        dist = 1 - (X * Ccur');
        [~, lnew] = min(dist, [], 2);
        if it > 1 && all(lnew == lcur)
            break;
        end
        lcur = lnew;
        for j = 1:k
            ii = find(lcur == j);
            if isempty(ii)
                Ccur(j,:) = X(randi(n),:);
            else
                c = mean(X(ii,:), 1);
                c = c / (norm(c) + eps);
                Ccur(j,:) = c;
            end
        end
    end
    dist = 1 - (X * Ccur');
    lin = sub2ind(size(dist), (1:n)', lcur);
    loss = sum(dist(lin));
    if loss < bestLoss
        bestLoss = loss;
        labels = lcur;
        C = Ccur;
    end
end
end

function makeSynergyPlot(W, outDir)
fig = figure('Visible','off','Color','w');
imagesc(W);
axis tight;
xlabel('Synergy');
ylabel('EMG Bipolar Channels');
title('W - Pesi Sinergie (10 bipolari)');
colorbar;
saveas(fig, fullfile(outDir, 'sinergie_W.png'));
close(fig);
end

function makeCoherenceChangePlot(cohMean, f, metrics, outDir)
nSyn = size(cohMean,1);
if nSyn < 1 || isempty(f)
    return;
end

band = (f >= 1 & f <= 45);
if ~any(band)
    return;
end
fb = f(band);
C = cohMean(:, band);
Cd = C - mean(C, 1, 'omitnan');
betaBand = (fb >= 13 & fb <= 30);
if ~any(betaBand)
    betaBand = true(size(fb));
end

% Ordina le sinergie in base alla coerenza media beta (piu' facile leggere differenze)
betaMeanRaw = mean(C(:, betaBand), 2, 'omitnan');
[~, ord] = sort(betaMeanRaw, 'descend');
C = C(ord, :);
Cd = Cd(ord, :);
betaMeanRaw = betaMeanRaw(ord);
if isfield(metrics,'betaPeakCoh') && numel(metrics.betaPeakCoh)==nSyn
    betaPeakRaw = metrics.betaPeakCoh(ord);
else
    betaPeakRaw = max(C(:, betaBand), [], 2);
end
synLab = arrayfun(@(k) sprintf('S%d', k), ord(:)', 'UniformOutput', false);

% Delta rispetto alla sinergia di riferimento (quella con betaMean piu' alta)
refIdx = 1;
Cref = C - C(refIdx, :);

% Distanza tra profili beta delle sinergie (1 - correlazione)
Cbeta = C(:, betaBand);
D = zeros(nSyn,nSyn);
for i = 1:nSyn
    for j = 1:nSyn
        xi = Cbeta(i,:)'; xj = Cbeta(j,:)';
        ok = isfinite(xi) & isfinite(xj);
        if nnz(ok) < 3
            D(i,j) = NaN;
        else
            cc = corrcoef(xi(ok), xj(ok));
            D(i,j) = 1 - cc(1,2);
        end
    end
end

fig = figure('Visible','off','Color','w','Position',[100 100 1350 980]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(fb, 1:nSyn, C);
axis tight;
xlabel('Hz');
ylabel('Synergy');
yticks(1:nSyn);
yticklabels(synLab);
title('Coerenza assoluta per sinergia (ordinata per beta mean)');
cb1 = colorbar;
ylabel(cb1, 'Coherence');

ax2 = nexttile;
imagesc(fb, 1:nSyn, Cd);
axis tight;
xlabel('Hz');
ylabel('Synergy');
yticks(1:nSyn);
yticklabels(synLab);
title('Delta rispetto alla media delle sinergie');
cb2 = colorbar;
ylabel(cb2, '\Delta coherence');
cmax = max(abs(Cd(:)));
if isfinite(cmax) && cmax > 0
    caxis([-cmax cmax]);
end
colormap(ax2, redBlueMap());

ax3 = nexttile;
imagesc(fb, 1:nSyn, Cref);
axis tight;
xlabel('Hz');
ylabel('Synergy');
yticks(1:nSyn);
yticklabels(synLab);
title(sprintf('Delta rispetto alla sinergia di riferimento (%s)', synLab{refIdx}));
cb3 = colorbar;
ylabel(cb3, '\Delta coherence vs ref');
cmax2 = max(abs(Cref(:)));
if isfinite(cmax2) && cmax2 > 0
    caxis(ax3, [-cmax2 cmax2]);
end
colormap(ax3, redBlueMap());

nexttile;
hold on;
plot(1:nSyn, betaMeanRaw(:), '-o', 'LineWidth', 1.7, 'Color', [0.15 0.45 0.85], ...
    'MarkerFaceColor', [0.15 0.45 0.85], 'DisplayName', 'beta mean (13-30 Hz)');
plot(1:nSyn, betaPeakRaw(:), '-s', 'LineWidth', 1.7, 'Color', [0.85 0.30 0.20], ...
    'MarkerFaceColor', [0.85 0.30 0.20], 'DisplayName', 'beta peak');
if isfinite(metrics.c95)
    yline(metrics.c95, '--k', 'DisplayName', 'C95');
end
xlabel('Synergy');
ylabel('Coherence');
xticks(1:nSyn);
xticklabels(synLab);
title('Profilo beta per sinergia (ordinate)');
grid on;
xlim([0.5 nSyn+0.5]);
legend('Location','best');

saveas(fig, fullfile(outDir, 'coerenza_variazione_sinergie.png'));
close(fig);

% Salva anche una matrice compatta di distanza tra sinergie in beta
fig2 = figure('Visible','off','Color','w','Position',[100 100 700 620]);
imagesc(D);
axis square;
xticks(1:nSyn); yticks(1:nSyn);
xticklabels(synLab); yticklabels(synLab);
xlabel('Synergy'); ylabel('Synergy');
title('Distanza tra sinergie in beta (1-corr)');
cb4 = colorbar; ylabel(cb4, 'distance');
saveas(fig2, fullfile(outDir, 'coerenza_distanza_sinergie_beta.png'));
close(fig2);
end

function cmap = redBlueMap()
n = 256;
r = [(0:(n/2-1))/(n/2) ones(1,n/2)];
b = [ones(1,n/2) ((n/2-1):-1:0)/(n/2)];
g = [((0:(n/2-1))/(n/2)) ((n/2-1):-1:0)/(n/2)];
cmap = [r(:) g(:) b(:)];
end

function writeSynergyMetricsCsv(pathCsv, metrics)
header = {'synergy','cluster_id','peak_coherence','peak_hz','mean_1_45_hz', ...
    'beta_peak_coherence','beta_peak_hz','beta_mean_13_30_hz','beta_auc_13_30_hz', ...
    'zcoh_beta_peak','zcoh_beta_peak_hz','zcoh_beta_mean','zcoh_threshold_95','zcoh_sigfrac_beta', ...
    'zcoh_surrogate_p','zcoh_surrogate_q','zcoh_surrogate_sig', ...
    'beta_mean_lagopt_13_30_hz','lag_opt_ms','psi_beta', ...
    'm1_left_beta_13_30_hz','m1_right_beta_13_30_hz','m1_delta_l_minus_r_13_30_hz', ...
    'm1_delta_l_minus_r_early','m1_delta_l_minus_r_middle','m1_delta_l_minus_r_late', ...
    'm1_delta_l_minus_r_13_20_hz','m1_delta_l_minus_r_20_30_hz', ...
    'im_beta_mean_13_30_hz','phase_beta_mean_rad', ...
    'surrogate_p','surrogate_q','surrogate_sig', ...
    'beta_mean_ci95_low','beta_mean_ci95_high','beta_stability_std','beta_stability_cv', ...
    'beta_early','beta_middle','beta_late','beta_delta_late_minus_early', ...
    'lateralization_li_motR_minus_motL','c95'};
rows = {};
for k = 1:numel(metrics.peakCoh)
    cl = NaN;
    if isfield(metrics,'synergyClusterId') && numel(metrics.synergyClusterId) >= k
        cl = metrics.synergyClusterId(k);
    end
    rows(end+1,:) = {num2str(k), fmt6(cl), sprintf('%.6f',metrics.peakCoh(k)), ... %#ok<AGROW>
        sprintf('%.6f',metrics.peakHz(k)), sprintf('%.6f',metrics.meanBand(k)), ...
        sprintf('%.6f',metrics.betaPeakCoh(k)), sprintf('%.6f',metrics.betaPeakHz(k)), ...
        sprintf('%.6f',metrics.betaMean(k)), sprintf('%.6f',metrics.betaAuc(k)), ...
        fmt6(metrics.zCohBetaPeak(k)), fmt6(metrics.zCohBetaPeakHz(k)), fmt6(metrics.zCohBetaMean(k)), ...
        fmt6(metrics.zCohThreshold(k)), fmt6(metrics.zCohSigFracBeta(k)), ...
        fmt6(metrics.zCohSurrogateP(k)), fmt6(metrics.zCohSurrogateQ(k)), num2str(metrics.zCohSurrogateSig(k)), ...
        sprintf('%.6f',metrics.betaMeanLagOpt(k)), sprintf('%.6f',metrics.lagOptMs(k)), sprintf('%.6f',metrics.psiBeta(k)), ...
        fmt6(metrics.m1LeftBeta(k)), fmt6(metrics.m1RightBeta(k)), fmt6(metrics.m1DeltaLminusR(k)), ...
        fmt6(metrics.m1DeltaLminusREarly(k)), fmt6(metrics.m1DeltaLminusRMiddle(k)), fmt6(metrics.m1DeltaLminusRLate(k)), ...
        fmt6(metrics.m1DeltaLminusRLowBeta(k)), fmt6(metrics.m1DeltaLminusRHighBeta(k)), ...
        sprintf('%.6f',metrics.imBetaMean(k)), sprintf('%.6f',metrics.phaseBetaMeanRad(k)), ...
        fmt6(metrics.surrogateP(k)), fmt6(metrics.surrogateQ(k)), num2str(metrics.surrogateSig(k)), ...
        sprintf('%.6f',metrics.betaMeanCiLow(k)), sprintf('%.6f',metrics.betaMeanCiHigh(k)), ...
        sprintf('%.6f',metrics.betaStabilityStd(k)), sprintf('%.6f',metrics.betaStabilityCv(k)), ...
        sprintf('%.6f',metrics.betaEarly(k)), sprintf('%.6f',metrics.betaMiddle(k)), ...
        sprintf('%.6f',metrics.betaLate(k)), sprintf('%.6f',metrics.gestureDeltaLateEarly(k)), ...
        sprintf('%.6f',metrics.lateralizationLi(k)), ...
        sprintf('%.6f',metrics.c95)};
end
writeCsv(pathCsv, header, rows);
end

function writeZcoherenceMetricsCsv(pathCsv, metrics)
header = {'synergy','zcoh_beta_mean','zcoh_beta_peak','zcoh_beta_peak_hz','zcoh_threshold_95','zcoh_sigfrac_beta', ...
    'zcoh_surrogate_p','zcoh_surrogate_q','zcoh_surrogate_sig'};
rows = {};
if ~isfield(metrics,'zCohBetaMean')
    writeCsv(pathCsv, header, rows);
    return;
end
for k = 1:numel(metrics.zCohBetaMean)
    rows(end+1,:) = {num2str(k), fmt6(metrics.zCohBetaMean(k)), fmt6(metrics.zCohBetaPeak(k)), ... %#ok<AGROW>
        fmt6(metrics.zCohBetaPeakHz(k)), fmt6(metrics.zCohThreshold(k)), fmt6(metrics.zCohSigFracBeta(k)), ...
        fmt6(metrics.zCohSurrogateP(k)), fmt6(metrics.zCohSurrogateQ(k)), num2str(metrics.zCohSurrogateSig(k))};
end
writeCsv(pathCsv, header, rows);
end

function writeSynergyClusterCsv(pathCsv, metrics)
header = {'synergy','cluster_id','zcoh_beta_mean','beta_mean_13_30_hz'};
rows = {};
if ~isfield(metrics,'synergyClusterId') || isempty(metrics.synergyClusterId)
    writeCsv(pathCsv, header, rows);
    return;
end
n = numel(metrics.synergyClusterId);
for k = 1:n
    rows(end+1,:) = {num2str(k), fmt6(metrics.synergyClusterId(k)), ... %#ok<AGROW>
        fmt6(metrics.zCohBetaMean(k)), fmt6(metrics.betaMean(k))};
end
rows(end+1,:) = {'all', fmt6(metrics.synergyClusterK), fmt6(metrics.zCohBetaMeanAll), fmt6(metrics.betaMeanAll)}; %#ok<AGROW>
writeCsv(pathCsv, header, rows);
end

function makeZcoherencePlot(f, metrics, outDir)
if ~isfield(metrics,'zCohCurve') || isempty(metrics.zCohCurve) || isempty(f)
    return;
end
Z = metrics.zCohCurve;
nSyn = size(Z,1);
if nSyn < 1
    return;
end
beta = (f >= 13 & f <= 30);
zOrder = metrics.zCohBetaMean;
if numel(zOrder) ~= nSyn
    zOrder = mean(Z(:,beta),2,'omitnan');
end
[~, ord] = sort(zOrder, 'descend');

fig = figure('Visible','off','Color','w','Position',[120 120 1200 820]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(f, 1:nSyn, Z(ord,:));
axis tight;
xlabel('Hz');
ylabel('Synergy');
yticks(1:nSyn);
yticklabels(arrayfun(@(k) sprintf('S%d',k), ord(:)', 'UniformOutput', false));
title('CSC Z-coherence (atanh(sqrt(C)))');
cb = colorbar; ylabel(cb, 'Z-coherence');

nexttile;
hold on;
for i = 1:nSyn
    k = ord(i);
    plot(f, Z(k,:), 'LineWidth', 1.2, 'DisplayName', sprintf('S%d',k));
    if isfield(metrics,'zCohThreshold') && numel(metrics.zCohThreshold) >= k && isfinite(metrics.zCohThreshold(k))
        yline(metrics.zCohThreshold(k), '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility','off');
    end
end
xlim([min(f) max(f)]);
xlabel('Hz');
ylabel('Z-coherence');
title('Profili Z-coherence per sinergia + soglia 95% empirica');
grid on;
legend('Location','eastoutside');

saveas(fig, fullfile(outDir, 'csc_zcoherence_sinergie.png'));
close(fig);
end

function writeRoiMetricsCsv(pathCsv, metrics)
if ~isfield(metrics,'roiLabels') || isempty(metrics.roiLabels) || ...
   ~isfield(metrics,'roiBetaMeanAll') || isempty(metrics.roiBetaMeanAll)
    writeCsv(pathCsv, {'roi_name','rank','beta_mean_13_30_hz_all'}, {});
    return;
end

roiNames = cellstr(string(metrics.roiLabels(:)));
valsAll = metrics.roiBetaMeanAll(:);
nRoi = numel(roiNames);
nSyn = 0;
if isfield(metrics,'roiBetaMeanBySynergy') && ~isempty(metrics.roiBetaMeanBySynergy)
    nSyn = size(metrics.roiBetaMeanBySynergy,2);
end

[~, ord] = sort(valsAll, 'descend');
header = {'roi_name','rank','beta_mean_13_30_hz_all'};
for k = 1:nSyn
    header{end+1} = sprintf('beta_mean_syn_%02d', k); %#ok<AGROW>
end
rows = {};
for r = 1:nRoi
    i = ord(r);
    row = {roiNames{i}, num2str(r), fmt6(valsAll(i))};
    for k = 1:nSyn
        row{end+1} = fmt6(metrics.roiBetaMeanBySynergy(i,k)); %#ok<AGROW>
    end
    rows(end+1,:) = row; %#ok<AGROW>
end
writeCsv(pathCsv, header, rows);
end

function writeM1LateralizationCsv(pathCsv, metrics)
header = {'scope','m1_left_beta_13_30_hz','m1_right_beta_13_30_hz','m1_delta_l_minus_r_13_30_hz', ...
    'm1_delta_l_minus_r_early','m1_delta_l_minus_r_middle','m1_delta_l_minus_r_late', ...
    'm1_delta_l_minus_r_13_20_hz','m1_delta_l_minus_r_20_30_hz'};

if ~isfield(metrics,'m1LeftBetaAll')
    writeCsv(pathCsv, header, {});
    return;
end

rows = {};
rows(end+1,:) = {'all', fmt6(metrics.m1LeftBetaAll), fmt6(metrics.m1RightBetaAll), fmt6(metrics.m1DeltaLminusRAll), ...
    fmt6(metrics.m1DeltaLminusREarlyAll), fmt6(metrics.m1DeltaLminusRMiddleAll), fmt6(metrics.m1DeltaLminusRLateAll), ...
    fmt6(metrics.m1DeltaLminusRLowBetaAll), fmt6(metrics.m1DeltaLminusRHighBetaAll)}; %#ok<AGROW>

nSyn = 0;
if isfield(metrics,'m1LeftBeta') && ~isempty(metrics.m1LeftBeta)
    nSyn = numel(metrics.m1LeftBeta);
end
for k = 1:nSyn
    rows(end+1,:) = {sprintf('syn_%02d',k), fmt6(metrics.m1LeftBeta(k)), fmt6(metrics.m1RightBeta(k)), fmt6(metrics.m1DeltaLminusR(k)), ... %#ok<AGROW>
        fmt6(metrics.m1DeltaLminusREarly(k)), fmt6(metrics.m1DeltaLminusRMiddle(k)), fmt6(metrics.m1DeltaLminusRLate(k)), ...
        fmt6(metrics.m1DeltaLminusRLowBeta(k)), fmt6(metrics.m1DeltaLminusRHighBeta(k))};
end
writeCsv(pathCsv, header, rows);
end

function writeRoiStateDependentCsv(pathCsv, metrics)
baseHeader = {'roi_name','rank_by_delta_all','activity_threshold_beta_power', ...
    'n_trials_active','n_trials_inactive','delta_all_synergies','ratio_all_synergies'};
if ~isfield(metrics,'roiLabels') || ~isfield(metrics,'roiStateDelta') || isempty(metrics.roiLabels)
    writeCsv(pathCsv, baseHeader, {});
    return;
end

roiNames = cellstr(string(metrics.roiLabels(:)));
nRoi = numel(roiNames);
nSyn = size(metrics.roiStateDelta,2);

header = baseHeader;
for k = 1:nSyn
    header{end+1} = sprintf('csc_active_syn_%02d',k); %#ok<AGROW>
    header{end+1} = sprintf('csc_inactive_syn_%02d',k); %#ok<AGROW>
    header{end+1} = sprintf('csc_delta_syn_%02d',k); %#ok<AGROW>
    header{end+1} = sprintf('csc_ratio_syn_%02d',k); %#ok<AGROW>
end

[~, ord] = sort(metrics.roiStateDeltaAll, 'descend');
rows = {};
for r = 1:nRoi
    i = ord(r);
    row = {roiNames{i}, num2str(r), fmt6(metrics.roiActivityThreshold(i)), ...
        num2str(metrics.roiStateNActive(i)), num2str(metrics.roiStateNInactive(i)), ...
        fmt6(metrics.roiStateDeltaAll(i)), fmt6(metrics.roiStateRatioAll(i))};
    for k = 1:nSyn
        row{end+1} = fmt6(metrics.roiStateCscActive(i,k)); %#ok<AGROW>
        row{end+1} = fmt6(metrics.roiStateCscInactive(i,k)); %#ok<AGROW>
        row{end+1} = fmt6(metrics.roiStateDelta(i,k)); %#ok<AGROW>
        row{end+1} = fmt6(metrics.roiStateRatio(i,k)); %#ok<AGROW>
    end
    rows(end+1,:) = row; %#ok<AGROW>
end

writeCsv(pathCsv, header, rows);
end

function writePrePostBetaSummary(pathCsv, summaryRows, summaryHeader)
if isempty(summaryRows)
    writeCsv(pathCsv, {'key','pair_pre','pair_post','beta_mean_pre','beta_mean_post','beta_mean_delta_post_minus_pre', ...
        'beta_peak_pre','beta_peak_post','beta_peak_delta_post_minus_pre', ...
        'beta_auc_pre','beta_auc_post','beta_auc_delta_post_minus_pre'}, {});
    return;
end

idx = headerIndexMap(summaryHeader);
rows = summaryRows;
ok = false(size(rows,1),1);
for i = 1:size(rows,1)
    ok(i) = strcmpi(rows{i,idx('status')}, 'OK');
end
rows = rows(ok,:);
if isempty(rows)
    writeCsv(pathCsv, {'key','pair_pre','pair_post','beta_mean_pre','beta_mean_post','beta_mean_delta_post_minus_pre', ...
        'beta_peak_pre','beta_peak_post','beta_peak_delta_post_minus_pre', ...
        'beta_auc_pre','beta_auc_post','beta_auc_delta_post_minus_pre'}, {});
    return;
end

keys = cell(size(rows,1),1);
for i = 1:size(rows,1)
    keys{i} = prePostKeyFromPair(lower(char(rows{i,idx('pair')})));
end
u = unique(keys);
out = {};
for i = 1:numel(u)
    k = u{i};
    ip = find(strcmp(keys,k) & strcmpi(cellstr(string(rows(:,idx('phase')))), 'pre'), 1, 'first');
    io = find(strcmp(keys,k) & strcmpi(cellstr(string(rows(:,idx('phase')))), 'post'), 1, 'first');
    if isempty(ip) || isempty(io)
        continue;
    end
    bmPre = toNum(rows{ip,idx('beta_mean_all')}); bmPost = toNum(rows{io,idx('beta_mean_all')});
    bpPre = toNum(rows{ip,idx('beta_peak_all')}); bpPost = toNum(rows{io,idx('beta_peak_all')});
    baPre = toNum(rows{ip,idx('beta_auc_all')});  baPost = toNum(rows{io,idx('beta_auc_all')});
    out(end+1,:) = {k, rows{ip,idx('pair')}, rows{io,idx('pair')}, ... %#ok<AGROW>
        fmt6(bmPre), fmt6(bmPost), fmt6(bmPost-bmPre), ...
        fmt6(bpPre), fmt6(bpPost), fmt6(bpPost-bpPre), ...
        fmt6(baPre), fmt6(baPost), fmt6(baPost-baPre)};
end

hdr = {'key','pair_pre','pair_post','beta_mean_pre','beta_mean_post','beta_mean_delta_post_minus_pre', ...
    'beta_peak_pre','beta_peak_post','beta_peak_delta_post_minus_pre', ...
    'beta_auc_pre','beta_auc_post','beta_auc_delta_post_minus_pre'};
writeCsv(pathCsv, hdr, out);
end

function m = headerIndexMap(h)
m = containers.Map();
for i = 1:numel(h)
    m(h{i}) = i;
end
end

function k = prePostKeyFromPair(s)
s = regexprep(s, 'ses[-_]?pre', 'ses');
s = regexprep(s, 'ses[-_]?post', 'ses');
s = regexprep(s, '(^|[_-])pre([_-]|$)', '$1$2');
s = regexprep(s, '(^|[_-])post([_-]|$)', '$1$2');
s = regexprep(s, '[_-]+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'pair'; end
k = s;
end

function v = toNum(x)
v = str2double(string(x));
if isnan(v), v = NaN; end
end

function s = fmt6(v)
if isnan(v), s = 'NaN'; else, s = sprintf('%.6f',v); end
end

function v = meanNoNan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = mean(x); end
end

function v = maxNoNan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end

function v = minNoNan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = min(x); end
end

function writeTrialInfoCsv(pathCsv, trialInfo, fs)
if nargin < 3 || isempty(fs) || ~isfinite(fs) || fs <= 0
    fs = NaN;
end

durSec = nan(size(trialInfo,1),1);
for i = 1:size(trialInfo,1)
    nS = trialInfo(i,3) - trialInfo(i,2) + 1;
    if isfinite(fs)
        durSec(i) = nS / fs;
    end
end

durFinite = durSec(isfinite(durSec));
if isempty(durFinite)
    runMedianSec = NaN;
    runQ1Sec = NaN;
    runQ3Sec = NaN;
else
    runMedianSec = median(durFinite);
    runQ1Sec = prctile(durFinite, 25);
    runQ3Sec = prctile(durFinite, 75);
end

header = {'trial_idx','start_sample','end_sample','n_samples','duration_sec', ...
    'duration_class', ...
    'run_q1_duration_sec','run_median_duration_sec','run_q3_duration_sec', ...
    'flag_lt_0p75s','flag_gt_1p25s','flag_below_run_q1','flag_above_run_q3'};
rows = {};
for i = 1:size(trialInfo,1)
    nS = trialInfo(i,3) - trialInfo(i,2) + 1;
    dS = durSec(i);
    if ~isfinite(dS)
        dSstr = '';
        durClass = 'unknown';
        isShort = '';
        isLong = '';
        isBelowQ1 = '';
        isAboveQ3 = '';
    else
        dSstr = sprintf('%.6f', dS);
        durClass = durationClassLabel(dS);
        isShort = logicalToCsv(dS < 0.75);
        isLong = logicalToCsv(dS > 1.25);
        isBelowQ1 = logicalToCsv(isfinite(runQ1Sec) && dS < runQ1Sec);
        isAboveQ3 = logicalToCsv(isfinite(runQ3Sec) && dS > runQ3Sec);
    end
    rows(end+1,:) = { ...
        num2str(trialInfo(i,1)), ...
        num2str(trialInfo(i,2)), ...
        num2str(trialInfo(i,3)), ...
        num2str(nS), ...
        dSstr, ...
        durClass, ...
        fmtNumOrBlank(runQ1Sec, 6), ...
        fmtNumOrBlank(runMedianSec, 6), ...
        fmtNumOrBlank(runQ3Sec, 6), ...
        isShort, ...
        isLong, ...
        isBelowQ1, ...
        isAboveQ3}; %#ok<AGROW>
end
writeCsv(pathCsv, header, rows);
end

function s = durationClassLabel(dSec)
if ~isfinite(dSec)
    s = 'unknown';
elseif dSec < 0.75
    s = 'short';
elseif dSec <= 1.25
    s = 'mid';
else
    s = 'long';
end
end

function s = logicalToCsv(tf)
if tf
    s = '1';
else
    s = '0';
end
end

function s = fmtNumOrBlank(x, nDec)
if nargin < 2 || isempty(nDec)
    nDec = 6;
end
if isfinite(x)
    s = sprintf(['%0.' num2str(nDec) 'f'], x);
else
    s = '';
end
end

function writeCsv(pathCsv, header, rows)
fid = fopen(pathCsv,'w');
fprintf(fid,'%s\n',strjoin(header,','));
for i = 1:size(rows,1)
    r = rows(i,:);
    for j = 1:numel(r), r{j} = esc(r{j}); end
    fprintf(fid,'%s\n',strjoin(r,','));
end
fclose(fid);
end

function s = esc(v)
if ~ischar(v), v = char(string(v)); end
if contains(v,'"'), v = strrep(v,'"','""'); end
if contains(v,',') || contains(v,'"') || contains(v,newline)
    s = ['"' v '"'];
else
    s = v;
end
end

function s = safeNameFromPath(p)
[~, s, ~] = fileparts(fileparts(p));
if isempty(s), s = 'pair'; end
end

function [pairName, phase] = safePairPhaseFromReportPath(p)
phaseDir = fileparts(p);
[pairDir, phase] = fileparts(phaseDir);
[~, pairName] = fileparts(pairDir);
if isempty(pairName), pairName = safeNameFromPath(p); end
if isempty(phase), phase = 'session'; end
end

function v = getFieldDef(S, field, def)
if isfield(S, field)
    v = S.(field);
else
    v = def;
end
end

function t = oneToken(txt, pat)
q = regexp(txt, pat, 'tokens', 'once');
if isempty(q), error('Pattern not found: %s', pat); end
t = q{1};
end

function epochSec = isoToEpoch(isoStr)
d = datetime(isoStr, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX', 'TimeZone','UTC');
epochSec = posixtime(d);
end

function p = prct(x, q, dim)
if nargin < 3, dim = 1; end
try
    p = prctile(x, q, dim);
catch
    if dim ~= 2
        x = x';
    end
    p = zeros(size(x,1),1);
    for i = 1:size(x,1)
        v = sort(x(i,:));
        if isempty(v)
            p(i) = NaN;
            continue;
        end
        pos = (q/100)*(numel(v)-1) + 1;
        lo = floor(pos); hi = ceil(pos);
        if lo == hi
            p(i) = v(lo);
        else
            p(i) = v(lo) + (pos-lo)*(v(hi)-v(lo));
        end
    end
    if dim ~= 2
        p = p';
    end
end
end

function cleanupTemp(d)
if isfolder(d)
    try
        rmdir(d, 's');
    catch
    end
end
end
