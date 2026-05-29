function validazione_eeg_emg_sync(inputRoot, alignedRoot, validationRoot, overwriteOutput, validationMode)
%VALIDAZIONE_EEG_EMG_SYNC
% Validate alignment quality and CSC-readiness with plots and criteria report.
%
% Supported input layouts:
% 1) Legacy pair structure
%   pair_xxx/pre/eeg + pair_xxx/pre/emg
%   pair_xxx/post/eeg + pair_xxx/post/emg
% 2) BIDS-like structure
%   sub-XX/ses-YY/...  (one .mff + one .otb+ per session)
%
% If one phase is empty/incomplete, phase is skipped.
%
% Validation logic note (aligned with sync alignment script):
% - trial gating relies on DIN4/DIN5 temporal order only;
% - DIN1 is never used as a trial-validity or CSC gating condition;
% - when a dense DIN1 train is available on both EEG and EMG, it is used
%   only to refine the affine fit after the DIN4/DIN5-based fit.

if nargin < 1 || isempty(inputRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    inputRoot = fullfile(scriptDir, 'batch_otb_mff', 'input');
end
if nargin < 2 || isempty(alignedRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    alignedRoot = fullfile(scriptDir, 'batch_otb_mff', 'output_matlab_sync');
end
if nargin < 3 || isempty(validationRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    validationRoot = fullfile(scriptDir, 'batch_otb_mff', 'validazione_sync');
end
if nargin < 4 || isempty(overwriteOutput)
    overwriteOutput = true;
end
if nargin < 5 || isempty(validationMode)
    validationMode = 'last_only'; % 'last_only' | 'full'
end

if ~isfolder(inputRoot), error('Input root not found: %s', inputRoot); end
if ~isfolder(alignedRoot), error('Aligned root not found: %s', alignedRoot); end
if isfolder(validationRoot) && overwriteOutput
    rmdir(validationRoot, 's');
end
if ~isfolder(validationRoot), mkdir(validationRoot); end

tolSec = 0.020; % 20 ms
movementOnly = true;
cscProfile = 'very_strict'; % 'standard' | 'strict' | 'very_strict'
criteria = getCriteria(cscProfile);

entries = collectInputEntries(inputRoot);
if isempty(entries), error('No valid input entries in %s', inputRoot); end
if movementOnly
    mask = false(1, numel(entries));
    for i = 1:numel(entries)
        mask(i) = isMovementEntry(entries(i));
    end
    entries = entries(mask);
    fprintf('Movement-only mode: %d entries found.\n', numel(entries));
    if isempty(entries)
        error('No movement entries found in %s', inputRoot);
    end
end

header = {'pair','phase','status','message','csc_ready','reason', ...
    'a','b','drift_ppm','mad_fit_ms','din1_recall','din4_recall','din5_recall', ...
    'din1_med_abs_ms','din4_med_abs_ms','din5_med_abs_ms', ...
    'din1_p95_abs_ms','din4_p95_abs_ms','din5_p95_abs_ms', ...
    'fit_stage','refit_on_valid_trials','refit_with_din1','n_din1_match','n_din1_used','din1_tol_ms','din1_match_ratio', ...
    'blocking_high','blocking_total', ...
    'csc_scarto_ms','csc_soglia_ms','p95_mov_max_ms','margine_csc_ms', ...
    'correzione_metodo','correzione_formula','correzione_mapping', ...
    'soglia_profilo','soglia_dettaglio', ...
    'report_json'};
rows = {};
idx = @(name) find(strcmp(header,name), 1);
sogliaProfilo = upper(cscProfile);
sogliaDettaglio = sprintf(['recall(DIN4/5)>=%.2f/%.2f; med(DIN4/5)<=%.1f/%.1f ms; ' ...
    'p95(DIN4/5)<=%.1f/%.1f ms; MAD<=%.1f ms'], ...
    criteria.recallDin4, criteria.recallDin5, ...
    criteria.medAbsDin4Ms, criteria.medAbsDin5Ms, ...
    criteria.p95Din4Ms, criteria.p95Din5Ms, ...
    criteria.maxMadFitMs);
criteriaHeader = {'pair','phase','file_tag','criterion','threshold','value','pass','note'};
criteriaRows = {};
plotHeader = {'pair','phase','file_tag','plot_type','plot_path'};
plotRows = {};

for i = 1:numel(entries)
    pairName = entries(i).pairName;
    phase = entries(i).phase;
    otbPath = entries(i).otbPath;
    state = entries(i).state;
    msg = entries(i).msg;
    try
        if ~strcmp(state,'OK')
            fprintf('[SKIP] %s/%s | %s\n', pairName, phase, msg);
            row = repmat({''}, 1, numel(header));
            row{idx('pair')} = pairName;
            row{idx('phase')} = phase;
            row{idx('status')} = 'SKIP';
            row{idx('message')} = msg;
            row{idx('soglia_profilo')} = sogliaProfilo;
            row{idx('soglia_dettaglio')} = sogliaDettaglio;
            rows(end+1,:) = row; %#ok<AGROW>
            continue;
        end

        phaseAlignedPath = fullfile(alignedRoot, pairName, phase);
        alignedMff = findOne(phaseAlignedPath, '*.mff', true);

        pairValOut = validationRoot;
        fileTag = makeFileTag(pairName, phase);

        ref = readMffReference(alignedMff, 'Events_8 DINs.xml');
        otb = readOtb(otbPath);
        [pred, meta] = buildAlignedEvents(ref, otb);

        metrics = computeMetrics(ref, pred, meta, tolSec);
        blockers = computeBlockingDifferences(ref, pred, tolSec);
        [isReady, reason, checkRows] = evaluateCscReadiness(metrics, meta, blockers, criteria, validationMode);
        sc = computeCscScarto(metrics, criteria);

        makePlots(pairValOut, fileTag, pairName, phase, ref, pred, metrics, tolSec, validationMode);

        for rr = 1:size(checkRows,1)
            criteriaRows(end+1,:) = {pairName, phase, fileTag, ... %#ok<AGROW>
                checkRows{rr,1}, checkRows{rr,2}, checkRows{rr,3}, checkRows{rr,4}, checkRows{rr,5}};
        end
        if strcmpi(validationMode,'last_only')
            finalPlot = fullfile(pairValOut,sprintf('%s_01_final_validation.png', fileTag));
            plotRows(end+1,:) = {pairName, phase, fileTag, 'final_dashboard', finalPlot}; %#ok<AGROW>
        else
            p1 = fullfile(pairValOut,sprintf('%s_01_overlay_events.png', fileTag));
            p2 = fullfile(pairValOut,sprintf('%s_02_residuals_vs_time.png', fileTag));
            p3 = fullfile(pairValOut,sprintf('%s_03_residual_hist.png', fileTag));
            plotRows(end+1,:) = {pairName, phase, fileTag, 'overlay', p1}; %#ok<AGROW>
            plotRows(end+1,:) = {pairName, phase, fileTag, 'residuals_time', p2}; %#ok<AGROW>
            plotRows(end+1,:) = {pairName, phase, fileTag, 'residuals_hist', p3}; %#ok<AGROW>
        end

        fprintf('[OK] %s/%s | CSC_READY=%d | reason=%s\n', pairName, phase, isReady, reason);

        row = repmat({''}, 1, numel(header));
        row{idx('pair')} = pairName;
        row{idx('phase')} = phase;
        row{idx('status')} = 'OK';
        row{idx('message')} = '';
        row{idx('csc_ready')} = tf(isReady);
        row{idx('reason')} = reason;
        row{idx('a')} = sprintf('%.12f',meta.a);
        row{idx('b')} = sprintf('%.12f',meta.b);
        row{idx('drift_ppm')} = sprintf('%.2f',meta.driftPpm);
        row{idx('mad_fit_ms')} = sprintf('%.6f',meta.mad*1000);
        row{idx('din1_recall')} = sprintf('%.4f',metrics.match.DIN1.recall);
        row{idx('din4_recall')} = sprintf('%.4f',metrics.match.DIN4.recall);
        row{idx('din5_recall')} = sprintf('%.4f',metrics.match.DIN5.recall);
        row{idx('din1_med_abs_ms')} = sprintf('%.3f',metrics.match.DIN1.medAbsMs);
        row{idx('din4_med_abs_ms')} = sprintf('%.3f',metrics.match.DIN4.medAbsMs);
        row{idx('din5_med_abs_ms')} = sprintf('%.3f',metrics.match.DIN5.medAbsMs);
        row{idx('din1_p95_abs_ms')} = sprintf('%.3f',metrics.match.DIN1.p95AbsMs);
        row{idx('din4_p95_abs_ms')} = sprintf('%.3f',metrics.match.DIN4.p95AbsMs);
        row{idx('din5_p95_abs_ms')} = sprintf('%.3f',metrics.match.DIN5.p95AbsMs);
        row{idx('fit_stage')} = meta.fit.stage;
        row{idx('refit_on_valid_trials')} = tf(meta.fit.refit_on_valid_trials);
        row{idx('refit_with_din1')} = tf(meta.fit.refit_with_din1);
        row{idx('n_din1_match')} = num2str(meta.fit.n_din1_match);
        row{idx('n_din1_used')} = num2str(meta.fit.n_din1_used);
        row{idx('din1_tol_ms')} = sprintf('%.3f',meta.fit.din1_tol_ms);
        row{idx('din1_match_ratio')} = sprintf('%.4f',meta.fit.din1_match_ratio);
        row{idx('blocking_high')} = num2str(blockers.highCount);
        row{idx('blocking_total')} = num2str(blockers.totalCount);
        row{idx('csc_scarto_ms')} = sprintf('%.3f',sc.cscScartoMs);
        row{idx('csc_soglia_ms')} = sprintf('%.3f',sc.sogliaMs);
        row{idx('p95_mov_max_ms')} = sprintf('%.3f',sc.p95MovMaxMs);
        row{idx('margine_csc_ms')} = sprintf('%.3f',sc.margineMs);
        if meta.fit.refit_with_din1
            row{idx('correzione_metodo')} = 'Affine fit on DIN4/DIN5 (minimum MAD) + optional DIN1 sync refinement';
        else
            row{idx('correzione_metodo')} = 'Affine fit on DIN4/DIN5 (minimum MAD)';
        end
        row{idx('correzione_formula')} = sprintf('t_corr = %.9f * t_otb + %.9f', meta.a, meta.b);
        row{idx('correzione_mapping')} = sprintf('AUX4->%s; AUX3->%s; AUX2->DIN1', meta.mapping.AUX4, meta.mapping.AUX3);
        row{idx('soglia_profilo')} = sogliaProfilo;
        row{idx('soglia_dettaglio')} = sogliaDettaglio;
        row{idx('report_json')} = '';
        rows(end+1,:) = row; %#ok<AGROW>

    catch ME
        fprintf('[ERR] %s/%s | %s\n', pairName, phase, ME.message);
        row = repmat({''}, 1, numel(header));
        row{idx('pair')} = pairName;
        row{idx('phase')} = phase;
        row{idx('status')} = 'ERR';
        row{idx('message')} = ME.message;
        row{idx('soglia_profilo')} = sogliaProfilo;
        row{idx('soglia_dettaglio')} = sogliaDettaglio;
        rows(end+1,:) = row; %#ok<AGROW>
    end
end

summaryXlsx = fullfile(validationRoot, 'summary_validazione_movimento.xlsx');
if isempty(rows)
    Tsum = cell2table(cell(0, numel(header)), 'VariableNames', header);
else
    Tsum = cell2table(rows, 'VariableNames', header);
end
writetable(Tsum, summaryXlsx, 'FileType', 'spreadsheet', 'Sheet', 'Summary');

Tdiff = buildPrePostDiffTable(Tsum);
writetable(Tdiff, summaryXlsx, 'FileType', 'spreadsheet', 'Sheet', 'PrePost_Diff');

if isempty(criteriaRows)
    Tcrit = cell2table(cell(0, numel(criteriaHeader)), 'VariableNames', criteriaHeader);
else
    Tcrit = cell2table(criteriaRows, 'VariableNames', criteriaHeader);
end
writetable(Tcrit, summaryXlsx, 'FileType', 'spreadsheet', 'Sheet', 'Criteri');

if isempty(plotRows)
    Tplot = cell2table(cell(0, numel(plotHeader)), 'VariableNames', plotHeader);
else
    Tplot = cell2table(plotRows, 'VariableNames', plotHeader);
end
writetable(Tplot, summaryXlsx, 'FileType', 'spreadsheet', 'Sheet', 'Grafici');
fprintf('\nValidation summary (Excel): %s\n', summaryXlsx);
end

function criteria = getCriteria(profile)
if nargin < 1 || isempty(profile)
    profile = 'standard';
end

criteria = struct();
switch lower(profile)
    case 'very_strict'
        criteria.minTrialsDin4 = 20;
        criteria.minTrialsDin5 = 20;
        criteria.recallDin4 = 0.98;
        criteria.recallDin5 = 0.98;
        criteria.medAbsDin4Ms = 3;
        criteria.medAbsDin5Ms = 3;
        criteria.p95Din4Ms = 8;
        criteria.p95Din5Ms = 8;
        criteria.maxMadFitMs = 1;
        criteria.maxHighBlockers = 0;
        criteria.driftWarningPpm = 2000;
    case 'strict'
        criteria.minTrialsDin4 = 20;
        criteria.minTrialsDin5 = 20;
        criteria.recallDin4 = 0.97;
        criteria.recallDin5 = 0.97;
        criteria.medAbsDin4Ms = 4;
        criteria.medAbsDin5Ms = 4;
        criteria.p95Din4Ms = 10;
        criteria.p95Din5Ms = 10;
        criteria.maxMadFitMs = 1.5;
        criteria.maxHighBlockers = 0;
        criteria.driftWarningPpm = 3000;
    otherwise % standard
        criteria.minTrialsDin4 = 20;
        criteria.minTrialsDin5 = 20;
        criteria.recallDin4 = 0.95;
        criteria.recallDin5 = 0.95;
        criteria.medAbsDin4Ms = 5;
        criteria.medAbsDin5Ms = 5;
        criteria.p95Din4Ms = 15;
        criteria.p95Din5Ms = 15;
        criteria.maxMadFitMs = 2;
        criteria.maxHighBlockers = 0;
        criteria.driftWarningPpm = 5000;
end
end

function sc = computeCscScarto(metrics, criteria)
p95MovMax = max([metrics.match.DIN4.p95AbsMs, metrics.match.DIN5.p95AbsMs]);
soglia = max([criteria.p95Din4Ms, criteria.p95Din5Ms]);
margine = soglia - p95MovMax;

% CSC scarto is intentionally based on movement markers (DIN4/DIN5) and fit quality.
excess = [ ...
    metrics.match.DIN4.medAbsMs - criteria.medAbsDin4Ms; ...
    metrics.match.DIN5.medAbsMs - criteria.medAbsDin5Ms; ...
    metrics.match.DIN4.p95AbsMs - criteria.p95Din4Ms; ...
    metrics.match.DIN5.p95AbsMs - criteria.p95Din5Ms; ...
    metrics.madFitMs - criteria.maxMadFitMs];
excess = excess(~isnan(excess));
if isempty(excess)
    cscScarto = NaN;
else
    cscScarto = max([0; excess]);
end

sc = struct('cscScartoMs', cscScarto, ...
            'sogliaMs', soglia, ...
            'p95MovMaxMs', p95MovMax, ...
            'margineMs', margine);
end

function [isReady, reason, checkRows] = evaluateCscReadiness(metrics, meta, blockers, criteria, validationMode)
c1 = metrics.match.DIN1; c4 = metrics.match.DIN4; c5 = metrics.match.DIN5;

checks = {};

% Gating follows the same trial logic as alignment: DIN4/DIN5 drive validity.
% DIN1 remains diagnostic-only and does not block CSC readiness.
okCounts = metrics.countsMff.DIN4 >= criteria.minTrialsDin4 && metrics.countsMff.DIN5 >= criteria.minTrialsDin5;
okRecall = c4.recall >= criteria.recallDin4 && c5.recall >= criteria.recallDin5;
okResidual = c4.medAbsMs <= criteria.medAbsDin4Ms && c5.medAbsMs <= criteria.medAbsDin5Ms && ...
             c4.p95AbsMs <= criteria.p95Din4Ms && c5.p95AbsMs <= criteria.p95Din5Ms;
okFit = metrics.madFitMs <= criteria.maxMadFitMs;
okBlock = blockers.highCount <= criteria.maxHighBlockers;

highDrift = abs(meta.driftPpm) > criteria.driftWarningPpm;
isReady = okCounts && okRecall && okResidual && okFit && okBlock;
if isReady
    if highDrift
        reason = 'alignment_good_for_csc_high_drift_corrected';
    else
        reason = 'alignment_good_for_csc';
    end
else
    why = {};
    if ~okCounts, why{end+1}='insufficient_trials'; end %#ok<AGROW>
    if ~okRecall, why{end+1}='low_recall'; end %#ok<AGROW>
    if ~okResidual, why{end+1}='high_residuals'; end %#ok<AGROW>
    if ~okFit, why{end+1}='high_fit_error'; end %#ok<AGROW>
    if ~okBlock, why{end+1}='blocking_differences'; end %#ok<AGROW>
    if highDrift, why{end+1}='high_clock_drift_warning'; end %#ok<AGROW>
    reason = strjoin(why,';');
end

if strcmpi(validationMode,'last_only')
    checks(end+1,:) = {'final_csc_ready','must be TRUE',tf(isReady),tf(isReady),reason};
    if highDrift
        checks(end+1,:) = mkCheck('drift_warning_ppm', sprintf('warning if > %.0f',criteria.driftWarningPpm), meta.driftPpm, true, ...
            'HIGH DRIFT WARNING (corrected in alignment)'); %#ok<AGROW>
    end
else
    checks(end+1,:) = mkCheck('min_trials_DIN4', sprintf('>= %d',criteria.minTrialsDin4), metrics.countsMff.DIN4, metrics.countsMff.DIN4 >= criteria.minTrialsDin4, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('min_trials_DIN5', sprintf('>= %d',criteria.minTrialsDin5), metrics.countsMff.DIN5, metrics.countsMff.DIN5 >= criteria.minTrialsDin5, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('recall_DIN4', sprintf('>= %.2f',criteria.recallDin4), c4.recall, c4.recall >= criteria.recallDin4, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('recall_DIN5', sprintf('>= %.2f',criteria.recallDin5), c5.recall, c5.recall >= criteria.recallDin5, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('recall_DIN1_diag', 'diagnostic only', c1.recall, true, 'not used for gating'); %#ok<AGROW>
    checks(end+1,:) = mkCheck('med_abs_DIN4_ms', sprintf('<= %.1f',criteria.medAbsDin4Ms), c4.medAbsMs, c4.medAbsMs <= criteria.medAbsDin4Ms, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('med_abs_DIN5_ms', sprintf('<= %.1f',criteria.medAbsDin5Ms), c5.medAbsMs, c5.medAbsMs <= criteria.medAbsDin5Ms, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('p95_abs_DIN4_ms', sprintf('<= %.1f',criteria.p95Din4Ms), c4.p95AbsMs, c4.p95AbsMs <= criteria.p95Din4Ms, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('p95_abs_DIN5_ms', sprintf('<= %.1f',criteria.p95Din5Ms), c5.p95AbsMs, c5.p95AbsMs <= criteria.p95Din5Ms, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('p95_abs_DIN1_ms_diag', 'diagnostic only', c1.p95AbsMs, true, 'not used for gating'); %#ok<AGROW>
    checks(end+1,:) = mkCheck('mad_fit_ms', sprintf('<= %.1f',criteria.maxMadFitMs), metrics.madFitMs, metrics.madFitMs <= criteria.maxMadFitMs, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('high_blockers', sprintf('<= %d',criteria.maxHighBlockers), blockers.highCount, blockers.highCount <= criteria.maxHighBlockers, ''); %#ok<AGROW>
    checks(end+1,:) = mkCheck('drift_warning_ppm', sprintf('warning if > %.0f',criteria.driftWarningPpm), meta.driftPpm, true, ...
        ternary(highDrift,'HIGH DRIFT WARNING (corrected in alignment)','')); %#ok<AGROW>
end

checkRows = checks;
end

function row = mkCheck(name, thr, value, pass, note)
if isnan(value)
    v = 'NaN';
elseif isfloat(value)
    v = sprintf('%.6f', value);
else
    v = num2str(value);
end
row = {name, thr, v, tf(pass), note};
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function blockers = computeBlockingDifferences(ref, pred, tolSec)
codes = {'DIN1','DIN4','DIN5'};
header = {'code','issue','severity','time_mff_s','time_pred_s','delta_ms','note'};
rows = {};
highCount = 0;

for i = 1:numel(codes)
    code = codes{i};
    tMff = ref.codeTimes.(code);
    tPred = pred.times(strcmp(pred.codes, code));
    [resMs,~,~,~,idxPred,idxMff] = pairResidualsMs(tPred,tMff,tolSec);

    miss = setdiff(1:numel(tMff), idxMff);
    extra = setdiff(1:numel(tPred), idxPred);

    for k = miss
        sev = sevRule(code,'missing');
        if strcmp(sev,'HIGH'), highCount = highCount + 1; end
        rows(end+1,:) = {code,'missing_event',sev,sprintf('%.6f',tMff(k)),'','','Event in MFF not matched'}; %#ok<AGROW>
    end
    for k = extra
        sev = sevRule(code,'extra');
        if strcmp(sev,'HIGH'), highCount = highCount + 1; end
        rows(end+1,:) = {code,'extra_event',sev,'',sprintf('%.6f',tPred(k)),'','Event in OTB not matched'}; %#ok<AGROW>
    end

    if strcmp(code,'DIN1'), hard = 30; else, hard = 20; end
    for k = 1:numel(resMs)
        if abs(resMs(k)) > hard
            sev = sevRule(code,'residual');
            if strcmp(sev,'HIGH'), highCount = highCount + 1; end
            rows(end+1,:) = {code,'large_residual',sev,sprintf('%.6f',tMff(idxMff(k))),sprintf('%.6f',tPred(idxPred(k))),sprintf('%.3f',resMs(k)),sprintf('|delta|>%d ms',hard)}; %#ok<AGROW>
        end
    end
end

blockers = struct('header',{header},'rows',{rows},'highCount',highCount,'totalCount',size(rows,1));
end

function sev = sevRule(code, issue)
if strcmp(issue,'residual') || strcmp(issue,'missing')
    if strcmp(code,'DIN1'), sev='MED'; else, sev='HIGH'; end
    return;
end
if strcmp(code,'DIN1'), sev='LOW'; else, sev='MED'; end
end

function metrics = computeMetrics(ref, pred, meta, tolSec)
codes = {'DIN1','DIN4','DIN5'};
countsMff = struct(); countsPred = struct(); match = struct();
for i = 1:numel(codes)
    code = codes{i};
    tMff = ref.codeTimes.(code);
    tPred = pred.times(strcmp(pred.codes, code));
    [resMs,nMatch,rec,prec] = pairResidualsMs(tPred,tMff,tolSec);
    a = abs(resMs);
    m = struct('nMatch',nMatch,'recall',rec,'precision',prec,'residualsMs',resMs);
    if isempty(a)
        m.medAbsMs = NaN; m.p95AbsMs = NaN; m.maxAbsMs = NaN;
    else
        m.medAbsMs = median(a); m.p95AbsMs = prct(a,95); m.maxAbsMs = max(a);
    end
    match.(code) = m;
    countsMff.(code) = numel(tMff);
    countsPred.(code) = numel(tPred);
end
metrics = struct('match',match,'countsMff',countsMff,'countsPred',countsPred,'madFitMs',meta.mad*1000,'driftPpm',meta.driftPpm);
end

function makePlots(outDir, fileTag, pairName, phase, ref, pred, metrics, tolSec, validationMode)
codes = {'DIN1','DIN4','DIN5'};
clr = [0.1 0.4 0.9; 0.1 0.7 0.2; 0.9 0.3 0.1];

if strcmpi(validationMode,'last_only')
    f = figure('Visible','off','Color','w','Position',[100 100 1300 760]);
    for i = 1:3
        code = codes{i};
        subplot(2,2,i);
        tMff = ref.codeTimes.(code);
        tPred = pred.times(strcmp(pred.codes,code));
        scatter(tMff, ones(size(tMff)), 12, clr(i,:), 'filled'); hold on;
        scatter(tPred, 2*ones(size(tPred)), 12, clr(i,:), 'o');
        yline(1,'-','MFF','Color',[0.4 0.4 0.4]); yline(2,'-','OTB->MFF','Color',[0.4 0.4 0.4]);
        ylim([0.5 2.5]); xlim([0 ref.durationS]); grid on;
        title(sprintf('%s | rec=%.3f | p95=%.2f ms', code, metrics.match.(code).recall, metrics.match.(code).p95AbsMs));
        xlabel('Time (s)'); ylabel('Events');
    end

    subplot(2,2,4);
    vals = [metrics.match.DIN1.medAbsMs, metrics.match.DIN4.medAbsMs, metrics.match.DIN5.medAbsMs; ...
            metrics.match.DIN1.p95AbsMs, metrics.match.DIN4.p95AbsMs, metrics.match.DIN5.p95AbsMs]';
    b = bar(vals); %#ok<NASGU>
    set(gca,'XTickLabel',codes);
    ylabel('Residual (ms)');
    legend({'Median |dt|','P95 |dt|'}, 'Location','northwest');
    grid on;
    title(sprintf('Final validation | tol %.1f ms', tolSec*1000));

    sgtitle(sprintf('%s/%s | Final alignment validation',pairName,phase),'Interpreter','none');
    saveas(f, fullfile(outDir,sprintf('%s_01_final_validation.png', fileTag)));
    close(f);
    return;
end

f1 = figure('Visible','off','Color','w','Position',[100 100 1200 700]);
for i = 1:3
    code = codes{i};
    subplot(3,1,i);
    tMff = ref.codeTimes.(code);
    tPred = pred.times(strcmp(pred.codes,code));
    scatter(tMff, ones(size(tMff)), 12, clr(i,:), 'filled'); hold on;
    scatter(tPred, 2*ones(size(tPred)), 12, clr(i,:), 'o');
    yline(1,'-','MFF','Color',[0.4 0.4 0.4]); yline(2,'-','OTB->MFF','Color',[0.4 0.4 0.4]);
    ylim([0.5 2.5]); xlim([0 ref.durationS]); grid on;
    title(sprintf('%s/%s | %s',pairName,phase,code),'Interpreter','none');
    xlabel('Time (s)'); ylabel('Events');
end
sgtitle(sprintf('Event Overlay | tol %.1f ms',tolSec*1000));
saveas(f1, fullfile(outDir,sprintf('%s_01_overlay_events.png', fileTag))); close(f1);

f2 = figure('Visible','off','Color','w','Position',[120 120 1200 700]);
for i = 1:3
    code = codes{i}; subplot(3,1,i);
    [tm, rs] = residualsVsTime(pred, ref, code, tolSec);
    if ~isempty(tm)
        scatter(tm, rs, 12, clr(i,:), 'filled'); hold on;
        yline(0,'k-'); yline(5,'r--'); yline(-5,'r--');
    end
    grid on; xlabel('MFF time (s)'); ylabel('\Delta t (ms)');
    title(sprintf('%s | med=%.2f p95=%.2f ms',code,metrics.match.(code).medAbsMs,metrics.match.(code).p95AbsMs));
end
sgtitle(sprintf('Residuals vs Time | %s/%s',pairName,phase),'Interpreter','none');
saveas(f2, fullfile(outDir,sprintf('%s_02_residuals_vs_time.png', fileTag))); close(f2);

f3 = figure('Visible','off','Color','w','Position',[140 140 1200 450]);
for i = 1:3
    code = codes{i}; subplot(1,3,i);
    rs = metrics.match.(code).residualsMs;
    if isempty(rs)
        text(0.5,0.5,'No matches','HorizontalAlignment','center'); axis off;
    else
        histogram(rs,40,'FaceColor',clr(i,:),'EdgeColor','none'); hold on;
        xline(0,'k-'); grid on; xlabel('\Delta t (ms)'); ylabel('Count');
    end
    title(code);
end
sgtitle(sprintf('Residual Histograms | %s/%s',pairName,phase),'Interpreter','none');
saveas(f3, fullfile(outDir,sprintf('%s_03_residual_hist.png', fileTag))); close(f3);
end

function tag = makeFileTag(pairName, phase)
raw = sprintf('%s_%s', pairName, phase);
tag = regexprep(raw, '[^A-Za-z0-9_-]', '_');
end

function [tMatch, rs] = residualsVsTime(pred, ref, code, tolSec)
tPred = pred.times(strcmp(pred.codes, code));
tMff = ref.codeTimes.(code);
[rs,~,~,~,~,idxMff] = pairResidualsMs(tPred,tMff,tolSec);
tMatch = tMff(idxMff);
if numel(tMatch) ~= numel(rs)
    n = min(numel(tMatch),numel(rs)); tMatch = tMatch(1:n); rs = rs(1:n);
end
end

function [resMs, nMatch, recall, precision, idxPred, idxMff] = pairResidualsMs(tPred, tMff, tolSec)
tPred = tPred(:)'; tMff = tMff(:)';
i=1; j=1; resMs=[]; idxPred=[]; idxMff=[];
while i<=numel(tPred) && j<=numel(tMff)
    d = tPred(i)-tMff(j);
    if abs(d)<=tolSec
        resMs(end+1)=d*1000; idxPred(end+1)=i; idxMff(end+1)=j; %#ok<AGROW>
        i=i+1; j=j+1;
    elseif d < -tolSec
        i=i+1;
    else
        j=j+1;
    end
end
nMatch = numel(resMs);
if isempty(tMff), recall=NaN; else, recall=nMatch/numel(tMff); end
if isempty(tPred), precision=NaN; else, precision=nMatch/numel(tPred); end
end

function p = prct(v, q)
v = sort(v(:));
if isempty(v), p=NaN; return; end
if q<=0, p=v(1); return; end
if q>=100, p=v(end); return; end
pos = (q/100)*(numel(v)-1)+1; lo=floor(pos); hi=ceil(pos);
if lo==hi, p=v(lo); else, p=v(lo)+(pos-lo)*(v(hi)-v(lo)); end
end

function tfm = isMovementEntry(e)
txt = lower(strjoin({e.pairName, e.phase, e.mffPath, e.otbPath}, ' '));
hasMove = contains(txt,'reach') || contains(txt,'reaching') || ...
          contains(txt,'movement') || contains(txt,'movimento') || contains(txt,'motor');
hasNonMove = contains(txt,'rest') || contains(txt,'resting') || ...
             contains(txt,'mmn') || contains(txt,'ssr') || ...
             contains(txt,'cmv') || contains(txt,'mvc');
tfm = hasMove && ~hasNonMove;
end

function entries = collectInputEntries(inputRoot)
entries = struct('pairName',{},'phase',{},'mffPath',{},'otbPath',{},'state',{},'msg',{},'layout',{});
d = dir(inputRoot);
d = d([d.isdir]);
d = d(~ismember({d.name},{'.','..'}));
if isempty(d)
    return;
end

names = {d.name};
isSub = startsWith(names, 'sub-');
isSes = startsWith(names, 'ses-');
if any(isSub) || any(isSes)
    entries = collectBidsEntries(inputRoot, d(isSub), d(isSes));
else
    entries = collectLegacyEntries(inputRoot, d);
end
end

function entries = collectLegacyEntries(inputRoot, pairDirs)
entries = struct('pairName',{},'phase',{},'mffPath',{},'otbPath',{},'state',{},'msg',{},'layout',{});
phases = {'pre','post'};
for i = 1:numel(pairDirs)
    pairName = pairDirs(i).name;
    pairPath = fullfile(inputRoot, pairName);
    for ip = 1:numel(phases)
        phase = phases{ip};
        [mffPath, otbPath, state, msg] = findPhaseFiles(pairPath, phase);
        entries(end+1) = struct( ... %#ok<AGROW>
            'pairName', pairName, ...
            'phase', phase, ...
            'mffPath', mffPath, ...
            'otbPath', otbPath, ...
            'state', state, ...
            'msg', msg, ...
            'layout', 'legacy_pair');
    end
end
end

function entries = collectBidsEntries(inputRoot, subDirs, rootSesDirs)
entries = struct('pairName',{},'phase',{},'mffPath',{},'otbPath',{},'state',{},'msg',{},'layout',{});
for i = 1:numel(subDirs)
    subName = subDirs(i).name;
    subPath = fullfile(inputRoot, subName);
    sesDirs = dir(fullfile(subPath, 'ses-*'));
    sesDirs = sesDirs([sesDirs.isdir]);

    if isempty(sesDirs)
        phase = inferPhaseLabel('', subPath, '', '');
        e = buildEntriesForSession(subPath, subName, phase, 'bids');
        entries = [entries e]; %#ok<AGROW>
        continue;
    end

    for s = 1:numel(sesDirs)
        sesName = sesDirs(s).name;
        sesPath = fullfile(subPath, sesName);
        phase = inferPhaseLabel(sesName, sesPath, '', '');
        e = buildEntriesForSession(sesPath, sprintf('%s_%s', subName, sesName), phase, 'bids');
        entries = [entries e]; %#ok<AGROW>
    end
end

for s = 1:numel(rootSesDirs)
    sesName = rootSesDirs(s).name;
    sesPath = fullfile(inputRoot, sesName);
    phase = inferPhaseLabel(sesName, sesPath, '', '');
    e = buildEntriesForSession(sesPath, ['bids_' sesName], phase, 'bids_root_ses');
    entries = [entries e]; %#ok<AGROW>
end

if isempty(entries)
    phase = inferPhaseLabel('', inputRoot, '', '');
    e = buildEntriesForSession(inputRoot, 'bids_root', phase, 'bids_flat');
    entries = [entries e]; %#ok<AGROW>
end
end

function entries = buildEntriesForSession(sessionPath, pairPrefix, phase, layoutLabel)
entries = struct('pairName',{},'phase',{},'mffPath',{},'otbPath',{},'state',{},'msg',{},'layout',{});
[mffPaths, otbPaths] = listBidsFiles(sessionPath);

if isempty(mffPaths) && isempty(otbPaths)
    entries(end+1) = mkEntry(pairPrefix, phase, '', '', 'SKIP', ...
        'empty BIDS session (no EEG/EMG files)', layoutLabel); %#ok<AGROW>
    return;
end
if isempty(mffPaths)
    entries(end+1) = mkEntry(pairPrefix, phase, '', '', 'SKIP', ...
        'missing EEG .mff in BIDS session', layoutLabel); %#ok<AGROW>
    return;
end
if isempty(otbPaths)
    entries(end+1) = mkEntry(pairPrefix, phase, '', '', 'SKIP', ...
        'missing EMG .otb+ in BIDS session', layoutLabel); %#ok<AGROW>
    return;
end

[matchIdx, usedOtb] = matchFilesByName(mffPaths, otbPaths);
for i = 1:numel(mffPaths)
    runTag = sanitizeName(stripExt(mffPaths{i}, '.mff'));
    pairName = sprintf('%s_%s', pairPrefix, runTag);
    if matchIdx(i) > 0
        entries(end+1) = mkEntry(pairName, phase, mffPaths{i}, otbPaths{matchIdx(i)}, ...
            'OK', '', layoutLabel); %#ok<AGROW>
    else
        entries(end+1) = mkEntry(pairName, phase, mffPaths{i}, '', ...
            'SKIP', 'no matching EMG .otb+ for this EEG file in BIDS session', layoutLabel); %#ok<AGROW>
    end
end

orph = find(~usedOtb);
for k = 1:numel(orph)
    j = orph(k);
    runTag = sanitizeName(stripExt(otbPaths{j}, '.otb+'));
    pairName = sprintf('%s_emgonly_%s', pairPrefix, runTag);
    entries(end+1) = mkEntry(pairName, phase, '', otbPaths{j}, ...
        'SKIP', 'no matching EEG .mff for this EMG file in BIDS session', layoutLabel); %#ok<AGROW>
end
end

function e = mkEntry(pairName, phase, mffPath, otbPath, state, msg, layoutLabel)
e = struct('pairName',pairName,'phase',phase,'mffPath',mffPath,'otbPath',otbPath, ...
    'state',state,'msg',msg,'layout',layoutLabel);
end

function [mffPaths, otbPaths] = listBidsFiles(rootPath)
mffList = dir(fullfile(rootPath,'**','*.mff'));
mffList = mffList([mffList.isdir]);
otbList = dir(fullfile(rootPath,'**','*.otb+'));
otbList = otbList(~[otbList.isdir]);

mffPaths = cell(numel(mffList),1);
for i = 1:numel(mffList)
    mffPaths{i} = fullfile(mffList(i).folder, mffList(i).name);
end
otbPaths = cell(numel(otbList),1);
for i = 1:numel(otbList)
    otbPaths{i} = fullfile(otbList(i).folder, otbList(i).name);
end
mffPaths = sort(mffPaths);
otbPaths = sort(otbPaths);
end

function [matchIdx, usedOtb] = matchFilesByName(mffPaths, otbPaths)
nM = numel(mffPaths);
nO = numel(otbPaths);
matchIdx = zeros(1, nM);
usedOtb = false(1, nO);
for i = 1:nM
    bestScore = -inf;
    bestJ = 0;
    for j = 1:nO
        if usedOtb(j), continue; end
        sc = nameSimilarity(mffPaths{i}, otbPaths{j});
        if sc > bestScore
            bestScore = sc;
            bestJ = j;
        end
    end
    if bestScore >= 1
        matchIdx(i) = bestJ;
        usedOtb(bestJ) = true;
    end
end
end

function sc = nameSimilarity(pathA, pathB)
tA = normalizeTokens(pathA);
tB = normalizeTokens(pathB);
if isempty(tA) || isempty(tB)
    sc = 0;
    return;
end
sc = numel(intersect(tA, tB));
sideA = extractSideToken(tA);
sideB = extractSideToken(tB);
if ~isempty(sideA) && ~isempty(sideB)
    if strcmp(sideA, sideB)
        sc = sc + 2;
    else
        sc = sc - 2;
    end
end
end

function toks = normalizeTokens(pathIn)
base = lower(stripExt(pathIn, ''));
base = regexprep(base, '[^a-z0-9]+', ' ');
raw = regexp(strtrim(base), '\s+', 'split');
stop = {'sub','ses','eeg','emg','task','run','acq','rec','mff','otb','pre','post','session','bids'};
toks = {};
for i = 1:numel(raw)
    t = mapToken(raw{i});
    if isempty(t), continue; end
    if any(strcmp(t, stop)), continue; end
    if ~isempty(regexp(t, '^\d+$', 'once')), continue; end
    if ~isempty(regexp(t, '^[a-z]+\d+$', 'once')), continue; end
    toks{end+1} = t; %#ok<AGROW>
end
toks = unique(toks);
end

function t = mapToken(t)
if isempty(t), return; end
if startsWith(t, 'reach')
    t = 'reach'; return;
end
switch t
    case {'sin','sx','left','lt'}
        t = 'sx';
    case {'des','dx','right','rt'}
        t = 'dx';
    case {'deltroide','deltoide'}
        t = 'deltoide';
    case {'biciite','bicipite','bic'}
        t = 'bicipite';
end
end

function s = extractSideToken(toks)
if any(strcmp(toks,'dx'))
    s = 'dx';
elseif any(strcmp(toks,'sx'))
    s = 'sx';
else
    s = '';
end
end

function s = sanitizeName(s)
s = lower(char(s));
s = regexprep(s, '[^a-z0-9]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'run'; end
end

function phase = inferPhaseLabel(sesName, sesPath, mffPath, otbPath)
txt = lower(strjoin({sesName, sesPath, mffPath, otbPath}, ' '));
if contains(txt, 'pre')
    phase = 'pre';
elseif contains(txt, 'post')
    phase = 'post';
elseif startsWith(lower(sesName), 'ses-')
    phase = erase(lower(sesName), 'ses-');
else
    phase = 'session';
end
phase = regexprep(phase, '[^A-Za-z0-9_-]', '_');
if isempty(phase), phase = 'session'; end
end

function [mffPath, otbPath, state, msg] = findPhaseFiles(pairPath, phase)
phaseRoot = fullfile(pairPath, phase);
eegRoot = fullfile(phaseRoot, 'eeg');
emgRoot = fullfile(phaseRoot, 'emg');
if ~isfolder(phaseRoot)
    mffPath=''; otbPath=''; state='SKIP'; msg='phase folder missing'; return;
end

mffList = dir(fullfile(eegRoot,'**','*.mff')); mffList = mffList([mffList.isdir]);
otbList = dir(fullfile(emgRoot,'**','*.otb+')); otbList = otbList(~[otbList.isdir]);

if isempty(mffList) && isempty(otbList)
    mffPath=''; otbPath=''; state='SKIP'; msg='phase empty (no EEG/EMG files)'; return;
end
if isempty(mffList)
    mffPath=''; otbPath=''; state='SKIP'; msg='missing EEG .mff in phase'; return;
end
if isempty(otbList)
    mffPath=''; otbPath=''; state='SKIP'; msg='missing EMG .otb+ in phase'; return;
end
if numel(mffList) > 1
    mffPath=''; otbPath=''; state='ERR'; msg='multiple EEG .mff in phase'; return;
end
if numel(otbList) > 1
    mffPath=''; otbPath=''; state='ERR'; msg='multiple EMG .otb+ in phase'; return;
end

mffPath = fullfile(mffList(1).folder, mffList(1).name);
otbPath = fullfile(otbList(1).folder, otbList(1).name);
state='OK'; msg='';
end

function filePath = findOne(rootPath, patt, wantDir)
if ~isfolder(rootPath), error('Folder not found: %s', rootPath); end
lst = dir(fullfile(rootPath,'**',patt));
if wantDir, lst = lst([lst.isdir]); else, lst = lst(~[lst.isdir]); end
if numel(lst)~=1
    error('Expected exactly 1 %s in %s, found %d', patt, rootPath, numel(lst));
end
filePath = fullfile(lst(1).folder, lst(1).name);
end

function ref = readMffReference(mffPath, eventFile)
infoTxt = fileread(fullfile(mffPath,'info.xml'));
epochsTxt = fileread(fullfile(mffPath,'epochs.xml'));
evPath = fullfile(mffPath,eventFile);
evTxt = fileread(evPath);
recordIso = oneToken(infoTxt,'<recordTime>([^<]+)</recordTime>');
recordEpoch = isoToEpoch(recordIso);
name = oneToken(evTxt,'<name>([^<]+)</name>');
trackType = oneToken(evTxt,'<trackType>([^<]+)</trackType>');
endUs = str2double(oneToken(epochsTxt,'<endTime>(\d+)</endTime>'));
durationS = endUs / 1e6;
allEv = regexp(evTxt,'<event>\s*<beginTime>([^<]+)</beginTime>\s*<duration>([^<]+)</duration>\s*<code>([^<]+)</code>[\s\S]*?</event>','tokens');
if isempty(allEv), error('No events in %s', evPath); end
codeTimes = struct('DIN1',[],'DIN4',[],'DIN5',[]);
for i=1:numel(allEv)
    tIso = allEv{i}{1}; code = allEv{i}{3};
    t = isoToEpoch(tIso) - recordEpoch;
    if isfield(codeTimes,code), codeTimes.(code)(end+1)=t; end %#ok<AGROW>
end
ref = struct('mffPath',mffPath,'eventPath',evPath,'recordEpoch',recordEpoch, ...
    'name',name,'trackType',trackType,'durationS',durationS,'codeTimes',codeTimes);
end

function otb = readOtb(otbPath)
tmpDir = tempname; mkdir(tmpDir); c = onCleanup(@() cleanupTemp(tmpDir)); %#ok<NASGU>
untar(otbPath,tmpDir);
xmlList = dir(fullfile(tmpDir,'*_08.xml')); if isempty(xmlList), xmlList=dir(fullfile(tmpDir,'*.xml')); end
sigList = dir(fullfile(tmpDir,'*.sig'));
if isempty(xmlList) || isempty(sigList), error('Cannot parse OTB: %s', otbPath); end
xmlPath = fullfile(xmlList(1).folder, xmlList(1).name);
sigPath = fullfile(sigList(1).folder, sigList(1).name);

doc = xmlread(xmlPath); root = doc.getDocumentElement();
fs = str2double(char(root.getAttribute('SampleFrequency')));
nChannels = str2double(char(root.getAttribute('DeviceTotalChannels')));
auxMap = containers.Map('KeyType','double','ValueType','double');
ad = root.getElementsByTagName('Adapter');
for i=0:(ad.getLength()-1)
    a = ad.item(i); st = str2double(char(a.getAttribute('ChannelStartIndex')));
    chNodes = a.getElementsByTagName('Channel');
    for j=0:(chNodes.getLength()-1)
        ch = chNodes.item(j);
        prefix = char(ch.getAttribute('Prefix'));
        idStr = char(ch.getAttribute('ID'));
        idx = str2double(char(ch.getAttribute('Index')));
        if contains(prefix,'AUX') && contains(idStr,'Trigger')
            tk = regexp(prefix,'AUX\s+(\d+)','tokens','once');
            if ~isempty(tk), auxMap(str2double(tk{1})) = st + idx; end
        end
    end
end
fid = fopen(sigPath,'r'); raw = fread(fid,inf,'int16=>double',0,'l'); fclose(fid);
N = floor(numel(raw)/nChannels); raw = raw(1:N*nChannels); data = reshape(raw,[nChannels,N]);
otb = struct('fs',fs,'nChannels',nChannels,'data',data,'auxMap',auxMap);
end

function [events, meta] = buildAlignedEvents(ref, otb)
needed = [2 3 4];
for k = needed
    if ~isKey(otb.auxMap,k), error('Missing AUX%d',k); end
end
edges = struct();
raw = struct();
for k = needed
    x = otb.data(otb.auxMap(k)+1,:);
    t = detectFallingEdges(x,otb.fs);
    edges.(sprintf('AUX%d',k)) = t;
    raw.(sprintf('AUX%d',k)) = numel(t);
end

aux4 = edges.AUX4;
aux3 = edges.AUX3;
aux2 = edges.AUX2;
din4 = ref.codeTimes.DIN4;
din5 = ref.codeTimes.DIN5;
if isempty(din4) || isempty(din5), error('Reference MFF must contain DIN4 and DIN5'); end

fit1 = bestLagFit(aux4, aux3, din4, din5, 8);
fit2 = bestLagFit(aux4, aux3, din5, din4, 8);
if ~fit1.valid && ~fit2.valid
    error('Unable to fit AUX to DIN (insufficient aligned trigger pairs)');
end
if fit1.valid && (~fit2.valid || fit1.mad <= fit2.mad)
    fitSel = fit1;
    map4 = 'DIN4'; map3 = 'DIN5';
else
    fitSel = fit2;
    map4 = 'DIN5'; map3 = 'DIN4';
end

a = fitSel.a;
b = fitSel.b;
mad = fitSel.mad;

keyO = fitSel.x4;
keyM = fitSel.y4;
beepO = fitSel.x3;
beepM = fitSel.y3;
nCommon = min(numel(keyO), numel(beepO));

fitInfo = struct('stage','lag_search_initial', ...
    'lag4',fitSel.lag4,'lag3',fitSel.lag3, ...
    'n_fit_trials_candidate',nCommon,'n_fit_trials_kept',0, ...
    'refit_on_valid_trials',false, ...
    'refit_with_din1',false, ...
    'n_din1_match',0, ...
    'n_din1_used',0, ...
    'din1_tol_ms',NaN, ...
    'din1_match_ratio',NaN);

idxKeep = [];
if nCommon >= 2
    tKeyPred = a*keyO + b;
    tBeepPred = a*beepO + b;
    % Trial-gating is based on keyboard/beep order only (DIN4/DIN5).
    idxKeep = findValidTrialIndicesKbBeep(tKeyPred, tBeepPred);
    fitInfo.n_fit_trials_kept = numel(idxKeep);
    if numel(idxKeep) >= 2
        xFit = [keyO(idxKeep) beepO(idxKeep)];
        yFit = [keyM(idxKeep) beepM(idxKeep)];
        [aV,bV,mV] = fitLinear(xFit, yFit);

        [xFitS,ordV] = sort(xFit);
        yFitS = yFit(ordV);
        if numel(xFitS) >= 2 && xFitS(end) > xFitS(1)
            aEpV = (yFitS(end)-yFitS(1))/(xFitS(end)-xFitS(1));
            bEpV = yFitS(1) - aEpV*xFitS(1);
            mEpV = median(abs(yFitS - (aEpV*xFitS+bEpV)));
            if mEpV <= mV
                aV = aEpV; bV = bEpV; mV = mEpV;
            end
        end

        a = aV; b = bV; mad = mV;
        fitInfo.stage = 'valid_trials_refit';
        fitInfo.refit_on_valid_trials = true;
    end
end

if ~isempty(idxKeep) && numel(idxKeep) >= 2
    xBase = [keyO(idxKeep) beepO(idxKeep)];
    yBase = [keyM(idxKeep) beepM(idxKeep)];
else
    xBase = [keyO beepO];
    yBase = [keyM beepM];
end

% Optional Sync refinement using dense DIN1 train. This improves clock fit
% but does not change trial validity logic, which remains DIN4/DIN5-only.
[xDin1, yDin1, din1Fit] = pairDin1ForRefit(aux2, ref.codeTimes.DIN1, a, b);
fitInfo.n_din1_match = din1Fit.nMatched;
fitInfo.n_din1_used = din1Fit.nUsed;
fitInfo.din1_tol_ms = din1Fit.tolS * 1000;
fitInfo.din1_match_ratio = din1Fit.matchRatio;
if din1Fit.used
    [aD, bD, mD] = fitLinear([xBase xDin1], [yBase yDin1]);
    if isfinite(mD)
        a = aD;
        b = bD;
        mad = mD;
        fitInfo.refit_with_din1 = true;
        if fitInfo.refit_on_valid_trials
            fitInfo.stage = 'valid_trials_plus_din1_refit';
        else
            fitInfo.stage = 'lag_search_plus_din1_refit';
        end
    end
end

t = [];
c = {};
for v = aux4
    tt = a*v+b;
    if tt>=0 && tt<=ref.durationS+0.1
        t(end+1) = tt; c{end+1} = map4; %#ok<AGROW>
    end
end
for v = aux3
    tt = a*v+b;
    if tt>=0 && tt<=ref.durationS+0.1
        t(end+1) = tt; c{end+1} = map3; %#ok<AGROW>
    end
end
for v = aux2
    tt = a*v+b;
    if tt>=0 && tt<=ref.durationS+0.1
        t(end+1) = tt; c{end+1} = 'DIN1'; %#ok<AGROW>
    end
end
[t,o] = sort(t); c = c(o);

events = struct('times',t,'codes',{c});
meta = struct('a',a,'b',b,'mad',mad,'driftPpm',(a-1)*1e6, ...
    'mapping',struct('AUX4',map4,'AUX3',map3,'AUX2','DIN1'), ...
    'rawCounts',raw, ...
    'fit',fitInfo);
end

function t = detectFallingEdges(x,fs)
head=sort(x(1:min(20000,numel(x)))); med=head(floor(numel(head)/2)+1); thr=(med+min(x))/2;
mask=x<thr; on=find(diff([false mask])==1); t=(on-1)/fs;
end

function [a,b,mad] = fitLinear(x,y)
x=double(x(:)); y=double(y(:)); p=polyfit(x,y,1); a=p(1); b=p(2); mad=median(abs(y-(a*x+b)));
end

function idxKeep = findValidTrialIndicesKbBeep(din4, din5)
nTrials = min(numel(din4), numel(din5));
idxKeep = [];
if nTrials == 0
    return;
end
for k = 1:nTrials
    tKey = din4(k);
    tBeep = din5(k);
    if ~(isfinite(tKey) && isfinite(tBeep) && (tKey < tBeep))
        continue;
    end
    if k < nTrials
        nextKey = din4(k+1);
    else
        nextKey = inf;
    end
    if ~(isfinite(nextKey) && (nextKey > tBeep)) && ~isinf(nextKey)
        continue;
    end
    idxKeep(end+1) = k; %#ok<AGROW>
end
end

function fit = bestLagFit(aux4, aux3, din4Target, din5Target, maxLag)
if nargin < 5 || isempty(maxLag)
    maxLag = 8;
end
fit = struct('valid',false,'a',NaN,'b',NaN,'mad',Inf, ...
    'lag4',0,'lag3',0,'x4',[],'y4',[],'x3',[],'y3',[]);
bestScore = Inf;
for lag4 = -maxLag:maxLag
    [x4,y4] = pairWithLag(aux4,din4Target,lag4);
    if numel(x4) < 10, continue; end
    for lag3 = -maxLag:maxLag
        [x3,y3] = pairWithLag(aux3,din5Target,lag3);
        if numel(x3) < 10, continue; end
        x = [x4 x3];
        y = [y4 y3];
        [aT,bT,mT] = fitLinear(x,y);
        sc = mT + 1e-9*abs(aT-1);
        if sc < bestScore
            bestScore = sc;
            fit.valid = true;
            fit.a = aT;
            fit.b = bT;
            fit.mad = mT;
            fit.lag4 = lag4;
            fit.lag3 = lag3;
            fit.x4 = x4;
            fit.y4 = y4;
            fit.x3 = x3;
            fit.y3 = y3;
        end
    end
end
end

function [xUse, yUse, info] = pairDin1ForRefit(aux2, din1, a, b)
xUse = [];
yUse = [];
info = struct('used',false,'nMatched',0,'nUsed',0,'tolS',NaN,'matchRatio',NaN);

if isempty(aux2) || isempty(din1)
    return;
end

aux2 = sort(double(aux2(:)'));
din1 = sort(double(din1(:)'));
if numel(aux2) < 10 || numel(din1) < 10
    return;
end

predDin1 = a * aux2 + b;
din1Step = median(diff(din1));
if ~isfinite(din1Step) || din1Step <= 0
    din1Step = 0.1;
end
tolS = min(0.05, max(0.008, 0.35 * din1Step));

[idxPred, idxRef] = pairSortedByTolerance(predDin1, din1, tolS);
nMatched = numel(idxPred);
minCount = min(numel(predDin1), numel(din1));
matchRatio = nMatched / max(minCount, 1);

info.nMatched = nMatched;
info.tolS = tolS;
info.matchRatio = matchRatio;
if nMatched < 20 || matchRatio < 0.30
    return;
end

x = aux2(idxPred);
y = din1(idxRef);
r = y - (a*x + b);
if ~isempty(r)
    rMed = median(r);
    rMad = median(abs(r - rMed));
    thr = max(3*rMad, 0.003);
    keep = abs(r - rMed) <= thr;
    if nnz(keep) >= 20
        x = x(keep);
        y = y(keep);
    end
end

maxUsed = 400;
if numel(x) > maxUsed
    sel = round(linspace(1, numel(x), maxUsed));
    x = x(sel);
    y = y(sel);
end

xUse = x;
yUse = y;
info.nUsed = numel(xUse);
info.used = info.nUsed >= 20;
end

function [idxA, idxB] = pairSortedByTolerance(aVals, bVals, tolS)
idxA = [];
idxB = [];
i = 1;
j = 1;
while i <= numel(aVals) && j <= numel(bVals)
    d = aVals(i) - bVals(j);
    if abs(d) <= tolS
        idxA(end+1) = i; %#ok<AGROW>
        idxB(end+1) = j; %#ok<AGROW>
        i = i + 1;
        j = j + 1;
    elseif d < -tolS
        i = i + 1;
    else
        j = j + 1;
    end
end
end

function [x,y] = pairWithLag(aux,din,lag)
x = [];
y = [];
if isempty(aux) || isempty(din), return; end
if lag >= 0
    n = min(numel(aux)-lag, numel(din));
    if n <= 0, return; end
    x = aux((1+lag):(lag+n));
    y = din(1:n);
else
    n = min(numel(aux), numel(din)+lag);
    if n <= 0, return; end
    k = -lag;
    x = aux(1:n);
    y = din((1+k):(k+n));
end
end

function s = oneToken(txt,pat)
t=regexp(txt,pat,'tokens','once'); if isempty(t), error('Pattern not found: %s', pat); end; s=t{1};
end

function epochSec = isoToEpoch(isoStr)
d=datetime(isoStr,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','UTC'); epochSec=posixtime(d);
end

function writeCsv(pathCsv, header, rows)
fid=fopen(pathCsv,'w'); fprintf(fid,'%s\n',strjoin(header,','));
for i=1:size(rows,1)
    r=rows(i,:);
    for j=1:numel(r), r{j}=esc(r{j}); end
    fprintf(fid,'%s\n',strjoin(r,','));
end
fclose(fid);
end

function s = esc(v)
if ~ischar(v), v=char(string(v)); end
if contains(v,'"'), v=strrep(v,'"','""'); end
if contains(v,',') || contains(v,'"') || contains(v,newline), s=['"' v '"']; else, s=v; end
end

function s = tf(x)
if x, s='true'; else, s='false'; end
end

function n = stripExt(pathIn, ext)
[~,n,e] = fileparts(pathIn);
if ~strcmpi(e,ext)
    % keep base name only
end
end

function cleanupTemp(d)
if isfolder(d), try, rmdir(d,'s'); catch, end, end
end

function Tdiff = buildPrePostDiffTable(Tsum)
varNames = {'key','pair_pre','pair_post','csc_ready_pre','csc_ready_post', ...
    'csc_scarto_ms_pre','csc_scarto_ms_post','csc_scarto_ms_delta_post_minus_pre', ...
    'p95_mov_max_ms_pre','p95_mov_max_ms_post','p95_mov_max_ms_delta_post_minus_pre', ...
    'margine_csc_ms_pre','margine_csc_ms_post','margine_csc_ms_delta_post_minus_pre', ...
    'din4_p95_abs_ms_pre','din4_p95_abs_ms_post','din4_p95_abs_ms_delta_post_minus_pre', ...
    'din5_p95_abs_ms_pre','din5_p95_abs_ms_post','din5_p95_abs_ms_delta_post_minus_pre', ...
    'mad_fit_ms_pre','mad_fit_ms_post','mad_fit_ms_delta_post_minus_pre', ...
    'drift_ppm_pre','drift_ppm_post','drift_ppm_delta_post_minus_pre'};

if isempty(Tsum) || height(Tsum) == 0
    Tdiff = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
    return;
end

pair = asCellStr(Tsum.pair);
phase = lower(asCellStr(Tsum.phase));
status = upper(asCellStr(Tsum.status));

ok = strcmp(status,'OK') & (strcmp(phase,'pre') | strcmp(phase,'post'));
if ~any(ok)
    Tdiff = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
    return;
end

pair = pair(ok);
phase = phase(ok);
S = Tsum(ok,:);

keys = cellfun(@prePostKeyFromPair, pair, 'UniformOutput', false);
u = unique(keys);
rows = {};
for i = 1:numel(u)
    k = u{i};
    iPre = find(strcmp(keys,k) & strcmp(phase,'pre'), 1, 'first');
    iPost = find(strcmp(keys,k) & strcmp(phase,'post'), 1, 'first');
    if isempty(iPre) || isempty(iPost)
        continue;
    end

    prePair = pair{iPre};
    postPair = pair{iPost};
    cscPre = asCellStr(S.csc_ready); cscPre = cscPre{iPre};
    cscPost = asCellStr(S.csc_ready); cscPost = cscPost{iPost};

    scPre = numCell(S.csc_scarto_ms, iPre); scPost = numCell(S.csc_scarto_ms, iPost);
    p95mPre = numCell(S.p95_mov_max_ms, iPre); p95mPost = numCell(S.p95_mov_max_ms, iPost);
    marPre = numCell(S.margine_csc_ms, iPre); marPost = numCell(S.margine_csc_ms, iPost);
    d4Pre = numCell(S.din4_p95_abs_ms, iPre); d4Post = numCell(S.din4_p95_abs_ms, iPost);
    d5Pre = numCell(S.din5_p95_abs_ms, iPre); d5Post = numCell(S.din5_p95_abs_ms, iPost);
    madPre = numCell(S.mad_fit_ms, iPre); madPost = numCell(S.mad_fit_ms, iPost);
    drPre = numCell(S.drift_ppm, iPre); drPost = numCell(S.drift_ppm, iPost);

    rows(end+1,:) = { ... %#ok<AGROW>
        k, prePair, postPair, cscPre, cscPost, ...
        fmt3(scPre), fmt3(scPost), fmt3(scPost-scPre), ...
        fmt3(p95mPre), fmt3(p95mPost), fmt3(p95mPost-p95mPre), ...
        fmt3(marPre), fmt3(marPost), fmt3(marPost-marPre), ...
        fmt3(d4Pre), fmt3(d4Post), fmt3(d4Post-d4Pre), ...
        fmt3(d5Pre), fmt3(d5Post), fmt3(d5Post-d5Pre), ...
        fmt3(madPre), fmt3(madPost), fmt3(madPost-madPre), ...
        fmt3(drPre), fmt3(drPost), fmt3(drPost-drPre)};
end

if isempty(rows)
    Tdiff = cell2table(cell(0, numel(varNames)), 'VariableNames', varNames);
else
    Tdiff = cell2table(rows, 'VariableNames', varNames);
end
end

function k = prePostKeyFromPair(pairName)
s = lower(char(pairName));
s = regexprep(s, 'ses[-_]?pre', 'ses');
s = regexprep(s, 'ses[-_]?post', 'ses');
s = regexprep(s, '(^|[_-])pre([_-]|$)', '$1$2');
s = regexprep(s, '(^|[_-])post([_-]|$)', '$1$2');
s = regexprep(s, '[_-]+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'pair'; end
k = s;
end

function c = asCellStr(x)
if iscell(x)
    c = cell(size(x));
    for i = 1:numel(x)
        if ismissingVal(x{i})
            c{i} = '';
        else
            c{i} = char(string(x{i}));
        end
    end
elseif isstring(x)
    c = cellstr(x);
else
    c = cell(size(x));
    for i = 1:numel(x)
        c{i} = char(string(x(i)));
    end
end
end

function tfm = ismissingVal(v)
tfm = isempty(v) || (isstring(v) && all(ismissing(v)));
end

function v = numCell(col, idx)
if iscell(col)
    v = str2double(string(col{idx}));
else
    v = str2double(string(col(idx)));
end
if isnan(v), v = NaN; end
end

function s = fmt3(v)
if isnan(v)
    s = 'NaN';
else
    s = sprintf('%.3f', v);
end
end
