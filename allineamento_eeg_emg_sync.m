function allineamento_eeg_emg_sync(inputRoot, outputRoot, overwriteOutput, templateMffPath)
%ALLINEAMENTO_EEG_EMG_SYNC
% Batch alignment of EMG (OTB+) to EEG (MFF) using shared trigger signals.
%
% The alignment relies on trigger signals that are simultaneously recorded
% by both instruments:
%   EEG (MFF)   |   EMG (OTB+)
%   ---------------------------
%   DIN4        |   AUX4 (falling edge) -> keyboard event
%   DIN5        |   AUX3 (falling edge) -> beep event
%   DIN1        |   AUX2 (falling edge) -> sync wave (dense train)
%
% The expected task order is: start < keyboard (DIN4) < beep (DIN5) < end,
% with DIN4 always preceding DIN5.
%
% Supported input layouts:
% 1) Legacy pair structure
%   pair_xxx/
%     pre/
%       eeg/   (contains one .mff)
%       emg/   (contains one .otb+)
%     post/
%       eeg/   (contains one .mff)
%       emg/   (contains one .otb+)
% 2) BIDS-like structure
%   sub-XX/ses-YY/...  (one .mff + one .otb+ per session)
%   Supports multi-subject roots and, if needed, one-level parent folders
%   containing multiple BIDS datasets.
%
% If PRE or POST is empty/incomplete, phase is skipped (no fatal error).
%
% Trial gating note (updated logic):
% - trial candidates are selected from DIN4/DIN5 order only;
% - DIN1 is preserved as trigger channel and is never used as a trial-validity gate;
% - if a dense DIN1 train is present on both EEG and EMG, it can refine the
%   affine drift fit after the DIN4/DIN5-based trial fit.

% --- Default parameter values ---
% inputRoot:    parent folder containing pair_xxx/ subfolders or BIDS structure
% outputRoot:   where aligned .mff, .csv, report.json and movement .otb+ are saved
% overwriteOutput: if true, delete existing output before writing
% templateMffPath: optional reference MFF with pre-labeled events (start, keyboard,
%                  beep, end) used to derive temporal offsets for event reconstruction
if nargin < 1 || isempty(inputRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    inputRoot = fullfile(scriptDir, 'batch_otb_mff', 'input');
end
if nargin < 2 || isempty(outputRoot)
    scriptDir = fileparts(mfilename('fullpath'));
    outputRoot = fullfile(scriptDir, 'batch_otb_mff', 'output_matlab_sync');
end
if nargin < 3 || isempty(overwriteOutput)
    overwriteOutput = true;
end
if nargin < 4 || isempty(templateMffPath)
    candidate = '/Users/martinaregazzetti/Desktop/dati per giocare/test.mff';
    if isfolder(candidate)
        templateMffPath = candidate;
    else
        templateMffPath = '';
    end
end

% XML event file name inside each .mff folder
eventFile = 'Events_8 DINs.xml';
if ~isfolder(inputRoot), error('Input root not found: %s', inputRoot); end
if ~isfolder(outputRoot), mkdir(outputRoot); end

% Discover all valid (MFF, OTB) pairs to process based on the input layout
entries = collectInputEntries(inputRoot);
if isempty(entries), error('No valid input entries in %s', inputRoot); end

% summary.csv header: one row per processed (pair, phase)
header = {'pair','phase','status','message','source_mff','source_otb','output_mff','output_csv', ...
    'a','b','drift_ppm','mad_ms','din1','din4','din5','start','end','keyboard','beep','total', ...
    'trial_candidates','trial_kept','trial_removed','trial_removed_pct','trial_filter_reason'};
rows = {};
% trial_sequence_report.csv header: concise trial filtering summary per run
trialHeader = {'pair','phase','status','trial_candidates','trial_kept','trial_removed','trial_removed_pct','trial_filter_reason'};
trialRows = {};

% --- Main processing loop: one iteration per (pair, phase) ---
for i = 1:numel(entries)
    pairName = entries(i).pairName;
    phase = entries(i).phase;
    mffSrc = entries(i).mffPath;
    otbSrc = entries(i).otbPath;
    state = entries(i).state;
    msg = entries(i).msg;

    try
        % Skip entries with missing or incomplete files
        if ~strcmp(state,'OK')
            fprintf('[SKIP] %s/%s | %s\n', pairName, phase, msg);
            rows(end+1,:) = {pairName,phase,'SKIP',msg,'','','','','','','','','','','','','','','','','','','','',''}; %#ok<AGROW>
            trialRows(end+1,:) = {pairName,phase,'SKIP','','','','',msg}; %#ok<AGROW>
            continue;
        end

        % Create output directory for this (pair, phase)
        pairOut = fullfile(outputRoot, pairName, phase);
        if ~isfolder(pairOut), mkdir(pairOut); end

        % Copy source MFF to output before modifying events
        outMffName = [stripExt(mffSrc,'.mff') '_aligned_otb.mff'];
        outMff = fullfile(pairOut, outMffName);
        if isfolder(outMff)
            if overwriteOutput
                rmdir(outMff,'s');
            else
                error('Output exists: %s', outMff);
            end
        end
        copyfile(mffSrc, outMff);

        % Step 1: read EEG reference (DIN events from MFF) and EMG data (AUX from OTB)
        ref = readMffReference(outMff, eventFile);
        otb = readOtb(otbSrc);
        % Step 2: align AUX falling edges to DIN events via affine fit,
        %         reconstruct task events (start, keyboard, beep, end)
        [events, meta] = buildAlignedEvents(ref, otb, templateMffPath);
        % Step 3: write aligned events back into the MFF XML
        writeEventXml(ref, events);
        % Step 4: export event list as CSV
        outCsv = fullfile(pairOut, [stripExt(outMff,'.mff') '_events.csv']);
        writeEventsCsv(outCsv, events);
        % Step 5: extract movement-only OTB segments from valid trial windows
        [outMovementOtb, movementSegCsv, movementInfo] = exportMovementOnlyOtb(otbSrc, pairOut, meta, otb);

        % Build and write JSON report with full alignment metadata
        report = struct();
        report.pair = pairName;
        report.phase = phase;
        report.source_mff = mffSrc;
        report.source_otb = otbSrc;
        report.output_mff = outMff;
        report.output_csv = outCsv;
        report.input_layout = entries(i).layout;
        report.transform = struct('a',meta.a,'b',meta.b,'drift_ppm',meta.driftPpm,'mad_ms',meta.mad*1000);
        report.fit = meta.fit;
        report.mapping_aux_to_din = meta.mapping;
        report.raw_otb_counts = meta.rawCounts;
        if isfield(meta,'cleanCounts'), report.clean_otb_counts = meta.cleanCounts; end
        if isfield(meta,'cleaning'), report.trigger_cleaning = meta.cleaning; end
        report.logical = meta.logical;
        report.logic_description_en = ['Trial gating uses DIN4/DIN5 only; DIN1 is never used as trial validity gate. ' ...
            'When a dense DIN1 train is available on both EEG and EMG, it is used only to refine clock-drift fit (a,b).'];
        report.output_movement_otb = outMovementOtb;
        report.output_movement_segments_csv = movementSegCsv;
        report.movement_only = movementInfo;
        report.movement_logic_description_en = ['Movement-only OTB is built by concatenating valid trial windows from reconstructed triggers. ' ...
            'Primary rule: beep -> start; automatic fallback: beep -> end when start does not occur after beep.'];
        report.final_counts = meta.finalCounts;
        report.total_events = numel(events.times);
        report.generated_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        reportPath = fullfile(pairOut,'report.json');
        fid = fopen(reportPath,'w'); fwrite(fid,jsonencode(report),'char'); fclose(fid);

        % Extract trial filtering statistics for summary output
        nCand = NaN; nKeep = NaN; nRem = NaN; remPct = NaN; remReason = '';
        if isfield(meta,'logical') && isstruct(meta.logical)
            lm = meta.logical;
            if isfield(lm,'n_trials_candidate'), nCand = double(lm.n_trials_candidate); end
            if isfield(lm,'n_trials'), nKeep = double(lm.n_trials); end
            if isfield(lm,'n_trials_removed'), nRem = double(lm.n_trials_removed); end
            if isfield(lm,'reason'), remReason = char(string(lm.reason)); end
        end
        if ~isnan(nCand) && nCand > 0 && ~isnan(nRem)
            remPct = 100 * nRem / nCand;
        end

        if isnan(nCand), nCandS = ''; else, nCandS = sprintf('%d', round(nCand)); end
        if isnan(nKeep), nKeepS = ''; else, nKeepS = sprintf('%d', round(nKeep)); end
        if isnan(nRem), nRemS = ''; else, nRemS = sprintf('%d', round(nRem)); end
        if isnan(remPct), remPctS = ''; else, remPctS = sprintf('%.2f', remPct); end

        fprintf('[OK] %s/%s | events=%d | drift=%.2f ppm | MAD=%.3f ms | trials kept=%s/%s removed=%s | move_seg=%d\n', ...
            pairName, phase, numel(events.times), meta.driftPpm, meta.mad*1000, nKeepS, nCandS, nRemS, movementInfo.n_segments);

        % Append row to summary.csv
        rows(end+1,:) = {pairName,phase,'OK','',mffSrc,otbSrc,outMff,outCsv, ...
            sprintf('%.12f',meta.a),sprintf('%.12f',meta.b),sprintf('%.2f',meta.driftPpm),sprintf('%.6f',meta.mad*1000), ...
            num2str(getCount(meta.finalCounts,'DIN1')),num2str(getCount(meta.finalCounts,'DIN4')),num2str(getCount(meta.finalCounts,'DIN5')), ...
            num2str(getCount(meta.finalCounts,'start')),num2str(getCount(meta.finalCounts,'end')), ...
            num2str(getCount(meta.finalCounts,'keyboard')),num2str(getCount(meta.finalCounts,'beep')), ...
            num2str(numel(events.times)),nCandS,nKeepS,nRemS,remPctS,remReason}; %#ok<AGROW>
        % Append row to trial_sequence_report.csv
        trialRows(end+1,:) = {pairName,phase,'OK',nCandS,nKeepS,nRemS,remPctS,remReason}; %#ok<AGROW>

    catch ME
        % Log error and continue with next entry
        fprintf('[ERR] %s/%s | %s\n', pairName, phase, ME.message);
        rows(end+1,:) = {pairName,phase,'ERR',ME.message,'','','','','','','','','','','','','','','','','','','','',''}; %#ok<AGROW>
        trialRows(end+1,:) = {pairName,phase,'ERR','','','','',ME.message}; %#ok<AGROW>
    end
end

% Write summary CSV with one row per processed (pair, phase)
summaryPath = fullfile(outputRoot,'summary.csv');
writeCsv(summaryPath, header, rows);
% Write trial sequence report with trial filtering counts
trialReportPath = fullfile(outputRoot,'trial_sequence_report.csv');
writeCsv(trialReportPath, trialHeader, trialRows);

% Build unified drift report with three row types:
% - 'session':       single drift measurement per (pair, phase)
% - 'pre_post':      delta drift (post - pre) for same-subject pairs
% - 'right_left':    delta drift (sx - dx) within the same phase
driftHeader = {'report_type','key','phase','pair_a','pair_b','side_a','side_b', ...
    'drift_ppm_a','drift_ppm_b','delta_ppm_b_minus_a','abs_delta_ppm','trend', ...
    'n_side_a_candidates','n_side_b_candidates','notes'};
driftRows = buildUnifiedDriftRows(header, rows);
driftPath = fullfile(outputRoot,'drift_report.csv');
writeCsv(driftPath, driftHeader, driftRows);

fprintf('\nSummary: %s\n', summaryPath);
fprintf('Trial report: %s\n', trialReportPath);
fprintf('Drift report: %s\n', driftPath);
end

function [mffPath, otbPath, state, msg] = findPhaseFiles(pairPath, phase)
% Find exactly one .mff and one .otb+ inside a legacy pair_xxx/(pre|post) folder.
% Returns:
%   mffPath, otbPath -- full file paths (empty if not found)
%   state            -- 'OK', 'SKIP', or 'ERR'
%   msg              -- description of the issue, if any
phaseRoot = fullfile(pairPath, phase);
eegRoot = fullfile(phaseRoot, 'eeg');
emgRoot = fullfile(phaseRoot, 'emg');

if ~isfolder(phaseRoot)
    mffPath=''; otbPath=''; state='SKIP'; msg='phase folder missing'; return;
end

mffList = dir(fullfile(eegRoot,'**','*.mff'));
mffList = mffList([mffList.isdir]);  % .mff is a directory
otbList = dir(fullfile(emgRoot,'**','*.otb+'));
otbList = otbList(~[otbList.isdir]);  % .otb+ is a file

% Check for missing or multiple files
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
state = 'OK';
msg = '';
end

function entries = collectInputEntries(inputRoot)
% Discover all (MFF, OTB) pairs to process.
% Supports two layouts:
%   - BIDS: sub-XX/ses-YY/... (auto-detected if 'sub-' or 'ses-' folders exist)
%   - Legacy: pair_xxx/(pre|post)/{eeg,emg}/
% Returns an array of entry structs with fields: pairName, phase, mffPath,
% otbPath, state, msg, layout.
entries = struct('pairName',{},'phase',{},'mffPath',{},'otbPath',{},'state',{},'msg',{},'layout',{});
bidsRoots = findBidsRoots(inputRoot);
if ~isempty(bidsRoots)
    % BIDS layout found: collect entries for each BIDS root
    for r = 1:numel(bidsRoots)
        rootPath = bidsRoots{r};
        d = dir(rootPath);
        d = d([d.isdir]);
        d = d(~ismember({d.name},{'.','..'}));
        names = {d.name};
        isSub = startsWith(names, 'sub-');
        isSes = startsWith(names, 'ses-');
        e = collectBidsEntries(rootPath, d(isSub), d(isSes));
        if numel(bidsRoots) > 1
            dsTag = sanitizeName(stripExt(rootPath, ''));
            for k = 1:numel(e)
                e(k).pairName = sprintf('%s_%s', dsTag, e(k).pairName);
                e(k).layout = [e(k).layout '_multiroot'];
            end
        end
        entries = [entries e]; %#ok<AGROW>
    end
    return;
end

% Fallback to legacy pair structure
d = dir(inputRoot);
d = d([d.isdir]);
d = d(~ismember({d.name},{'.','..'}));
if isempty(d)
    return;
end
entries = collectLegacyEntries(inputRoot, d);
end

function roots = findBidsRoots(inputRoot)
% Detect BIDS-like directory structure.
% Returns a cell array of paths containing 'sub-*' or 'ses-*' folders.
% Supports both direct BIDS roots and parent folders with multiple BIDS datasets.
roots = {};

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
    roots = {inputRoot};
    return;
end

% Also support passing a parent folder that contains multiple BIDS roots.
for i = 1:numel(d)
    p = fullfile(inputRoot, d(i).name);
    dd = dir(p);
    dd = dd([dd.isdir]);
    dd = dd(~ismember({dd.name},{'.','..'}));
    if isempty(dd)
        continue;
    end
    nn = {dd.name};
    if any(startsWith(nn,'sub-')) || any(startsWith(nn,'ses-'))
        roots{end+1} = p; %#ok<AGROW>
    end
end

if ~isempty(roots)
    roots = sort(unique(roots));
end
end

function entries = collectLegacyEntries(inputRoot, pairDirs)
% Collect entries for legacy pair_xxx/(pre|post)/{eeg,emg}/ structure.
% Each pair is processed for both 'pre' and 'post' phases.
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
% Collect entries from BIDS structure: sub-XX/ses-YY/...
% Handles subjects with/without sessions, root-level sessions, and flat BIDS.
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
% Scan a BIDS session directory for .mff and .otb+ files and create
% (pairName, phase) entries, matching EEG to EMG by filename similarity.
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

% Match EEG (.mff) to EMG (.otb+) by name similarity
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

% Report orphan EMG files without matching EEG
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
% Utility: create a single entry struct with consistent field order.
e = struct('pairName',pairName,'phase',phase,'mffPath',mffPath,'otbPath',otbPath, ...
    'state',state,'msg',msg,'layout',layoutLabel);
end

function [mffPaths, otbPaths] = listBidsFiles(rootPath)
% Recursively list .mff directories and .otb+ files under rootPath.
mffList = dir(fullfile(rootPath,'**','*.mff'));
mffList = mffList([mffList.isdir]);  % .mff is a directory
otbList = dir(fullfile(rootPath,'**','*.otb+'));
otbList = otbList(~[otbList.isdir]);  % .otb+ is a tar file

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
% Match EEG (.mff) to EMG (.otb+) files by name similarity.
% Uses greedy one-to-one assignment: each EMG is used at most once.
% Returns matchIdx(i) = j (1-based) if mffPaths{i} matches otbPaths{j},
% or 0 if no match.
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
% Compute a similarity score between two file paths based on shared
% meaningful tokens (after stripping stop words, numbers, sides).
% Side agreement adds +2; side mismatch subtracts 2.
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
% Extract meaningful tokens from a file path for name matching.
% Strips: extension, stop words, bare numbers, alphanumeric codes.
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
% Normalize synonym tokens to a canonical form.
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
% Return 'dx', 'sx', or '' if a side token is present in the token list.
if any(strcmp(toks,'dx'))
    s = 'dx';
elseif any(strcmp(toks,'sx'))
    s = 'sx';
else
    s = '';
end
end

function s = sanitizeName(s)
% Clean a string for use as a filename or identifier:
% lowercase, collapse non-alphanumeric to underscores.
s = lower(char(s));
s = regexprep(s, '[^a-z0-9]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'run'; end
end

function phase = inferPhaseLabel(sesName, sesPath, mffPath, otbPath)
% Infer whether this session is 'pre', 'post', or other from its name/path.
% Checks for 'pre' or 'post' substrings in session name or path.
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

function ref = readMffReference(mffPath, eventFile)
% Read EEG reference metadata and DIN trigger events from a .mff folder.
%
% Parses:
%   info.xml    -> recording start time and timezone
%   epochs.xml  -> total recording duration
%   Events_*.xml -> DIN events (DIN1, DIN4, DIN5) with absolute timestamps
%
% Returns a struct with recording timing and DIN trigger times (in seconds
% relative to recording start).
infoTxt = fileread(fullfile(mffPath,'info.xml'));
epochsTxt = fileread(fullfile(mffPath,'epochs.xml'));
evPath = fullfile(mffPath,eventFile);
evTxt = fileread(evPath);

recordIso = oneToken(infoTxt,'<recordTime>([^<]+)</recordTime>');
recordEpoch = isoToEpoch(recordIso);
offTok = regexp(recordIso,'([+-]\d\d:\d\d)$','tokens','once');
if isempty(offTok), error('Cannot parse timezone offset in %s', recordIso); end
offsetStr = offTok{1};

endUs = str2double(oneToken(epochsTxt,'<endTime>(\d+)</endTime>'));
durationS = endUs / 1e6;

name = oneToken(evTxt,'<name>([^<]+)</name>');
trackType = oneToken(evTxt,'<trackType>([^<]+)</trackType>');

% Parse all events from the event XML
allEv = regexp(evTxt,'<event>\s*<beginTime>([^<]+)</beginTime>\s*<duration>([^<]+)</duration>\s*<code>([^<]+)</code>[\s\S]*?</event>','tokens');
if isempty(allEv), error('No events in %s', evPath); end

templateDuration = str2double(allEv{1}{2});
codeTimes = struct('DIN1',[],'DIN4',[],'DIN5',[]);
for i = 1:numel(allEv)
    tIso = allEv{i}{1};
    code = allEv{i}{3};
    t = isoToEpoch(tIso) - recordEpoch;
    if isfield(codeTimes,code)
        codeTimes.(code)(end+1) = t; %#ok<AGROW>
    end
end

ref = struct('mffPath',mffPath,'eventPath',evPath,'recordEpoch',recordEpoch, ...
    'offsetStr',offsetStr,'durationS',durationS,'name',name,'trackType',trackType, ...
    'codeTimes',codeTimes,'templateDuration',templateDuration);
end

function otb = readOtb(otbPath)
% Read EMG data and AUX trigger channel mapping from an .otb+ file.
%
% .otb+ is a tar archive containing:
%   - a .xml config file (contains channel metadata, sample rate, AUX mapping)
%   - a .sig binary file (raw int16 samples, all channels interleaved)
%
% Returns a struct with:
%   fs        -- sampling frequency (Hz)
%   nChannels -- total number of channels
%   data      -- [nChannels x N] matrix of raw EMG samples
%   auxMap    -- containers.Map: AUX channel number -> row index in data
tmpDir = tempname; mkdir(tmpDir);
cleanup = onCleanup(@() cleanupTemp(tmpDir)); %#ok<NASGU>
untar(otbPath, tmpDir);

% Find the XML metadata file (prefer *_08.xml pattern)
xmlList = dir(fullfile(tmpDir,'*_08.xml'));
if isempty(xmlList), xmlList = dir(fullfile(tmpDir,'*.xml')); end
if isempty(xmlList), error('No xml in %s', otbPath); end
xmlPath = fullfile(xmlList(1).folder, xmlList(1).name);

% Find the binary signal file
sigList = dir(fullfile(tmpDir,'*.sig'));
if isempty(sigList), error('No sig in %s', otbPath); end
sigPath = fullfile(sigList(1).folder, sigList(1).name);

% Parse XML to get sampling rate, channel count, and AUX trigger mapping
doc = xmlread(xmlPath);
root = doc.getDocumentElement();
fs = str2double(char(root.getAttribute('SampleFrequency')));
nChannels = str2double(char(root.getAttribute('DeviceTotalChannels')));

% Build a map from AUX channel number (1-based) to data row index (0-based)
% Only AUX channels with 'Trigger' in their ID string are mapped.
auxMap = containers.Map('KeyType','double','ValueType','double');
adapters = root.getElementsByTagName('Adapter');
for i = 0:(adapters.getLength()-1)
    ad = adapters.item(i);
    startIdx = str2double(char(ad.getAttribute('ChannelStartIndex')));
    chNodes = ad.getElementsByTagName('Channel');
    for j = 0:(chNodes.getLength()-1)
        ch = chNodes.item(j);
        prefix = char(ch.getAttribute('Prefix'));
        idStr = char(ch.getAttribute('ID'));
        idx = str2double(char(ch.getAttribute('Index')));
        if contains(prefix,'AUX') && contains(idStr,'Trigger')
            tk = regexp(prefix,'AUX\s+(\d+)','tokens','once');
            if ~isempty(tk)
                auxMap(str2double(tk{1})) = startIdx + idx;
            end
        end
    end
end

% Read raw int16 samples (little-endian) and reshape into [nChannels x N]
fid = fopen(sigPath,'r');
raw = fread(fid,inf,'int16=>double',0,'l');
fclose(fid);
N = floor(numel(raw)/nChannels);
raw = raw(1:N*nChannels);
data = reshape(raw,[nChannels,N]);

otb = struct('fs',fs,'nChannels',nChannels,'data',data,'auxMap',auxMap);
end

function [events, meta] = buildAlignedEvents(ref, otb, templateMffPath)
% Core alignment: maps OTB (AUX) trigger falling edges to MFF (DIN) events.
%
% Trigger mapping (hardware configuration):
%   AUX4 (OTB, falling edge) <-> DIN4 (EEG) -> keyboard event
%   AUX3 (OTB, falling edge) <-> DIN5 (EEG) -> beep event
%   AUX2 (OTB, falling edge) <-> DIN1 (EEG) -> sync wave (optional refinement)
%
% Algorithm steps:
%   1) Detect falling edges on AUX channels
%   2) Clean trigger edges (remove noise spikes via minimum-gap filter)
%   3) Plausibility check on cleaned trigger counts
%   4) Affine fit between AUX and DIN times (with lag search)
%   5) Optional DIN1 refinement to improve clock drift estimate
%   6) Guardrail: reject unstable fits (|a-1| > 2%, MAD > 500 ms)
%   7) Map all AUX edges to MFF timebase using the affine transform
%   8) Reconstruct task events (start, keyboard, beep, end)
%
% Returns:
%   events -- struct with .times (s) and .codes (cell) of all aligned events
%   meta   -- struct with fit parameters, cleaning stats, logical metadata
%
% Required AUX channels in OTB: 2 (DIN1), 3 (DIN5), 4 (DIN4)
needed = [2 3 4];
for k = needed
    if ~isKey(otb.auxMap,k), error('OTB missing AUX%d', k); end
end

din4 = ref.codeTimes.DIN4;  % keyboard events from EEG
din5 = ref.codeTimes.DIN5;  % beep events from EEG
din1 = ref.codeTimes.DIN1;  % sync wave from EEG
if isempty(din4) || isempty(din5), error('Reference MFF must contain DIN4 and DIN5'); end

% Step 1-2: Detect and clean falling edges on each AUX channel
% AUX4 (DIN4) and AUX3 (DIN5) use the corresponding DIN count as expected
% reference; AUX2 (DIN1) uses an adaptive gap based on the DIN1 step size.
edges = struct();
rawCounts = struct();
cleanCounts = struct();
cleaning = struct();
for k = needed
    x = otb.data(otb.auxMap(k)+1,:);
    tRaw = detectFallingEdges(x, otb.fs);
    rawCounts.(sprintf('AUX%d',k)) = numel(tRaw);

    if k == 4
        [tClean, cInfo] = cleanTriggerEdges(tRaw, numel(din4), 0.12);
    elseif k == 3
        [tClean, cInfo] = cleanTriggerEdges(tRaw, numel(din5), 0.12);
    else
        % AUX2: adaptive gap based on DIN1 density
        din1Step = NaN;
        if numel(din1) >= 2
            din1Step = median(diff(din1));
        end
        if ~isfinite(din1Step) || din1Step <= 0
            minGapAux2 = 0.008;
        else
            minGapAux2 = min(0.05, max(0.004, 0.25 * din1Step));
        end
        [tClean, cInfo] = cleanTriggerEdges(tRaw, numel(din1), minGapAux2);
    end

    edges.(sprintf('AUX%d',k)) = tClean;
    cleanCounts.(sprintf('AUX%d',k)) = numel(tClean);
    cleaning.(sprintf('AUX%d',k)) = cInfo;
end

% Step 3: Verify that cleaned trigger counts are in a plausible range
% (60%-180% of expected DIN count)
checkTriggerCountPlausibility(cleanCounts.AUX4, numel(din4), 'AUX4', 'DIN4');
checkTriggerCountPlausibility(cleanCounts.AUX3, numel(din5), 'AUX3', 'DIN5');

aux4 = edges.AUX4;
aux3 = edges.AUX3;
aux2 = edges.AUX2;

% Step 4: Affine fit between AUX and DIN times.
% Try both mapping orders (AUX4<->DIN4/AUX3<->DIN5 and the reverse) and
% select the one with lower MAD.
% NOTE: The hardware mapping is fixed (AUX4->DIN4, AUX3->DIN5), but the
% dual-fit approach handles possible cable swaps.
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

keyO = fitSel.x4;    % AUX4 times (OTB timebase)
keyM = fitSel.y4;    % paired DIN4 times (MFF timebase)
beepO = fitSel.x3;   % AUX3 times (OTB timebase)
beepM = fitSel.y3;   % paired DIN5 times (MFF timebase)
nCommon = min(numel(keyO), numel(beepO));

% Refit the affine transform using only trials that satisfy the
% keyboard-before-beep logical order (trial gating on DIN4/DIN5).
% Also try an endpoint-based fit as a robust alternative.
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
    tKeyPred = a*keyO + b;    % predicted keyboard times in MFF base
    tBeepPred = a*beepO + b;  % predicted beep times in MFF base
    % Trial-gating uses DIN4/DIN5 order only (keyboard must precede beep).
    idxKeep = findValidTrialIndicesKbBeep(tKeyPred, tBeepPred);
    fitInfo.n_fit_trials_kept = numel(idxKeep);
    if numel(idxKeep) >= 2
        xFit = [keyO(idxKeep) beepO(idxKeep)];
        yFit = [keyM(idxKeep) beepM(idxKeep)];
        [aV,bV,mV] = fitLinear(xFit, yFit);

        % Alternative: endpoint-based fit (first/last point slope)
        % as a robust estimator; use it if it gives equal or better MAD.
        [xFitS,ordV] = sort(xFit);
        yFitS = yFit(ordV);
        if numel(xFitS)>=2 && xFitS(end)>xFitS(1)
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

% Base point set for optional DIN1 refinement
if ~isempty(idxKeep) && numel(idxKeep) >= 2
    xBase = [keyO(idxKeep) beepO(idxKeep)];
    yBase = [keyM(idxKeep) beepM(idxKeep)];
else
    xBase = [keyO beepO];
    yBase = [keyM beepM];
end

% Step 5: Optional Sync refinement using the dense DIN1 train.
% DIN1 (AUX2) provides many additional point pairs that can improve
% the clock-drift estimate without affecting trial validity logic
% (which remains DIN4/DIN5-only).
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

% Step 6: Guardrail -- reject clearly unstable fits.
% Criteria: non-finite parameters, |a-1| > 2% (excessive clock drift),
% or MAD > 500 ms (unreliable alignment).
if ~isfinite(a) || ~isfinite(b) || ~isfinite(mad) || abs(a - 1) > 0.02 || mad > 0.5
    error('Unstable affine fit after trigger cleaning: a=%.6f, b=%.6f, MAD=%.3f s', a, b, mad);
end

% Step 7: Map all AUX falling edges to MFF timebase using the affine
% transform t_mff = a * t_otb + b. Discard events outside recording bounds.
times = []; codes = {};
for t = aux4
    tt = a*t+b; if tt>=0 && tt<=ref.durationS+0.1, times(end+1)=tt; codes{end+1}=map4; end %#ok<AGROW>
end
for t = aux3
    tt = a*t+b; if tt>=0 && tt<=ref.durationS+0.1, times(end+1)=tt; codes{end+1}=map3; end %#ok<AGROW>
end
for t = aux2
    tt = a*t+b; if tt>=0 && tt<=ref.durationS+0.1, times(end+1)=tt; codes{end+1}='DIN1'; end %#ok<AGROW>
end
[times,ord2] = sort(times); codes = codes(ord2);

% Step 8: Reconstruct task events (start, keyboard, beep, end) from DIN4/DIN5
% using temporal offsets derived from a template recording.
events = struct('times',times,'codes',{codes});
[events, logicalMeta] = addReconstructedTaskEvents(events, ref.durationS, templateMffPath);

meta = struct('a',a,'b',b,'mad',mad,'driftPpm',(a-1)*1e6, ...
    'mapping',struct('AUX4',map4,'AUX3',map3,'AUX2','DIN1'), ...
    'rawCounts',rawCounts, ...
    'cleanCounts',cleanCounts, ...
    'cleaning',cleaning, ...
    'fit',fitInfo, ...
    'logical',logicalMeta, ...
    'finalCounts',countCodes(events.codes));
end

function [outOtbPath, segCsvPath, info] = exportMovementOnlyOtb(otbSrc, pairOut, meta, otb)
% Create a new .otb+ file containing only movement segments extracted from
% valid trial windows.
%
% The movement windows are built from reconstructed task events:
%   primary rule:  beep -> start
%   fallback rule: beep -> end (when start is not after beep)
%
% Returns paths to the movement-only OTB and the segment CSV, plus an info
% struct with statistics.
[~, otbBase, ~] = fileparts(otbSrc);
outOtbPath = fullfile(pairOut, [otbBase '_movement_only.otb+']);
segCsvPath = fullfile(pairOut, [otbBase '_movement_segments.csv']);

info = struct('enabled',false, ...
    'rule_requested','beep_to_start', ...
    'rule_used','', ...
    'reason','', ...
    'n_segments',0, ...
    'n_samples',0, ...
    'duration_s',0);

% Clean up any pre-existing output files
if exist(outOtbPath,'file')
    delete(outOtbPath);
end
legacyTarPath = [outOtbPath '.tar'];
if exist(legacyTarPath,'file')
    delete(legacyTarPath);
end
if exist(segCsvPath,'file')
    delete(segCsvPath);
end

% Determine sample ranges for movement windows
[sampleRanges, segRows, ruleUsed, reason] = buildMovementWindows(meta, otb);
info.rule_used = ruleUsed;
info.reason = reason;
if isempty(sampleRanges)
    return;
end

% Concatenate the selected sample ranges and write to new OTB
moveData = concatenateSampleRanges(otb.data, sampleRanges);
writeMovementOtbFromSource(otbSrc, outOtbPath, moveData);
writeMovementSegmentsCsv(segCsvPath, segRows);

info.enabled = true;
info.n_segments = size(sampleRanges,1);
info.n_samples = size(moveData,2);
info.duration_s = size(moveData,2) / otb.fs;
end

function [sampleRanges, segRows, ruleUsed, reason] = buildMovementWindows(meta, otb)
% Determine sample ranges in OTB timebase corresponding to movement windows.
%
% Movement windows are defined by trial events:
%   Primary: beep -> start  (start is the onset of movement)
%   Fallback: beep -> end   (when start does not occur after beep)
%
% Returns:
%   sampleRanges -- [N x 2] matrix with [start_sample, end_sample] per segment
%   segRows      -- cell array with metadata rows for the segments CSV
%   ruleUsed     -- which windowing rule was applied
%   reason       -- description if no windows could be built
sampleRanges = zeros(0,2);
segRows = {};
ruleUsed = '';
reason = '';

% Validate required metadata
if ~isfield(meta,'logical') || ~isstruct(meta.logical)
    reason = 'logical_metadata_missing';
    return;
end
lm = meta.logical;
if ~isfield(lm,'enabled') || ~lm.enabled
    reason = 'logical_reconstruction_disabled';
    return;
end
if ~isfield(lm,'trial_times_s') || ~isstruct(lm.trial_times_s)
    reason = 'trial_times_missing';
    return;
end
if ~isfield(meta,'a') || ~isfield(meta,'b') || ~isfinite(meta.a) || ~isfinite(meta.b) || abs(meta.a) < eps
    reason = 'invalid_time_transform';
    return;
end

tt = lm.trial_times_s;
if ~isfield(tt,'beep')
    reason = 'beep_times_missing';
    return;
end

% Try primary rule: beep -> start
tBeepAll = double(tt.beep(:)');
tBeep = tBeepAll;
tStop = [];
ruleUsed = 'beep_to_start';

if isfield(tt,'start')
    tStart = double(tt.start(:)');
    n = min(numel(tBeepAll), numel(tStart));
    tBeep = tBeepAll(1:n);
    tStop = tStart(1:n);
end
if isempty(tStop)
    tStop = nan(size(tBeep));
end

valid = isfinite(tBeep) & isfinite(tStop) & (tStop > tBeep);

% In the reconstructed event order, 'start' typically precedes 'beep'.
% If beep->start windows are invalid (start occurs before beep),
% fall back to beep->end to preserve movement segments.
if ~any(valid)
    if ~isfield(tt,'end')
        reason = 'beep_to_start_invalid_and_end_missing';
        return;
    end
    tEnd = double(tt.end(:)');
    n = min(numel(tBeepAll), numel(tEnd));
    tBeep = tBeepAll(1:n);
    tStop = tEnd(1:n);
    valid = isfinite(tBeep) & isfinite(tStop) & (tStop > tBeep);
    ruleUsed = 'beep_to_end_fallback';
    if ~any(valid)
        reason = 'no_valid_movement_windows_after_fallback';
        return;
    end
    reason = 'beep_to_start_not_valid_in_reconstructed_sequence';
end

% Filter to valid trials and map to OTB sample indices
trialIdx = 1:numel(tBeep);
if isfield(lm,'trial_indices')
    tri = double(lm.trial_indices(:)');
    if numel(tri) >= numel(trialIdx)
        trialIdx = tri(1:numel(trialIdx));
    end
end

tBeep = tBeep(valid);
tStop = tStop(valid);
trialIdx = trialIdx(valid);

% Convert MFF timestamps to OTB sample indices via inverse affine transform
nTot = size(otb.data,2);
tBeepOtb = (tBeep - meta.b) / meta.a;
tStopOtb = (tStop - meta.b) / meta.a;
i1 = floor(tBeepOtb * otb.fs) + 1;
i2 = ceil(tStopOtb * otb.fs);
i1 = max(i1, 1);
i2 = min(i2, nTot);
ok = (i2 > i1);
if ~any(ok)
    sampleRanges = zeros(0,2);
    segRows = {};
    if isempty(reason), reason = 'valid_time_windows_outside_otb_bounds'; end
    return;
end

i1 = i1(ok);
i2 = i2(ok);
tBeep = tBeep(ok);
tStop = tStop(ok);
trialIdx = trialIdx(ok);
[i1, ord] = sort(i1);
i2 = i2(ord);
tBeep = tBeep(ord);
tStop = tStop(ord);
trialIdx = trialIdx(ord);

% Build output arrays
sampleRanges = [i1(:) i2(:)];
segRows = cell(size(sampleRanges,1), 11);
for k = 1:size(sampleRanges,1)
    nS = sampleRanges(k,2) - sampleRanges(k,1) + 1;
    segRows{k,1} = sprintf('%d', k);
    segRows{k,2} = sprintf('%d', round(trialIdx(k)));
    segRows{k,3} = ruleUsed;
    segRows{k,4} = sprintf('%.6f', tBeep(k));
    segRows{k,5} = sprintf('%.6f', tStop(k));
    segRows{k,6} = sprintf('%.6f', (sampleRanges(k,1)-1)/otb.fs);
    segRows{k,7} = sprintf('%.6f', sampleRanges(k,2)/otb.fs);
    segRows{k,8} = sprintf('%d', sampleRanges(k,1));
    segRows{k,9} = sprintf('%d', sampleRanges(k,2));
    segRows{k,10} = sprintf('%d', nS);
    segRows{k,11} = sprintf('%.6f', nS/otb.fs);
end
end

function dataOut = concatenateSampleRanges(dataIn, ranges)
% Concatenate selected sample ranges from a multichannel data matrix.
% dataIn: [nChannels x N]  ranges: [nSeg x 2] with [start, end] indices
nChan = size(dataIn,1);
nSeg = size(ranges,1);
len = ranges(:,2) - ranges(:,1) + 1;
tot = sum(len);
dataOut = zeros(nChan, tot);
p = 1;
for k = 1:nSeg
    idx = ranges(k,1):ranges(k,2);
    q = p + numel(idx) - 1;
    dataOut(:,p:q) = dataIn(:,idx);
    p = q + 1;
end
end

function writeMovementSegmentsCsv(pathCsv, rows)
% Write movement segment metadata to CSV with predefined header.
header = {'segment_id','trial_index','rule_used','mff_t_start_s','mff_t_stop_s', ...
    'otb_t_start_s','otb_t_stop_s','sample_start','sample_stop','n_samples','duration_s'};
writeCsv(pathCsv, header, rows);
end

function writeMovementOtbFromSource(srcOtbPath, outOtbPath, data)
% Create a new .otb+ file by replacing the .sig binary in a copy of the
% source OTB archive with the movement-only data.
%
% The input .otb+ is a tar archive; we untar it, overwrite the .sig file,
% and re-tar to produce the output.
tmpDir = tempname;
mkdir(tmpDir);
cleanup = onCleanup(@() cleanupTemp(tmpDir)); %#ok<NASGU>
untar(srcOtbPath, tmpDir);

sigList = dir(fullfile(tmpDir,'*.sig'));
if isempty(sigList)
    error('Cannot build movement OTB: no .sig file in %s', srcOtbPath);
end
if numel(sigList) > 1
    error('Cannot build movement OTB: multiple .sig files in %s', srcOtbPath);
end
sigPath = fullfile(sigList(1).folder, sigList(1).name);

% Clamp to int16 range and write
x = round(double(data(:)));
x(x > 32767) = 32767;
x(x < -32768) = -32768;

fid = fopen(sigPath, 'w');
if fid < 0
    error('Cannot open .sig for writing: %s', sigPath);
end
fwrite(fid, int16(x), 'int16', 0, 'l');
fclose(fid);

% Re-tar the archive with the updated .sig
if exist(outOtbPath,'file')
    delete(outOtbPath);
end
pack = dir(tmpDir);
pack = pack(~ismember({pack.name},{'.','..'}));
packNames = {pack.name};
tmpTarPath = [tempname '.tar'];
tar(tmpTarPath, packNames, tmpDir);
movefile(tmpTarPath, outOtbPath, 'f');
end

function writeEventXml(ref, events)
% Write aligned events into the MFF event XML file.
% Creates a backup of the original event XML before overwriting.
%
% Each event includes ISO-8601 timestamps, a gidx (global index), and a
% cidx (per-code counter) as custom keys.
backup = strrep(ref.eventPath,'.xml','_original.xml');
if ~exist(backup,'file'), copyfile(ref.eventPath, backup); end

cidx = containers.Map('KeyType','char','ValueType','double');
L = {};
L{end+1} = '<?xml version="1.0" encoding="UTF-8"?>'; %#ok<AGROW>
L{end+1} = '<eventTrack xmlns="http://www.egi.com/event_mff" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'; %#ok<AGROW>
L{end+1} = sprintf('    <name>%s</name>', ref.name); %#ok<AGROW>
L{end+1} = sprintf('    <trackType>%s</trackType>', ref.trackType); %#ok<AGROW>

for i = 1:numel(events.times)
    code = events.codes{i};
    if ~isKey(cidx, code), cidx(code) = 0; end
    cidx(code) = cidx(code) + 1;
    beginIso = epochToIso(ref.recordEpoch + events.times(i), ref.offsetStr);
    L{end+1} = '    <event>'; %#ok<AGROW>
    L{end+1} = sprintf('        <beginTime>%s</beginTime>', beginIso); %#ok<AGROW>
    L{end+1} = sprintf('        <duration>%d</duration>', ref.templateDuration); %#ok<AGROW>
    L{end+1} = sprintf('        <code>%s</code>', code); %#ok<AGROW>
    L{end+1} = sprintf('        <label>%s</label>', code); %#ok<AGROW>
    L{end+1} = sprintf('        <description>%d</description>', i); %#ok<AGROW>
    L{end+1} = '        <sourceDevice></sourceDevice>'; %#ok<AGROW>
    L{end+1} = '        <keys>'; %#ok<AGROW>
    L{end+1} = '            <key>'; %#ok<AGROW>
    L{end+1} = '                <keyCode>gidx</keyCode>'; %#ok<AGROW>
    L{end+1} = sprintf('                <data dataType="string">%d</data>', i); %#ok<AGROW>
    L{end+1} = '            </key>'; %#ok<AGROW>
    L{end+1} = '            <key>'; %#ok<AGROW>
    L{end+1} = '                <keyCode>cidx</keyCode>'; %#ok<AGROW>
    L{end+1} = sprintf('                <data dataType="string">%d</data>', cidx(code)); %#ok<AGROW>
    L{end+1} = '            </key>'; %#ok<AGROW>
    L{end+1} = '        </keys>'; %#ok<AGROW>
    L{end+1} = '    </event>'; %#ok<AGROW>
end
L{end+1} = '</eventTrack>'; %#ok<AGROW>

fid = fopen(ref.eventPath,'w');
for i = 1:numel(L), fprintf(fid,'%s\n',L{i}); end
fclose(fid);
end

function writeEventsCsv(pathCsv, events)
% Write aligned events to CSV with columns:
% Event_type, Event_name, Event_time, Event_duration, Offset.
% Original triggers (DIN1/4/5) are labeled 'OTB_trigger';
% reconstructed events (start/keyboard/beep/end) are 'OTB_reconstructed'.
fid = fopen(pathCsv,'w');
fprintf(fid,'Event_type,Event_name,Event_time,Event_duration,Offset\n');
for i = 1:numel(events.times)
    code = events.codes{i};
    if strcmp(code,'DIN1') || strcmp(code,'DIN4') || strcmp(code,'DIN5')
        et = 'OTB_trigger';
    else
        et = 'OTB_reconstructed';
    end
    fprintf(fid,'%s,%s,%.6f,0.001,0\n', et, code, events.times(i));
end
fclose(fid);
end

function [eventsOut, info] = addReconstructedTaskEvents(eventsIn, durationS, templateMffPath)
% Reconstruct task-level events (start, keyboard, beep, end) from DIN4/DIN5.
%
% The temporal offsets between DIN4/DIN5 and the task events are derived
% from a template MFF recording that contains pre-labeled events.
%
% Trial validity is enforced through the logical sequence:
%   start < keyboard < beep < end < next_keyboard
%
% DIN4 maps to 'keyboard', DIN5 maps to 'beep'. The 'start' and 'end'
% events are placed relative to these anchors using the template offsets.
%
% Returns updated events struct (with added task events) and an info struct
% containing trial filtering statistics and per-trial event times.
times = eventsIn.times(:)';
codes = eventsIn.codes;

din4 = times(strcmp(codes,'DIN4'));  % keyboard anchors
din5 = times(strcmp(codes,'DIN5'));  % beep anchors
nTrials = min(numel(din4), numel(din5));
if nTrials == 0
    eventsOut = eventsIn;
    info = struct('enabled',false,'reason','no_din4_or_din5');
    return;
end

% Filter trials: keep only valid DIN4 < DIN5 sequences
idxKeep = findValidTrialIndicesKbBeep(din4(1:nTrials), din5(1:nTrials));
tKeyKeep = din4(idxKeep);
tBeepKeep = din5(idxKeep);

nKeep = numel(tKeyKeep);
if nKeep == 0
    eventsOut = eventsIn;
    info = struct('enabled',false,'reason','no_valid_trials_with_keyboard_beep_sequence', ...
        'n_trials_candidate',nTrials,'n_trials_removed',nTrials,'n_trials',0);
    return;
end

% Apply temporal offsets from template to derive start/keyboard/beep/end
[offStart, offKey, offBeep, offEnd, offSource] = getOffsetsProfile(templateMffPath, nKeep);
tStart = tKeyKeep + offStart;
tKey = tKeyKeep + offKey;
tBeep = tBeepKeep + offBeep;
tEnd = tBeepKeep + offEnd;

% Enforce logical sequence and remove invalid trials.
% Required order: start < keyboard < beep < end < next_keyboard
trialIdx = idxKeep;
nextKeyBound = inf(1, numel(trialIdx));
hasNext = trialIdx < nTrials;
nextKeyBound(hasNext) = din4(trialIdx(hasNext) + 1);

validSeq = isfinite(tStart) & isfinite(tKey) & isfinite(tBeep) & isfinite(tEnd) & ...
    (tStart < tKey) & (tKey < tBeep) & (tBeep < tEnd) & (tEnd < nextKeyBound);
if any(~validSeq)
    tStart = tStart(validSeq);
    tKey = tKey(validSeq);
    tBeep = tBeep(validSeq);
    tEnd = tEnd(validSeq);
    offStart = offStart(validSeq);
    offKey = offKey(validSeq);
    offBeep = offBeep(validSeq);
    offEnd = offEnd(validSeq);
    trialIdx = trialIdx(validSeq);
    nextKeyBound = nextKeyBound(validSeq);
end

% Clamp event times to recording bounds
tStart = clampTimes(tStart, durationS);
tKey = clampTimes(tKey, durationS);
tBeep = clampTimes(tBeep, durationS);
tEnd = clampTimes(tEnd, durationS);

% Re-check sequence after clamping
validAfterClamp = (tStart < tKey) & (tKey < tBeep) & (tBeep < tEnd) & (tEnd < nextKeyBound);
if any(~validAfterClamp)
    tStart = tStart(validAfterClamp);
    tKey = tKey(validAfterClamp);
    tBeep = tBeep(validAfterClamp);
    tEnd = tEnd(validAfterClamp);
    offStart = offStart(validAfterClamp);
    offKey = offKey(validAfterClamp);
    offBeep = offBeep(validAfterClamp);
    offEnd = offEnd(validAfterClamp);
end

nKeepFinal = numel(tKey);
if nKeepFinal == 0
    eventsOut = eventsIn;
    info = struct('enabled',false,'reason','no_valid_trials_after_start_keyboard_beep_end_reconstruction', ...
        'n_trials_candidate',nTrials,'n_trials_removed',nTrials,'n_trials',0);
    return;
end

% Merge reconstructed events with original DIN events and sort
addTimes = [tStart tKey tBeep tEnd];
addCodes = [repmat({'start'},1,nKeepFinal) repmat({'keyboard'},1,nKeepFinal) ...
            repmat({'beep'},1,nKeepFinal) repmat({'end'},1,nKeepFinal)];

allTimes = [times addTimes];
allCodes = [codes addCodes];
[allTimes, ord] = sort(allTimes);
allCodes = allCodes(ord);

eventsOut = struct('times',allTimes,'codes',{allCodes});
info = struct();
info.enabled = true;
info.n_trials_candidate = nTrials;
info.n_trials = nKeepFinal;
info.n_trials_removed = nTrials - nKeepFinal;
info.rule = 'valid trial if DIN4(keyboard) < DIN5(beep) and reconstructed order is start < keyboard < beep < end with end before next keyboard';
info.offset_source = offSource;
info.offset_median_ms = struct( ...
    'start',median(offStart)*1000, ...
    'keyboard',median(offKey)*1000, ...
    'beep',median(offBeep)*1000, ...
    'end',median(offEnd)*1000);
info.trial_indices = trialIdx;
info.trial_times_s = struct( ...
    'start', tStart, ...
    'keyboard', tKey, ...
    'beep', tBeep, ...
    'end', tEnd);
end

function idxKeep = findValidTrialIndicesKbBeep(din4, din5)
% Find valid trial pairs (DIN4, DIN5) that satisfy:
%   - keyboard (DIN4) precedes beep (DIN5) by at least 50 ms
%   - beep (DIN5) precedes the next keyboard (DIN4) if one exists
% This implements the trial-gating rule: only DIN4/DIN5 order matters.
nTrials = min(numel(din4), numel(din5));
idxKeep = [];
minKbBeepGapS = 0.05; % minimum gap to reject near-simultaneous pairs
if nTrials == 0
    return;
end

for k = 1:nTrials
    tKey = din4(k);
    tBeep = din5(k);
    if ~(isfinite(tKey) && isfinite(tBeep) && ((tKey + minKbBeepGapS) < tBeep))
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
% Search for the best lag combination (lag4, lag3) that minimizes MAD of
% the affine fit between AUX and DIN trigger pairs.
%
% aux4/aux3: OTB falling edge times
% din4Target/din5Target: expected EEG DIN times to match against
%
% Returns a fit struct with affine parameters (a, b), MAD, the optimal
% lags, and the paired (x, y) values used in the fit.
if nargin < 5 || isempty(maxLag)
    maxLag = 8;
end
fit = struct('valid',false,'a',NaN,'b',NaN,'mad',Inf, ...
    'lag4',0,'lag3',0,'x4',[],'y4',[],'x3',[],'y3',[]);

bestScore = Inf;
for lag4 = -maxLag:maxLag
    [x4, y4] = pairWithLag(aux4, din4Target, lag4);
    if numel(x4) < 10  % need at least 10 pairs for a reliable fit
        continue;
    end
    for lag3 = -maxLag:maxLag
        [x3, y3] = pairWithLag(aux3, din5Target, lag3);
        if numel(x3) < 10
            continue;
        end
        x = [x4 x3];
        y = [y4 y3];
        [aT, bT, mT] = fitLinear(x, y);
        sc = mT + 1e-9*abs(aT-1);  % MAD + tiny regularizer on drift
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

function [x, y] = pairWithLag(aux, din, lag)
% Pair AUX edges with DIN events at a given lag offset.
% Positive lag: skip the first `lag` AUX edges before pairing.
% Negative lag: skip the first `-lag` DIN events before pairing.
x = [];
y = [];
if isempty(aux) || isempty(din)
    return;
end

if lag >= 0
    n = min(numel(aux)-lag, numel(din));
    if n <= 0
        return;
    end
    x = aux((1+lag):(lag+n));
    y = din(1:n);
else
    n = min(numel(aux), numel(din)+lag);
    if n <= 0
        return;
    end
    k = -lag;
    x = aux(1:n);
    y = din((1+k):(k+n));
end
end

function [offStart, offKey, offBeep, offEnd, source] = getOffsetsProfile(templateMffPath, nTrials)
% Get temporal offsets to reconstruct task events from DIN4/DIN5.
%
% Offsets define where 'start', 'keyboard', 'beep', and 'end' events
% should be placed relative to DIN4 (keyboard) and DIN5 (beep).
%
% If a template MFF with pre-labeled events is available, offsets are
% derived from it; otherwise hard-coded defaults are used.
% Fallback profile based on observed template behavior:
% start before DIN4, keyboard approx DIN4, beep approx DIN5, end after DIN5.
baseStart = -2.2;
baseKey = 0;
baseBeep = 0;
baseEnd = 0.64;

source = 'defaults';
offStartT = [];
offKeyT = [];
offBeepT = [];
offEndT = [];

if ~isempty(templateMffPath) && isfolder(templateMffPath)
    try
        [offStartT, offKeyT, offBeepT, offEndT] = deriveOffsetsFromTemplate(templateMffPath);
        if ~isempty(offStartT) && ~isempty(offKeyT) && ~isempty(offBeepT) && ~isempty(offEndT)
            source = 'template_mff';
        end
    catch
        % keep defaults
    end
end

if isempty(offStartT), offStartT = baseStart; end
if isempty(offKeyT), offKeyT = baseKey; end
if isempty(offBeepT), offBeepT = baseBeep; end
if isempty(offEndT), offEndT = baseEnd; end

offStart = resampleOffsets(offStartT, nTrials);
offKey = resampleOffsets(offKeyT, nTrials);
offBeep = resampleOffsets(offBeepT, nTrials);
offEnd = resampleOffsets(offEndT, nTrials);
end

function [offStart, offKey, offBeep, offEnd] = deriveOffsetsFromTemplate(templateMffPath)
% Derive temporal offsets from a template MFF recording that contains
% pre-labeled events (start, keyboard, beep, end) alongside DIN triggers.
%
% The offsets are computed as:
%   start_offset = start_time - DIN4_time
%   keyboard_offset = keyboard_time - DIN4_time
%   beep_offset = beep_time - DIN5_time
%   end_offset = end_time - DIN5_time
%
% Warm-up DIN pulses are handled by aligning task labels to tail DIN sequences.
offStart = [];
offKey = [];
offBeep = [];
offEnd = [];

evXml = '';
cand = dir(fullfile(templateMffPath,'Events*.xml'));
if isempty(cand), return; end

% Prefer EEGLAB-exported track when present (contains logical labels).
for i = 1:numel(cand)
    if contains(lower(cand(i).name), 'eeglab')
        evXml = fullfile(cand(i).folder, cand(i).name);
        break;
    end
end
if isempty(evXml), evXml = fullfile(cand(1).folder, cand(1).name); end

txt = fileread(evXml);
allEv = regexp(txt,'<event>\s*<beginTime>([^<]+)</beginTime>[\s\S]*?<code>([^<]+)</code>[\s\S]*?</event>','tokens');
if isempty(allEv), return; end

t = zeros(1,numel(allEv));
c = cell(1,numel(allEv));
for i = 1:numel(allEv)
    t(i) = isoToEpoch(allEv{i}{1});
    c{i} = allEv{i}{2};
end
t = t - min(t);

din4 = t(strcmp(c,'DIN4'));
din5 = t(strcmp(c,'DIN5'));
start = t(strcmp(c,'start'));
key = t(strcmp(c,'keyboard'));
beep = t(strcmp(c,'beep'));
en = t(strcmp(c,'end'));
if isempty(din4) || isempty(din5) || isempty(start) || isempty(key) || isempty(beep) || isempty(en)
    return;
end

% Skip warm-up DIN pulses to align with the labeled trials
skip4 = max(0, numel(din4) - numel(key));
skip5 = max(0, numel(din5) - numel(beep));
n = min([numel(start), numel(key), numel(beep), numel(en), numel(din4)-skip4, numel(din5)-skip5]);
if n <= 0, return; end

din4u = din4(skip4 + (1:n));
din5u = din5(skip5 + (1:n));
offStart = start(1:n) - din4u;
offKey = key(1:n) - din4u;
offBeep = beep(1:n) - din5u;
offEnd = en(1:n) - din5u;
end

function y = resampleOffsets(x, nOut)
% Resample offset values to match the number of trials.
% If a single value is given, it is replicated.
% If multiple values are given, linear interpolation is used.
x = x(:)';
if isempty(x)
    y = zeros(1,nOut);
elseif numel(x) == 1
    y = repmat(x,1,nOut);
else
    xi = linspace(1,numel(x),nOut);
    y = interp1(1:numel(x),x,xi,'linear','extrap');
end
end

function t = clampTimes(t, durationS)
% Clamp event times to the recording bounds [0, durationS).
t = max(t, 0);
t = min(t, max(durationS - 1e-6, 0));
end

function counts = countCodes(codes)
% Count occurrences of each event code.
% Returns a struct with field names derived from the code strings.
counts = struct();
u = unique(codes);
for i = 1:numel(u)
    f = matlab.lang.makeValidName(u{i});
    counts.(f) = sum(strcmp(codes,u{i}));
end
end

function n = getCount(s, fieldName)
% Safely get a field value from a struct, defaulting to 0 if missing.
f = matlab.lang.makeValidName(fieldName);
if isfield(s,f), n = s.(f); else, n = 0; end
end

function t = detectFallingEdges(x, fs)
% Detect falling edges in a trigger signal.
% Threshold is set at the midpoint between the median and minimum of the
% first 20000 samples. Returns edge times in seconds.
head = sort(x(1:min(20000,numel(x))));
med = head(floor(numel(head)/2)+1);
thr = (med + min(x))/2;
mask = x < thr;
on = find(diff([false mask])==1);
t = (on-1)/fs;
end

function [tClean, info] = cleanTriggerEdges(tRaw, expectedCount, minGapS)
% Clean trigger edges by removing noise spikes via minimum-gap filtering.
%
% If an expected DIN count is provided, the gap is adaptively increased
% until the number of remaining edges is at most 1.6x the expected count
% (plus a margin of 8).
%
% Returns cleaned edge times and an info struct with cleaning statistics.
tRaw = sort(double(tRaw(:)'));
tClean = tRaw;
if nargin < 3 || ~isfinite(minGapS) || minGapS <= 0
    minGapS = 0.01;
end

if ~isempty(tRaw)
    gap = minGapS;
    tClean = applyMinGap(tRaw, gap);
    if expectedCount > 0
        maxExpected = max(ceil(1.6 * expectedCount), expectedCount + 8);
        while numel(tClean) > maxExpected && gap < 0.5
            gap = gap * 1.5;
            tClean = applyMinGap(tRaw, gap);
        end
        minGapS = gap;
    end
end

info = struct('n_raw',numel(tRaw), ...
    'n_clean',numel(tClean), ...
    'n_removed',max(0, numel(tRaw)-numel(tClean)), ...
    'min_gap_s',minGapS, ...
    'ratio_to_expected',NaN);
if expectedCount > 0
    info.ratio_to_expected = numel(tClean) / expectedCount;
end
end

function tOut = applyMinGap(tIn, minGapS)
% Remove trigger edges that are closer than minGapS to the previous edge.
% Keeps the first edge in each valid group.
tIn = sort(double(tIn(:)'));
if isempty(tIn)
    tOut = tIn;
    return;
end
keep = true(size(tIn));
last = tIn(1);
for i = 2:numel(tIn)
    if (tIn(i) - last) < minGapS
        keep(i) = false;
    else
        last = tIn(i);
    end
end
tOut = tIn(keep);
end

function checkTriggerCountPlausibility(nAux, nDin, auxName, dinName)
% Verify that the number of cleaned AUX edges is within a plausible range
% relative to the expected DIN count (60%-180%).
% Throws an error if the count is outside this range.
if nDin <= 0
    return;
end
low = max(10, floor(0.60 * nDin));
high = ceil(1.8 * nDin);
if nAux < low || nAux > high
    error('Trigger count mismatch after cleaning: %s=%d vs %s=%d (allowed [%d,%d]).', ...
        auxName, nAux, dinName, nDin, low, high);
end
end

function [a,b,mad] = fitLinear(x,y)
% Fit an affine model y = a*x + b using least squares (polyfit degree 1).
% Returns slope a, intercept b, and median absolute deviation (MAD).
x = double(x(:)); y = double(y(:));
p = polyfit(x,y,1); a = p(1); b = p(2);
mad = median(abs(y-(a*x+b)));
end

function [xUse, yUse, info] = pairDin1ForRefit(aux2, din1, a, b)
% Pair AUX2 (DIN1 sync wave) edges with EEG DIN1 events to refine the
% affine clock-drift fit.
%
% The initial fit (from DIN4/DIN5) is used to predict DIN1 times from AUX2.
% Matching is done with an adaptive tolerance based on the DIN1 step size.
% Outliers are removed using a 3*MAD threshold, and the result is
% subsampled to at most 400 pairs.
%
% Returns the paired points and an info struct with matching statistics.
% The refinement is only applied if at least 20 pairs are found with a
% match ratio >= 30%.
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

% Predict DIN1 times from AUX2 using current affine transform
predDin1 = a * aux2 + b;
din1Step = median(diff(din1));
if ~isfinite(din1Step) || din1Step <= 0
    din1Step = 0.1;
end
tolS = min(0.05, max(0.008, 0.35 * din1Step));

% Match predicted to actual DIN1 events within tolerance
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

% Remove outlier pairs using 3*MAD on residuals
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

% Subsample to prevent overweighting the dense DIN1 train vs DIN4/DIN5
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
% Pair two sorted sequences within a given tolerance.
% Uses a simple two-pointer merge: if |a - b| <= tol, they are paired;
% otherwise the smaller value advances.
% Returns indices of paired elements in each sequence.
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

function s = oneToken(txt, pat)
% Extract the first capture group matching a regex pattern from a string.
% Throws if the pattern is not found.
t = regexp(txt, pat, 'tokens','once');
if isempty(t), error('Pattern not found: %s', pat); end
s = t{1};
end

function epochSec = isoToEpoch(isoStr)
% Convert an ISO-8601 timestamp string to Unix epoch seconds.
d = datetime(isoStr,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','UTC');
epochSec = posixtime(d);
end

function isoStr = epochToIso(epochSec, offsetStr)
% Convert Unix epoch seconds to an ISO-8601 timestamp string with the
% given timezone offset. Rounds to microseconds.
offSec = offsetToSeconds(offsetStr);
localEpoch = round((epochSec + offSec)*1e6)/1e6;
d = datetime(localEpoch,'ConvertFrom','posixtime','TimeZone','UTC');
yy = year(d); mm = month(d); dd = day(d);
HH = hour(d); MN = minute(d);
ssF = second(d); ss = floor(ssF); mic = round((ssF-ss)*1e6);
if mic >= 1e6
    d = d + seconds(1);
    yy = year(d); mm = month(d); dd = day(d);
    HH = hour(d); MN = minute(d); ss = floor(second(d)); mic = 0;
end
isoStr = sprintf('%04d-%02d-%02dT%02d:%02d:%02d.%06d%s',yy,mm,dd,HH,MN,ss,mic,offsetStr);
end

function s = offsetToSeconds(off)
% Convert a timezone offset string ('+HH:MM' or '-HH:MM') to seconds.
if numel(off)~=6 || off(4)~=':'
    error('Invalid offset: %s', off);
end
hh = str2double(off(2:3)); mm = str2double(off(5:6));
sgn = 1; if off(1)=='-', sgn=-1; end
s = sgn*(hh*3600+mm*60);
end

function writeCsv(pathCsv, header, rows)
% Write a cell array of rows to a CSV file with proper escaping.
% header: cell array of column names
% rows: M x N cell array where each cell is a string to be escaped
fid = fopen(pathCsv,'w');
fprintf(fid,'%s\n', strjoin(header,','));
for i = 1:size(rows,1)
    row = rows(i,:);
    for j = 1:numel(row), row{j} = csvEscape(row{j}); end
    fprintf(fid,'%s\n', strjoin(row,','));
end
fclose(fid);
end

function s = csvEscape(v)
% Escape a value for CSV: double-quote if it contains commas, quotes,
% or newlines. Internal quotes are doubled.
if ~ischar(v), v = char(string(v)); end
if contains(v,'"'), v = strrep(v,'"','""'); end
if contains(v,',') || contains(v,'"') || contains(v,newline)
    s = ['"' v '"'];
else
    s = v;
end
end

function rowsOut = buildUnifiedDriftRows(summaryHeader, summaryRows)
% Build a unified drift report with three types of rows:
%   'session':    single drift measurement per (pair, phase)
%   'pre_post':   delta drift (post - pre) for matching subject keys
%   'right_left': delta drift (sx - dx) within the same phase
%
% This provides a compact overview of clock-drift stability across sessions,
% pre/post intervention changes, and left/right differences.
records = collectDriftRecords(summaryHeader, summaryRows);
rowsOut = {};
if isempty(records)
    return;
end

% 1) Per-session drift rows.
for i = 1:numel(records)
    rowsOut(end+1,:) = { ... %#ok<AGROW>
        'session', ...
        records(i).prepostKey, ...
        records(i).phase, ...
        records(i).pair, ...
        '', ...
        records(i).side, ...
        '', ...
        fmtNum(records(i).driftPpm, 3), ...
        '', ...
        '', ...
        '', ...
        '', ...
        '', ...
        '', ...
        'single_session_measure'};
end

% 2) Pre/Post drift deltas (post - pre).
keys = unique({records.prepostKey});
for i = 1:numel(keys)
    k = keys{i};
    iPre = find(strcmp({records.prepostKey}, k) & strcmp({records.phase}, 'pre'), 1, 'first');
    iPost = find(strcmp({records.prepostKey}, k) & strcmp({records.phase}, 'post'), 1, 'first');
    if isempty(iPre) || isempty(iPost)
        continue;
    end

    dPre = records(iPre).driftPpm;
    dPost = records(iPost).driftPpm;
    dDelta = dPost - dPre;
    if dDelta > 0
        trend = 'increase_post';
    elseif dDelta < 0
        trend = 'decrease_post';
    else
        trend = 'stable';
    end

    rowsOut(end+1,:) = { ... %#ok<AGROW>
        'pre_post', ...
        k, ...
        'pre_post', ...
        records(iPre).pair, ...
        records(iPost).pair, ...
        'pre', ...
        'post', ...
        fmtNum(dPre, 3), ...
        fmtNum(dPost, 3), ...
        fmtNum(dDelta, 3), ...
        fmtNum(abs(dDelta), 3), ...
        trend, ...
        '', ...
        '', ...
        'post_minus_pre'};
end

% 3) Right/Left drift deltas (sx - dx) within each phase.
hasSide = strcmp({records.side}, 'dx') | strcmp({records.side}, 'sx');
recordsRL = records(hasSide);
combo = cell(1, numel(recordsRL));
for i = 1:numel(recordsRL)
    combo{i} = sprintf('%s__%s', recordsRL(i).phase, recordsRL(i).sideKey);
end
u = unique(combo);

for i = 1:numel(u)
    idx = find(strcmp(combo, u{i}));
    idxDx = idx(strcmp({recordsRL(idx).side}, 'dx'));
    idxSx = idx(strcmp({recordsRL(idx).side}, 'sx'));
    if isempty(idxDx) || isempty(idxSx)
        continue;
    end

    iDx = idxDx(1);
    iSx = idxSx(1);
    dDx = recordsRL(iDx).driftPpm;
    dSx = recordsRL(iSx).driftPpm;
    dDelta = dSx - dDx;
    if abs(dDelta) < eps
        trend = 'equal';
    elseif dDelta > 0
        trend = 'sx_higher';
    else
        trend = 'dx_higher';
    end

    rowsOut(end+1,:) = { ... %#ok<AGROW>
        'right_left', ...
        recordsRL(iDx).sideKey, ...
        recordsRL(iDx).phase, ...
        recordsRL(iDx).pair, ...
        recordsRL(iSx).pair, ...
        'dx', ...
        'sx', ...
        fmtNum(dDx, 3), ...
        fmtNum(dSx, 3), ...
        fmtNum(dDelta, 3), ...
        fmtNum(abs(dDelta), 3), ...
        trend, ...
        sprintf('%d', numel(idxDx)), ...
        sprintf('%d', numel(idxSx)), ...
        'sx_minus_dx'};
end
end

function records = collectDriftRecords(summaryHeader, summaryRows)
% Parse the summary CSV rows into drift record structs for report generation.
% Only rows with status 'OK' and valid drift_ppm values are included.
% Side information is extracted from the pair name.
records = struct('pair',{},'phase',{},'prepostKey',{},'side',{},'sideKey',{}, ...
    'driftPpm',{},'madMs',{},'a',{},'b',{});
if isempty(summaryRows)
    return;
end

idxPair = find(strcmp(summaryHeader, 'pair'), 1);
idxPhase = find(strcmp(summaryHeader, 'phase'), 1);
idxStatus = find(strcmp(summaryHeader, 'status'), 1);
idxDrift = find(strcmp(summaryHeader, 'drift_ppm'), 1);
idxMad = find(strcmp(summaryHeader, 'mad_ms'), 1);
idxA = find(strcmp(summaryHeader, 'a'), 1);
idxB = find(strcmp(summaryHeader, 'b'), 1);

for i = 1:size(summaryRows, 1)
    st = upper(strCell(summaryRows{i, idxStatus}));
    if ~strcmp(st, 'OK')
        continue;
    end

    driftPpm = str2double(strCell(summaryRows{i, idxDrift}));
    if ~isfinite(driftPpm)
        continue;
    end

    pairName = strCell(summaryRows{i, idxPair});
    phase = lower(strCell(summaryRows{i, idxPhase}));
    side = detectSideFromName(pairName);
    if isempty(side)
        side = 'unknown';
    end

    records(end+1) = struct( ... %#ok<AGROW>
        'pair', pairName, ...
        'phase', phase, ...
        'prepostKey', prePostKeyFromPairLocal(pairName), ...
        'side', side, ...
        'sideKey', rightLeftKeyFromPair(pairName), ...
        'driftPpm', driftPpm, ...
        'madMs', str2double(strCell(summaryRows{i, idxMad})), ...
        'a', str2double(strCell(summaryRows{i, idxA})), ...
        'b', str2double(strCell(summaryRows{i, idxB})));
end
end

function s = prePostKeyFromPairLocal(pairName)
% Derive a key for matching pre/post pairs by stripping phase/session labels.
% E.g., "sub01_sespre_dx" -> "sub01_dx" so pre and post can be compared.
s = lower(char(pairName));
s = regexprep(s, 'ses[-_]?pre', 'ses');
s = regexprep(s, 'ses[-_]?post', 'ses');
s = regexprep(s, '(^|[_-])pre([_-]|$)', '$1$2');
s = regexprep(s, '(^|[_-])post([_-]|$)', '$1$2');
s = regexprep(s, '[_-]+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'pair'; end
end

function s = rightLeftKeyFromPair(pairName)
% Derive a key for matching right/left pairs by stripping side tags.
% E.g., "sub01_pre_dx" -> "sub01_pre" so dx/sx can be compared.
s = lower(char(pairName));
s = regexprep(s, '(^|[_-])(dx|sx|right|left|rt|lt|des|sin|destro|sinistro)([_-]|$)', '$1$3');
s = regexprep(s, '[_-]+', '_');
s = regexprep(s, '^_|_$', '');
if isempty(s), s = 'pair'; end
end

function side = detectSideFromName(pairName)
% Detect whether a pair name refers to right ('dx') or left ('sx') side.
% Uses token normalization first, then falls back to regex on the raw name.
toks = normalizeTokens(pairName);
side = extractSideToken(toks);
if isempty(side)
    s = lower(char(pairName));
    if ~isempty(regexp(s, '(^|[_-])(dx|right|rt|des|destro)([_-]|$)', 'once'))
        side = 'dx';
    elseif ~isempty(regexp(s, '(^|[_-])(sx|left|lt|sin|sinistro)([_-]|$)', 'once'))
        side = 'sx';
    else
        side = '';
    end
end
end

function s = strCell(v)
% Convert a cell entry (char, string, numeric, empty) to a plain char.
if isempty(v)
    s = '';
elseif ischar(v)
    s = v;
elseif isstring(v)
    s = char(v);
else
    s = char(string(v));
end
end

function s = fmtNum(v, nDec)
% Format a numeric value as a string with nDec decimal places.
% Returns 'NaN' for non-finite values.
if nargin < 2 || isempty(nDec)
    nDec = 3;
end
if ~isfinite(v)
    s = 'NaN';
else
    s = sprintf(['%0.' num2str(nDec) 'f'], v);
end
end

function n = stripExt(pathIn, ext)
% Strip extension from a file path. If ext is provided and matches,
% the extension is removed; otherwise the base name is returned.
% Note: this function currently returns only the base name due to the
% trailing `if` block that does nothing.
[~,n,e] = fileparts(pathIn);
if ~strcmpi(e,ext)
    % keep base name only
end
end

function cleanupTemp(tmpDir)
% Safely remove a temporary directory (with all contents).
if isfolder(tmpDir)
    try, rmdir(tmpDir,'s'); catch, end
end
end
