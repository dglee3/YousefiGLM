function [bigStruct, D, W, samps] = generateCovariates(datmat, includeWhisker)
% David Lee
% Last update: 6/17/19

% Change log:
% 6/17/19 - Removed extra code for adding time-shifted feature vectors; old
% code is now generateCovariates_defshift.m
% 6/17/19 - Added Hilbert transform to extract phase, amplitude, and set
% point from whisker data
% 6/14/19 - Added code to determine if whisker data has too many NaNs or
% too many consecutive NaNs and ignore the data
% 6/11/19 - Added eventKey, generateEventVector, and includeWhisker option
% 6/10/19 - Implemented lick vector and whisker data import
% tic

% Name of file to import & parse
% datmat = 'jn029-10.mat';
D = load(datmat);

if isfield(D,'licks')
    licks = D.licks;
    lick_exist = 1;
else
    lick_exist = 0;
end

trials = D.trials;
summary = D.summary;
CaA0 = D.CaA0;

if isfield(D,'CaA2')
    CaA2 = D.CaA2;
elseif isfield(D,'CaA1')
    CaA2 = D.CaA1;
else
    disp('Second area (CaA1 or CaA2) not found')
    return
end

if contains(datmat,'.mat')
    datmat = datmat(1:end-4);
end

wdatmat = [datmat '_whisker.mat'];

if nargin == 1
    includeWhisker = 0;
end

if exist(wdatmat,'file') == 2 && includeWhisker
    whiskerdat_available = 1;
    W = load(wdatmat);
    %     no_whiskers = W.no_whiskers;
    whisker_dat = W.whisker_dat;
    whiskSampRate = 500; % Sampling rate of whisker data, Hz
    lpcutoff = 15; % Hz, cutoff for lowpass whisker filter
    w_lp = lpcutoff*2*pi/whiskSampRate;
    [b_lp,a_lp] = butter(2,w_lp,'low');
else
    whiskerdat_available = 0;
    W = [];
end

stimlength = 1.2 * 1000; % Length (ms) of stimulus
samps = 400; % Number of samples to include in each neuron trial
Ca_offset = 0.667 * 1000; % milliseconds from behavior start to Ca imaging start

cal = CaA0;
cal2 = CaA2;
samprate = cal.sampling_rate;

matlist = summary.table(:,contains(summary.table(1,:),'A0_Name'));

% -----------------------------------------------------------

% Do not edit
% Used for adjusting which window of data to select
minlen = samps;
samprange = 1:samps;

bigStruct = struct('trialID', [], 'Ca1_spikes',[],'Ca2_spikes',[],'timevec',[],...
    'stim1',[],'stim2',[],'CW',[],'CCW',[],'match',[],'nonmatch',[],...
    'response_HFCM', []);
%     'def_shift',def_shift);

% Events and numbers referring to each event (used in timevec)
eventKey = {'begin','present_time','direction_1_time',...
    'delay_time','delay_withdraw_time','delay_present_time','direction_2_time',...
    'report_time','reward_time','withdraw_time'};
eKnum = num2cell(1:length(eventKey));
eKdesc = {'Trial Start','Present','Direction 1 Start','Delay',...
    'Delay Withdraw','Delay Present','Direction 2 Start','Report',...
    'Reward','Withdraw'};
bigStruct.eventKey = {eventKey{:}; eKnum{:}; eKdesc{:}};

if lick_exist
    bigStruct.licks = [];
end

if whiskerdat_available
    bigStruct.whisker_touch = [];
    bigStruct.whisker_angle = [];
    bigStruct.whisk_reconstruction = [];
    bigStruct.whisk_phase = [];
    bigStruct.whisk_amp = [];
    bigStruct.whisk_setpt = [];
end

bS = 1; % bigStruct index

% Make second calcium trial info more useful (as a single struct)
cal2_ti = struct('motion_metric',[],'time_stamp',[],'fileloc',[],'mat_file',[]);
for struct_ind = 1:length(cal2.trial_info)
    cal2_ti(struct_ind) = cal2.trial_info{struct_ind};
end

% Make a list of mat files in second calcium dataset
cal2_matfiles = {cal2_ti.mat_file}';

for c1ti = 1:length(cal.trial_info)
    trialName = cal.trial_info{c1ti}.mat_file(4:end-4);
    
    % Proceed if the trial was recorded on both channels
    if any(contains(cal2_matfiles,trialName))
        
        % Add neuron data
        Ca1_spikes = cal.deconv{c1ti};
        %         Ca1_F = cal.F_dF{c1ti};
        %         ca2_ind = find(contains(cal2_matfiles,trialName));
        %         Ca2_F = cal2.F_dF{contains(cal2_matfiles,trialName)};
        Ca2_spikes = cal2.deconv{contains(cal2_matfiles,trialName)};
        
        [~, c1col] = size(Ca1_spikes);
        [~, c2col] = size(Ca2_spikes);
        
        % Some calcium data is missing a sample. Alert if more than one is
        % missing; may indicate something is misaligned
        
        % Also alert if trial is shorter than minimum trial length (samps)
        
        if abs(c1col-c2col) <= 1 && c1col >= minlen && c2col >= minlen
            
            % Begin by finding stimulus onset times (used to build feature vectors)
            sRow = find(contains(matlist,trialName),1,'last'); % Row in summary table containing trial info
            trialnum = summary.table{sRow,1};
            bigStruct(bS).trialID = trialnum;
            
            % Find times of each stimulus
            % 1 - Stim1 begin; 2 - Stim1 end; 3 - Stim2 begin; 4 - Stim2 end
            % 5 - Decision time
            sTimes = zeros(1,5);
            sTimes(1) = trials(trialnum).direction_1_time-Ca_offset;
            sTimes(2) = sTimes(1) + stimlength;
            sTimes(3) = trials(trialnum).direction_2_time-Ca_offset;
            sTimes(4) = sTimes(3) + stimlength;
            sTimes(5) = trials(trialnum).decision_time;
            
            % Cut off extra samples so all vectors have same length
            Ca1 = Ca1_spikes(:,samprange);
            Ca2 = Ca2_spikes(:,samprange);
            
            bigStruct(bS).Ca1_spikes = Ca1;
            bigStruct(bS).Ca2_spikes = Ca2;
            
            % Create time vector
            timevec = (0:samps-1)/samprate * 1000; % use ms for time unit
            eventvec = generateEventVector(timevec,trials(trialnum),bigStruct(1).eventKey,Ca_offset);
            bigStruct(bS).timevec = [timevec;eventvec];
            
            % Import lick vector
            if lick_exist
                lick_ind = find([licks(:).trial] == trialnum); % Index in licks struct for trial of interest
                lickvec = licks(lick_ind).lick_vector(1,:); % Binary lick vector
                lto = licks(lick_ind).lick_vector(2,:) * 1000 - Ca_offset; % Timestamps to ms & time-shifted
                lickTimeSamps = [find(lto<0,1,'last') find(lto>=0 & lto<=max(timevec))]; % Find samples within time range and one before t=0 for binning
                lt = lto(lickTimeSamps); % Lick time vector
                lickvec = lickvec(lickTimeSamps);
                lickvec_ds = zeros(1,samps); % Downsampled lick vector
                
                % If any licks occurred in time bin, consider as one lick
                lickvec_ds(1) = any(lickvec(lt<=0));
                for dsi = 2:samps
                    lickvec_ds(dsi) = any(lickvec(lt>timevec(dsi-1) & lt<=(timevec(dsi))));
                end
                
                bigStruct(bS).licks = lickvec_ds;
            end
            
            % Import whisker data if available
            if whiskerdat_available
                tRow_whisk = 1 + trialnum; % whisker_dat has blank info in first field; assuming all trials present in summary table are present in whisker file
                
                % Grab relevant data
                touchVec = whisker_dat(tRow_whisk).touch_vector;
                meanAng = whisker_dat(tRow_whisk).mean_angle;
                
                % Determine relevant time window
                nWsamps = min([length(touchVec) length(meanAng)]); % Number of samples in whisker data (smallest of the two)
                whiskTimevec_O = (0:nWsamps-1)./whiskSampRate * 1000 - Ca_offset; % Time in ms, beginning at behavior start, relative to 2P start
                whiskTimeSamps = [find(whiskTimevec_O<0,1,'last') find(whiskTimevec_O>=0 & whiskTimevec_O<=max(timevec))]; % Indices for samples within relevant range (decided by timevec)
                whiskTimevec = whiskTimevec_O(whiskTimeSamps);
                
                % Downsample relevant touch data
                touchVec = touchVec(whiskTimeSamps);
                touch_ds = [touchVec(1) zeros(1,samps-1)]; % Pre-allocate downsampled touch vector
                for dsi = 2:samps
                    touch_ds(dsi) = any(touchVec(whiskTimevec>=timevec(dsi-1) & whiskTimevec<timevec(dsi))); % If whisker touched during time bin, record 1
                end
                bigStruct(bS).whisker_touch = touch_ds;
                
                % Downsample relevant whisking data
                meanAng = meanAng(whiskTimeSamps); % Get relevant whisking data
                
                % meanAng may have NaN, which is a problem for filtering & GLM
                tf = isnan(meanAng);
                maxConsNaNs = maxConsecutive(tf);
                ConsNaN_Limit = 15; % Limit of consecutive NaNs before ignoring whisker data
                NaN_Limit = 0.1*numel(meanAng); % Limit of total NaNs before ignoring whisker data
                
                if sum(tf) < NaN_Limit && maxConsNaNs < ConsNaN_Limit
                    ix = 1:numel(meanAng);
                    meanAng(tf) = interp1(ix(~tf),meanAng(~tf),ix(tf)); % Interpolate NaNs
                    
                    % Low-pass filter the meanAngle signal for Hilbert
                    meanAng = filtfilt(b_lp,a_lp,meanAng);
                    
                    % Find phase, set point, and amplitude using Hilbert transform
                    bp = [4 20]; % bandpass filter for finding phase
                    setpt_func = inline( '(max(x) + min (x)) / 2');
                    amp_func  = @range ;
                    
                    phase = phase_from_hilbert( meanAng, whiskSampRate, bp );
                    [amp,tops] = get_slow_var(meanAng, phase, amp_func );
                    setpt = get_slow_var(meanAng, phase, setpt_func );
                    reconstruction = setpt + (amp/2).*cos(phase);
                    
                    % Prea-allocate vectors for filtered whisker signal,
                    % reconstruction, phase, setpoint, and amplitude vectors
                    ang_ds = [meanAng(1) zeros(1,samps-1)];
                    rec_ds = [reconstruction(1) zeros(1,samps-1)];
                    phase_ds = [phase(1) zeros(1,samps-1)];
                    amp_ds = [amp(1) zeros(1,samps-1)];
                    setpt_ds = [setpt(1) zeros(1,samps-1)];
                    
                    for dsi = 2:samps
                        ds_samp = find(whiskTimevec<=timevec(dsi),1,'last'); % Flexible downsampling: grab data point from nearest time before each point in timevec
                        ang_ds(dsi) = meanAng(ds_samp);
                        rec_ds(dsi) = reconstruction(ds_samp);
                        phase_ds(dsi) = phase(ds_samp);
                        amp_ds(dsi) = amp(ds_samp);
                        setpt_ds(dsi) = setpt(ds_samp);
                    end
                    bigStruct(bS).whisker_angle = ang_ds;
                    bigStruct(bS).whisk_reconstruction = rec_ds;
                    bigStruct(bS).whisk_phase = phase_ds;
                    bigStruct(bS).whisk_amp = amp_ds;
                    bigStruct(bS).whisk_setpt = setpt_ds;
                else
                    bigStruct(bS).whisker_angle = [];
                end
            end
            
            % Find indices of stim start & end
            sInds = zeros(1,length(sTimes));
            for sI = 1:length(sTimes)
                tv = timevec-sTimes(sI);
                sInds(sI) = find(tv<=0,1,'last');
            end
            
            % Time indices where rotation occurs
            pulse1 = sInds(1):sInds(2);
            pulse2 = sInds(3):sInds(4);
            % Pulses are at true location relative to trial start
            
            % Label stimuli directions
            delim1 = strfind(trials(trialnum).direction_1,' ');
            stim1 = trials(trialnum).direction_1(1:delim1-1);
            delim2 = strfind(trials(trialnum).direction_2,' ');
            stim2 = trials(trialnum).direction_2(1:delim2-1);
            bigStruct(bS).stim1 = stim1;
            bigStruct(bS).stim2 = stim2;
            
            % Create feature vectors for CW & CCW
            bigStruct(bS).CW = sparse(zeros(1,samps));
            bigStruct(bS).CCW = sparse(zeros(1,samps));
            
            if strcmp(stim1,'CW')
                % Label time indices for rotation
                bigStruct(bS).CW(1,pulse1) = 1;
            elseif strcmp(stim1,'CCW')
                bigStruct(bS).CCW(1,pulse1) = 1;
            else
                disp('Error')
            end
            
            if strcmp(stim2,'CW')
                bigStruct(bS).CW(1,pulse2) = 1;
            elseif strcmp(stim2,'CCW')
                bigStruct(bS).CCW(1,pulse2) = 1;
            else
                disp('Error')
            end
            
            % Create feature vectors for match & nonmatch
            bigStruct(bS).match = zeros(1,samps);
            bigStruct(bS).nonmatch = zeros(1,samps);
            
            if strcmp(stim1,stim2)
                bigStruct(bS).match(1,pulse2) = 1;
            else
                bigStruct(bS).nonmatch(1,pulse2) = 1;
            end
            
            % Code result and decision time in bigStruct.response_HFCM
            trialResult = trials(trialnum).decision;
            bigStruct(bS).response_HFCM = sparse(4,samps);
            
            if strcmp(trialResult,'Hit')
                bigStruct(bS).response_HFCM(1,sInds(5)) = 1;
            elseif strcmp(trialResult,'FA')
                bigStruct(bS).response_HFCM(2,sInds(5)) = 1;
            elseif strcmp(trialResult,'CR')
                bigStruct(bS).response_HFCM(3,sInds(5)) = 1;
            elseif strcmp(trialResult,'Miss')
                bigStruct(bS).response_HFCM(4,sInds(5)) = 1;
            end
            
            % Increment structure index
            bS = bS + 1;
            
            % Make a note and skip the data if more than one sample is missing
        elseif c1col < samps || c2col < samps
            disp(['Ignoring Trial ' num2str(c1ti) '. Trial is too short.'])
        else
            disp(['Large discrepancy in F_dF size, ignoring trial: ' trialName])
            disp(['Ca1: ' num2str(c1col) '; Ca2: ' num2str(c2col)])
            disp(' ')
        end
        
        % Make a note if a trial is missing in the second calcium dataset
    else
        str = ['Note: ' trialName ' does not exist in cal2'];
        disp(str)
        disp(' ')
    end
end
end
% toc
%%
function [eventTiming] = generateEventVector(timeVector,trialInfo,eventList,offset)
% Generates a vector containing zeros and numbers which represent events
% aligned with timeVector
% timeVector: Timing (ms) of each sample from experiment start
% trialInfo: Struct for trial of interest with time stamps for events
% eventList: Cell containing fields to extract from trialInfo in row 1 and
%              number representing each event in row 2
% offset: Timing difference between imaging start and trial start

eventTiming = zeros(1,length(timeVector));

for eventNo = 1:length(eventList)
    event = eventList{1,eventNo};
    
    if isfield(trialInfo,event)
        timestamp = trialInfo.(event) - offset;
        eInd = find(timeVector<timestamp,1,'last');
        eventTiming(eInd) = eventNo;
    elseif strcmp(event,'begin')
        eventTiming(1) = eventNo;
    end
end
end
%%
function [maxCons] = maxConsecutive(x)
% Finds max consecutive elements in a vector that are >0
y = find(~(x>0)); % Indices of elements >0
diffs = diff([0 y numel(x) + 1])-1; % Distance between zero elements
maxCons = max(diffs);
end