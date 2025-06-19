function features_data = lv_feature_extractor(cfg)
% shortened version of the main function in the toolbox
% takes struct with trial (trl_ch_time), label and time and returns the same struct but
% now with features not normal EEG

if isfield(cfg,'data')
    data = cfg.data;
    if ~isfield(cfg,'method'), cfg.method = 'mean'; cfg.window=80; end
    if strcmp(cfg.method,'erp')==1, cfg.method='mean';
        if ~isfield(cfg,'fs') % uncomment for time domain smoothing
            sampling_rate = length(nearest(data.time,0):nearest(data.time,1)) - 1; warning(['sampling rate is set to: ' num2str(sampling_rate)]);
        else
            sampling_rate = cfg.fs; warning(['sampling rate is set to: ' num2str(sampling_rate)]);
        end
    end

    features_data = data;
end


switch cfg.method
    case {'phase','power'} % hilbert at every time pt .. be careful of edges
        warning('be careful of edges');
        hilbert_dat=[];
        cfg_preprocessing                 = [];
        cfg_preprocessing.bpfilter        = 'yes';
        %range = lv_tune_params('frequency band','4 8');
        cfg_preprocessing.bpfreq          = [4 8]; %str2double(split(range)');
        data_bp= ft_preprocessing(cfg_preprocessing, data);
        % matlab's hilbert
        for h=1:size(data_bp.trial,1) %trls
            for j=1:size(data_bp.trial,2) %ch
                hilbert_dat(h,j,:) = hilbert( squeeze(data_bp.trial(h,j,:)) ); % time should be in the first dimension.
            end
        end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
        if strcmp(cfg.method,'power')==1
            features_data.trial = abs(hilbert_dat).^2;
        else
            phase_vals = angle(hilbert_dat);
            features_data.trial = [cos(phase_vals) sin(phase_vals)];
        end

    case 'spectrum'
        % the third dim should be time .. keeps the dims and get the
        % fourier coeffs of the signals
        dat = data.trial;
        nyquistVar = cfg.sampling_rate/2; N = size(dat,3);
        pts_hz = linspace(0,nyquistVar,(N/2)+1);
        coeff = fft(dat,[],3);
        features_data.trial = (abs(coeff(:,:,1:length(pts_hz)))*2)/N;  % spectrums
        features_data.pts_hz = pts_hz;

    otherwise
        window_samples = (sampling_rate*cfg.window)/1000;
        ws = floor(window_samples/2);
        time_idx = ws+1 : size(data.trial,3)-ws; % from the first possible window in time
        for tt=1:length(time_idx)
            % any function that acts on a dimension as SECOND argument can be used: mean, median.. etc,
            features_data.trial(:,:,time_idx(tt)) = eval([cfg.method '( data.trial(:,:,time_idx(tt)-ws : time_idx(tt)+ws) ,3)']);
        end
        features_data.trial = features_data.trial(:,:,time_idx);
        features_data.time = features_data.time(1,time_idx);

end



end


% helping functions
function [x,y,z] = hilbert3(n)
% Hilbert 3D curve.
%
% function [x,y,z] = hilbert3(n) gives the vector coordinates of points
% in n-th order Hilbert curve of area 1.
%
% Example: plot the 3-rd order curve
%
% [x,y,z] = hilbert3(3); plot3(x,y,z)
%   Copyright (c) by Ivan Martynov
%   Inspired by function [x,y] = hilbert(n) made by Federico Forte
%   Date: September 17, 2009
if nargin ~= 1
    n = 2;
end

if n <= 0
    x = 0;
    y = 0;
    z = 0;
else
    [xo,yo,zo] = hilbert3(n-1);
    x = .5*[.5+zo .5+yo -.5+yo -.5-xo -.5-xo -.5-yo .5-yo .5+zo];
    y = .5*[.5+xo .5+zo .5+zo .5+yo -.5+yo -.5-zo -.5-zo -.5-xo];
    z = .5*[.5+yo -.5+xo -.5+xo .5-zo .5-zo -.5+xo -.5+xo .5-yo];
end
end


function cov_data_pooled=pooled_cov(data)
% gets the pooled covariance of different trials (trials in cells) or 3d
if length(size(data))==3, for i=1:size(data,1), data_temp{1,i} = squeeze(data(i,:,:)); end, data=data_temp; end % if 3d make it 2d
data_centered = cellfun(@(x) (x-repmat(mean(x,2),1,size(x,2)) ), data,'Un',0);
cov_data = cellfun(@(x) ((x*x')./size(x,2)), data_centered,'Un',0);
cov_data_pooled=zeros(size(cov_data{1,1}));
for i=1:length(data), cov_data_pooled = cov_data_pooled + cov_data{1,i}; end

cov_data_pooled = cov_data_pooled./length(data);

end



function topodata = topomaps(lv_layout,dat)
% from Mike's https://www.youtube.com/watch?v=VyAAW_fWo6M&t=1232s&ab_channel=MikeXCohen
% takes one vector and does interpoaltion with respect to the channels' locations
% topographical maps
% mikeXcohen@gmail.com

% get cartesian coordinates
elocsX = lv_layout.pos(:,1); elocsY = lv_layout.pos(:,2);

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX),max(elocsX),interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);

% now interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsX,elocsY,double(dat'), 'natural');
topodata = interpFunction(gridX,gridY);
topodata = single(topodata);
% figure, imagesc(topodata)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION from fieldtrip's crossfrequency for PAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacdata,lv_phase_angle] = data2mvl(LFsigtemp,HFsigtemp) % Canolty 2006
% calculate  mean vector length (complex value) per trial
% mvldata dim: LF*HF
% each row is a different frequency for phase inside LF for which you want
% PAC and also each row is a different frequency for amplitude inside HF
LFphas   = LFsigtemp;
HFamp    = HFsigtemp;
mvldata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));    % mean vector length

for i = 1:size(LFsigtemp,1)
    for j = 1:size(HFsigtemp,1)
        mvldata(i,j) = nanmean(HFamp(j,:).*exp(1i*LFphas(i,:)));
    end
end

pacdata = abs(mvldata); % the magnitude of the modulation index is the coupling strength
lv_phase_angle = angle(mvldata); % angle is phase
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacdata,lv_phase_angle] = data2pac(LFsigtemp,HFsigtemp,nbin)
% calculate phase amplitude distribution per trial
% pacdata dim: LF*HF*Phasebin

pacdata = zeros(size(LFsigtemp,1),size(HFsigtemp,1),nbin);

Ang  = LFsigtemp;
Amp  = HFsigtemp;
[dum,bin] = histc(Ang, linspace(-pi,pi,nbin));  % binned low frequency phase
binamp = zeros (size(HFsigtemp,1),nbin);      % binned amplitude

for i = 1:size(Ang,1)
    for k = 1:nbin
        idx = (bin(i,:)==k);
        binamp(:,k) = mean(Amp(:,idx),2);
    end
    pacdata(i,:,:) = binamp;
end

% lv: getting the phase of the amplitudes
bins = linspace(-pi,pi,nbin)';
% we convert the amplitudes into weights for the binned angles so that the circular mean of their amplitude weighted values is the phase of coupling
for i=1:size(pacdata,1) %freqlow
    for j=1:size(pacdata,2) %freqhigh
        bins_temp = bins;
        w = squeeze(pacdata(i,j,:))./nansum(squeeze(pacdata(i,j,:)));
        bins_temp(isnan(w))=[]; w(isnan(w))=[];
        lv_phase_angle(i,j) = circ_mean(bins_temp, w);
    end
end

% now the KL distances from the uniform distribution for every point
nlf = size(Ang, 1); nhf = size(Amp, 1);
% from fieldtrip, the part of getting the KL distance
Q =ones(nbin,1)/nbin; % uniform distribution
mi = zeros(nlf,nhf);
for i=1:nlf
    for j=1:nhf
        P = squeeze(pacdata(i,j,:))/ nansum(pacdata(i,j,:));  % normalized distribution
        mi(i,j) = nansum(P.* log2(P./Q))./log2(nbin); % KL distance
    end
end
pacdata = mi;

end % function