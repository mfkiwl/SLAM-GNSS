function acqResults = BDS_acquisition(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Condition input signal to speed up acquisition =========================

% If input IF signal freq. is too high, a resampling strategy is applied
% to speed up the acquisition. This is user selectable.
if (settings.samplingFreq > settings.resamplingThreshold && ...
                                       settings.resamplingflag == 1)
   
    %--- Filiter out signal power outside the main lobe of CM code --------
    fs = settings.samplingFreq;
    IF = settings.IF;
    % Bandwidth of CA mian lobe: 0.5e6 is a margin to make sure most of CA
    % code power will not filtered out
    BW = settings.codeFreqBasis*2 + 0.5e6;
    % Filter parameter
    w1 = (IF)-BW/2;
    w2 = (IF)+BW/2;
    wp = [w1*2/fs w2*2/fs];
    % Filter coefficients
    b  = fir1(700,wp);
    % Filter operation
    longSignal = filtfilt(b,1,longSignal);
    
    % --- Find resample frequency -----------------------------------------
    % Refer to bandpass sampling theorem(Yi-Ran Sun,Generalized Bandpass
    % Sampling Receivers for Software Defined Radio)
    
    % Upper boundary frequency of the bandpass IF signal
    fu = settings.IF + BW/2;
    
    % Lower freq. of the acceptable sampling Freq. range
    n = floor(fu/BW);
    if (n<1)
        n = 1;
    end
    lowerFreq = 2*fu/n;
    
    % Lower boundary frequency of the bandpass IF signal
    fl = settings.IF - BW/2;
    
    % Upper boundary frequency of the acceptable sampling Freq. range
    if(n>1)
        upperFreq = 2*fl/(n-1);
    else
        upperFreq = lowerFreq;
    end
    
    % Save orignal Freq. for later use
    oldFreq = settings.samplingFreq;
    
    % Take the center of the acceptable sampling Freq. range as
    % resampling frequency. As settings are used to generate local
    % CM code samples, so assign the resampling freq. to settings.
    % This can not change the settings.samplingFreq outside this 
    % acquisition function.
    settings.samplingFreq = ceil((lowerFreq + upperFreq)/2);
    
    %--- Downsample input IF signal -------------------------------------------------
    % Signal length after resampling
    signalLen = floor((length(longSignal)-1) /oldFreq * settings.samplingFreq);
    % Resampled signal index
    index = ceil((0:signalLen-1)/settings.samplingFreq * oldFreq);
    index(1) = 1;
    % Resampled signal
    longSignal = longSignal(index);
    
    % For later use
    oldIF = settings.IF;
    
    % Resampling is equivalent to down-converting the original IF by integer
    % times of resampling freq.. So the IF after resampling is equivalent to:
    settings.IF = rem(settings.IF,settings.samplingFreq); 
end % resampling input IF signals

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

AcqCodeLength = 4;                    
% Create two 1msec vectors of data to correlate with and one with zero DC
signal1 = longSignal(1 : AcqCodeLength*samplesPerCode);
signal2 = longSignal(AcqCodeLength*samplesPerCode+1   : 2*AcqCodeLength*samplesPerCode);
signal3 = longSignal(AcqCodeLength*samplesPerCode*2+1 : 3*AcqCodeLength*samplesPerCode);
signal4 = longSignal(AcqCodeLength*samplesPerCode*3+1 : 4*AcqCodeLength*samplesPerCode);
signal5 = longSignal(AcqCodeLength*samplesPerCode*4+1 : 5*AcqCodeLength*samplesPerCode);

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (AcqCodeLength*samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
% the frequency step is 500/AcqCodeLength
numberOfFrqBins = round(settings.acqSearchBand * AcqCodeLength *2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings);

save BDS_caCodesTable.txt -ascii caCodesTable

%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode*AcqCodeLength);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 37);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 37);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 37);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
	caCodeFreqDom = conj(fft([caCodesTable(PRN, :), caCodesTable(PRN, :) zeros(1, 2*samplesPerCode)]));
                           
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins %*AcqCodeLength

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 / AcqCodeLength * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(1i*frqBins(frqBinIndex) * phasePoints);
        
        %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* signal1);
        Q1      = imag(sigCarr .* signal1);
        I2      = real(sigCarr .* signal2);
        Q2      = imag(sigCarr .* signal2);
        I3      = real(sigCarr .* signal3);
        Q3      = imag(sigCarr .* signal3);
        I4      = real(sigCarr .* signal4);
        Q4      = imag(sigCarr .* signal4);
        I5      = real(sigCarr .* signal5);
        Q5      = imag(sigCarr .* signal5);
        
        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1i*Q1);
        IQfreqDom2 = fft(I2 + 1i*Q2);
        IQfreqDom3 = fft(I3 + 1i*Q3);
        IQfreqDom4 = fft(I4 + 1i*Q4);
        IQfreqDom5 = fft(I5 + 1i*Q5);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
        convCodeIQ3 = IQfreqDom3 .* caCodeFreqDom;
        convCodeIQ4 = IQfreqDom4 .* caCodeFreqDom;
        convCodeIQ5 = IQfreqDom5 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        acqRes1 = abs(ifft(convCodeIQ1));
        acqRes2 = abs(ifft(convCodeIQ2));
        acqRes3 = abs(ifft(convCodeIQ3));
        acqRes4 = abs(ifft(convCodeIQ4));
        acqRes5 = abs(ifft(convCodeIQ5));
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        results(frqBinIndex, :) = acqRes1;
        
        [~, MaxIndex] = max([max(acqRes1) max(acqRes2) max(acqRes3) ... 
            max(acqRes4) max(acqRes5)]);

        switch MaxIndex
            case 1
                results(frqBinIndex, :) = acqRes1;
            case 2
                results(frqBinIndex, :) = acqRes2;
            case 3
                results(frqBinIndex, :) = acqRes3;
            case 4
                results(frqBinIndex, :) = acqRes4;
            case 5
                results(frqBinIndex, :) = acqRes5;
        end                
                
    end % frqBinIndex = 1:numberOfFrqBins

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [~, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(results));

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        acqResults.codePhase(PRN) = codePhase;
        acqResults.carrFreq(PRN)  = settings.IF - ...
                                    (settings.acqSearchBand/2) * 1000 + ...
                                    0.5e3 / AcqCodeLength * (frequencyBinIndex - 1);
        
%% Downsampling recovery ============================================
        % Find acquisition results corresponding to orignal sampling freq
        if (exist('oldFreq', 'var') && settings.resamplingflag == 1)
            % Find code phase
            acqResults.codePhase(PRN) = floor((codePhase - 1)/ ...
                settings.samplingFreq * oldFreq)+1;

            % Doppler frequency
            doppler = acqResults.carrFreq(PRN) - settings.IF;
                                    
            % Carrier freq. corresponding to orignal sampling freq
            acqResults.carrFreq(PRN) = doppler + oldIF;
        end
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
