function yixin(signal)
% Fs is recording frequency
  Fs=25000; 
% Butterworth Bandstop filter designed using FDESIGN.BANDSTOP.
% All frequency values are in Hz.
%     Fpass1 = 50;          % First Passband Frequency
%     Fstop1 = 55;          % First Stopband Frequency
%     Fstop2 = 65;          % Second Stopband Frequency
%     Fpass2 = 70;          % Second Passband Frequency
%     Apass1 = 0.05;        % First Passband Ripple (dB)
%     Astop  = 60;          % Stopband Attenuation (dB)
%     Apass2 = 0.05;        % Second Passband Ripple (dB)
%     match  = 'stopband';  % Band to match exactly
% 
% % Construct an FDESIGN object and call its BUTTER method.
%     h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
%                           Apass2, Fs);
%     Hd_bandstop = design(h, 'butter', 'MatchExactly', match);
%     signal_bandstop=filter(Hd_bandstop,signal);
%   
%     figure;
%     plot((0:length(signal)-1)/Fs,signal)
%     hold on
%     plot((0:length(signal_bandstop)-1)/Fs,signal_bandstop,'r')
%     hold off
% % [EOF]

Fstop = 1;           % Stopband Frequency
Fpass = 3;           % Passband Frequency
Astop = 80;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd_highpass = design(h, 'butter', 'MatchExactly', match);

  signal_highpass=filter(Hd_highpass,signal);
  
  figure;
  plot((0:length(signal)-1)/Fs,signal,'b')
  hold on
  plot((0:length(signal_highpass)-1)/Fs,signal_highpass,'r')
  hold off
  
%----------------Lowpass digital Filter----------------%
  Fpass = 2000;            % Passband Frequency
  Fstop = 2250;            % Stopband Frequency
  Dpass = 0.057501127785;  % Passband Ripple
  Dstop = 0.0001;          % Stopband Attenuation
  dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
  [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
  coeff_lowpass  = firpm(N, Fo, Ao, W, {dens});
  Hd = dfilt.dffir(coeff_lowpass);
  signal_lowpass = filter(Hd, signal_highpass);

  
% The following plot is to plot the raw signal and the filtered signal
  figure
    plot((0:length(signal_highpass)-1)/Fs,signal,'b');
    hold on
  plot((0:length(signal_highpass)-1)/Fs,signal_highpass,'r');
  xlabel('Time (s)');
  ylabel('Amplitude (uV)');
  hold on          % to compare raw data and filtered data in one figure
  plot((0:length(signal_lowpass)-1)/Fs,signal_lowpass,'g');
  xlabel('Time (s)');
  ylabel('Amplitude (uV)');
  hold off

  figure
  pwelch(signal,[],[],[],Fs,'onesided');
  hold on
  pwelch(signal_lowpass,[],[],[],Fs,'onesided');
  hold off
  
% ---------------Signal to Noise ratio calculation----------------%
% The following function can be used for spike sorting
% To find the amplitude below the signal baseline
%   signal_f = signal_f_temp(length(signal_f_temp)-Fs*10+1:length(signal_f_temp));
  signal_f = signal_lowpass(10*Fs:20*Fs-1);
  poss_low = signal_f<-150;                                   
  left_low = find(diff([poss_low])==1); 
  right_low = find(diff([poss_low])== -1); 
 min_low_spikes (length(left_low), 1)=0;
    for min_index = 1:length(left_low)
         [minvalue(min_index) minlocation(min_index)] = min(signal_f(left_low(min_index):right_low(min_index)));
          minlocation(min_index) = minlocation(min_index)+left_low(min_index)-2; 
          min_value_location(1,min_index) = minlocation(min_index);
          min_value(1,min_index) = minvalue(min_index);
          min_low_spikes (min_index, 1) = min(signal_f(left_low(min_index):right_low(min_index)));
    end
    figure
    plot((0:length(signal_f)-1)/Fs,signal_f)   %generate X value
    hold on
    plot(min_value_location/Fs, min_value, '*r'); % locate the red stars
    xlabel('Time (s)');
    ylabel('Amplitude (uV)');
   

% To find the amplitude upper the signal baseline
  poss_high = signal_f>50;                                   
  left_high = find(diff([poss_high])==1); 
  right_high = find(diff([poss_high])== -1); 
  max_high_spikes (length(right_high), 1) = 0;
    for max_index = 1:length(left_high)
         [maxvalue(max_index) maxlocation(max_index)] = max(signal_f(left_high(max_index):right_high(max_index)));
          maxlocation(max_index) = maxlocation(max_index)+left_high(max_index)-2; 
          max_value_location(1,max_index) = maxlocation(max_index);
          max_value(1,max_index) = maxvalue(max_index);
          max_high_spikes (max_index, 1) = max(signal_f(left_high(max_index):right_high(max_index)));
    end
  plot(max_value_location/Fs, max_value, '*g');
  xlabel('Time (s)');
  ylabel('Amplitude (uV)');
  hold off
  
  
  
end

