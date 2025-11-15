
[x1] = audioread("Short_QuranPalestine.wav");
[x2] = audioread("Short_BBCArabic2.wav");
[x3] = audioread("Short_FM9090.wav");
[x4] = audioread("Short_RussianVoice.wav");
[x5, fs] = audioread("Short_SkyNewsArabia.wav");

audioarray = {x1, x2, x3, x4, x5};

% creating 1 single channel by averaging both channels
for i = 1:5
    audioarray{i} = (audioarray{i}(:,1) + audioarray{i}(:,2)) / 2;
end

% padding all function with zero for them to be the same length
maxN = max(cellfun(@length, audioarray));
for i = 1:5
    audioarray{i} = [audioarray{i}; zeros(maxN - length(audioarray{i}), 1)];
end

f = (-maxN/2:maxN/2-1) * (fs/maxN); % Frequency Axis

figure; hold on;
for k = 1:5
    X = fft(audioarray{k});
    plot(f, fftshift(abs(X)));
end
xlabel('Frequency (Hz)');
ylabel('FFT');
title('FFT of Message Signals');
hold off;

 % We can observe BW to be equal to around 10KHZ

  for i = 1:5
    audioarray{i} = interp(audioarray{i}, 15);
 end
 fs_new = 15*fs;
 maxN_new = length(audioarray{1});

omega = 2*pi*100e3;  % carrier radian freq
delta_f = 55e3;
Ts = 1/fs_new;
t = (0:maxN_new-1)' * Ts;

carrier_arr = cell(1,5);

for n = 0:4
    carrier_arr{n+1} = cos(omega*t + 2*pi*delta_f*n*t);
end

f_new = (-maxN_new/2:maxN_new/2 - 1) * (fs_new/maxN_new);

figure; hold on;
for k = 1:5
    A = fft(carrier_arr{k});
    plot(f_new, fftshift(abs(A)));
end
xlabel('Frequency (Hz)');
ylabel('FFT');
legend('carrier 1','carrier 2','carrier 3','carrier 4','carrier 5');
title('FFT of Carrier Signals');
hold off;

modulated_signals = cell(1, 5);

for i = 1:5
    modulated_signals{i} = audioarray{i} .* carrier_arr{i};
end

figure; hold on;
for k = 1:5
    B = fft(modulated_signals{k});
    plot(f_new, fftshift(abs(B)));
end
xlabel('Frequency (Hz)');
ylabel('FFT');
legend('signal 1','signal 2','signal 3','signal 4','signal 5');
title('FFT of Modulated Signals');
hold off;

Mixed_signal = 0;
for i = 1:5
    Mixed_signal = Mixed_signal + modulated_signals{i};
end
BPF1_out = filter(HDBPF1, Mixed_signal);
BPF1_out_FFT = fft(BPF1_out);
figure; hold on;
plot(f_new, fftshift(abs(BPF1_out_FFT)));
xlabel('Frequency (Hz)');
ylabel('RF Filter Output');
hold off;

omega_IF = 2*pi*27.5e3;
IF_freq_arr = cell(1, 5);
for n = 0:4
    IF_freq_arr{n+1} = cos(omega*t + 2*pi*delta_f*n*t + omega_IF*t);
end

IF_BPF1_out = filter(HDIFBPF1, BPF1_out .* IF_freq_arr{1});
figure; hold on;
plot(f_new, fftshift(abs(fft(IF_BPF1_out))));
xlabel('Frequency (Hz)');
ylabel('IF Stage Output');
hold off;
Demodulated_Signal1 = filter(LPF1, IF_BPF1_out .* cos(omega_IF*t));
figure;
plot(f_new, fftshift(abs(fft(Demodulated_Signal1))));
xlabel('Freqeuncy (Hz)');
ylabel('Bandpass detection output');
Demodulated_Signal1 = downsample(Demodulated_Signal1, 15);
soundsc(Demodulated_Signal1, fs);




























