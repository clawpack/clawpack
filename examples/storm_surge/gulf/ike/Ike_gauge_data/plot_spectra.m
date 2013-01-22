% Simple program to help with plotting spectra
clear
close all

% Change the station number to the one you want

gauge = 'X'

figure(1)
clf

% loads yd_processed, Hs_wave, mean_water, peak_freq
str1 = ['load result_',gauge]
eval(str1)

% Steps through time - hit any key to move to next time step
for m = 1:length(Hs_wave)

    yd_processed(m)
    Hs = sqrt(sum(S_q(m,:))*(F(2)-F(1)))*4
    
    clf
plot(F, S_q(m,:), 'k-')
xlabel('Frequency (Hz)')
ylabel('S (m^2/Hz)')

pause

end