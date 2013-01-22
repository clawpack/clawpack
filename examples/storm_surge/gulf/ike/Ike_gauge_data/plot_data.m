% Simple program to help with loading data
clear
close all

stations = ['RSUVWXYZ'];
L1 = length(stations);

for j  = 1:L1
figure(1)
clf
j

% loads yd_processed, Hs_wave, mean_water, peak_freq
str1 = ['load result_',stations(j)]
eval(str1)

subplot(3,1,1)
plot(yd_processed, Hs_wave, 'k-')
xlabel('Yearday')
ylabel('Hs (m)')

subplot(3,1,2)
plot(yd_processed, mean_water, 'k-')
xlabel('Yearday')
ylabel('Total Depth (m)')

subplot(3,1,3)
plot(yd_processed, 1./peak_freq, 'k-')
xlabel('Yearday')
ylabel('Tp (s)')


pause
end
