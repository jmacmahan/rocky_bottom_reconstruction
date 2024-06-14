clear all
close all
clc


%% Apply calculation to the measured mean spectral results
% Measured Mean Gzz
data = [0.0021  235.4612
    0.0024  190.7869
    0.0027  167.0137
    0.0030  134.7467
    0.0034  101.4884
    0.0038   86.8570
    0.0043   75.9809
    0.0049   62.5006
    0.0055   55.1436
    0.0061   48.7867
    0.0069   41.8597
    0.0078   34.7688
    0.0087   31.1040
    0.0098   26.0298
    0.0111   22.0649
    0.0125   18.3490
    0.0141   16.2768
    0.0158   13.1796
    0.0178   11.1871
    0.0200    9.5240
    0.0225    7.8055
    0.0254    6.3853
    0.0285    5.4089
    0.0321    4.4193
    0.0361    3.5707
    0.0407    2.9037
    0.0458    2.3132
    0.0515    1.8247
    0.0580    1.4315
    0.0653    1.0940
    0.0734    0.8365
    0.0826    0.6200
    0.0930    0.4589
    0.1047    0.3319
    0.1178    0.2368
    0.1325    0.1687
    0.1491    0.1200
    0.1678    0.0879
    0.1889    0.0662];

% Define columns
k_m = data(:,1);
Gzz_m = data(:,2);
clear data

% Desired frequency resolution
dk = 1/750; % This defines the spacing in wavenumber (1/meters)

% Interpolate the 1-D slope spectrum
max_k = 0.5;
k = dk:dk:max_k;
Gzz = interp1(log10(k_m), log10(Gzz_m), log10(k), 'linear', 'extrap');
Gzz = 10.^Gzz;

% Normalize 1D spectrum
Gzz = Gzz ./ trapz(k, Gzz);
Gzz_meas = Gzz;

% Plot spectrum
subplot(231)
loglog(k, Gzz, 'ro', 'LineWidth', 1, 'MarkerFaceColor', 'r')
hold on

% Convert variance to amplitude
Azz = sqrt(2 .* (Gzz .* dk));
Azz = [0 Azz];
Azz_meas = Azz;
k = [0 k];

% Generate random phase spectrum
N = length(Azz);
random_phases_1d = exp(1i * 2 * pi * rand(1, N));

% Generate FFT spectrum
k_2side = [k -fliplr(k(2:end-1))];
Azz_2side = [Azz fliplr(Azz(2:end-1))] ./ 2;
phases_2side = [random_phases_1d conj(fliplr(random_phases_1d(2:end-1)))];
complex_spectrum = Azz_2side .* phases_2side;

% Perform inverse FFT to get the spatial domain profile
profile1 = ifft(complex_spectrum .* (2 .* N));
profile1 = real(profile1); % Take the real part

% Scale the profile to the physical length of 750 meters
x = linspace(0, 750, 750);

% Plot profile
profile1 = profile1 ./ max(abs(profile1));
subplot(232)
hold on
plot(x, profile1, 'r-');

% Create the 2-D isotropic spectrum from the interpolated 1-D spectrum
M = 2 * (length(Azz) - 1);
[fx, fy] = meshgrid(linspace(-M/2, M/2-1, M) / 750, linspace(-M/2, M/2-1, M) / 750);
radius = sqrt(fx.^2 + fy.^2);

% Interpolate the isotropic spectrum to match the radius values
radius_flat = radius(:);
radius_flat(radius_flat > max(k)) = max(k);
Azz_2d = interp1(k, Azz, radius_flat, 'linear', 'extrap');
Azz_2d = reshape(Azz_2d, M, M);

% Normalize 2D spectrum
Azz_2d = Azz_2d ./ sum(sum(Azz_2d));

% Define wavenumber axes
kn = [-fliplr(k(2:end)) k(1:end-1)];

% Generate random phase spectrum
random_phases_2d = exp(1i * 2 * pi * rand(M, M));

% Create the complex spectrum (amplitude * phase)
complex_spectrum = Azz_2d .* random_phases_2d;

% Perform inverse FFT to get the spatial domain surface
surface = ifft2(ifftshift(complex_spectrum));
surface = real(surface); % Take the real part

% Normalize the surface to have zero mean and unit variance
surface1 = surface ./ max(abs(surface(:)));

% Scale the surface to the physical length of 750 meters
y = linspace(0, 750, 750);

% Display the generated surface
subplot(223)
pcolor(x, y, surface1); shading interp
colormap(cmocean('balance', 'pivot', 0))
xlabel('$x~\mathrm{(m)}$', 'Interpreter', 'latex');
ylabel('$y~\mathrm{(m)}$', 'Interpreter', 'latex');
title('Measurements')
text(25,40,'d)','FontSize',36)
%% Computations based on the simple spectral model
% Given 1-D slope spectrum data
k_i = [0.0013, 0.012, 0.1, 0.5]; % in 1/meters
Gzz_i = [330.7617, 19.0108, 0.3854, 0.0046]; % slope spectrum values

% Interpolate the 1-D slope spectrum
max_k = 0.5;
k = dk:dk:max_k;
Gzz = interp1(log10(k_i), log10(Gzz_i), log10(k), 'linear', 'extrap');
Gzz = 10.^Gzz;

% Normalize 1D spectrum
Gzz = Gzz ./ trapz(k, Gzz);
Gzz_model = Gzz;

% Plot spectrum
subplot(231)
hold on
loglog(k, Gzz, 'b')
xlabel('$k~\mathrm{(m^{-1})}$', 'Interpreter', 'latex')
ylabel('$\overline{G}_{z^\prime}~\mathrm{(m)}$', 'Interpreter', 'latex')
legend('Measurements', 'Model', 'Location', 'best')
text(2*10^-3,.01,'a)','FontSize',36)

% Convert variance to amplitude
Azz = sqrt(2 .* (Gzz .* dk));
Azz = [0 Azz];
Azz_model = Azz;
k = [0 k];

% Generate FFT spectrum
phases_2side = [random_phases_1d conj(fliplr(random_phases_1d(2:end-1)))];
complex_spectrum = Azz_2side .* phases_2side;

% Perform inverse FFT to get the spatial domain profile
profile2 = ifft(complex_spectrum .* (2 .* N));
profile2 = real(profile2); % Take the real part

profile2 = profile2 ./ max(abs(profile2));
subplot(232)
hold on
plot(x, profile2, 'b');
xlabel('$x~\mathrm{(m)}$', 'Interpreter', 'latex');
ylabel('$z^\prime/z_{max}^\prime$', 'Interpreter', 'latex')
ylim([-1 1])
xlim([0 750])
text(10,-.9,'b)','FontSize',36)
% Interpolate the isotropic spectrum to match the radius values
radius_flat = radius(:);
radius_flat(radius_flat > max(k)) = max(k);
Azz_2d = interp1(k, Azz, radius_flat, 'linear', 'extrap');
Azz_2d = reshape(Azz_2d, M, M);

% Normalize 2D spectrum
Azz_2d = Azz_2d ./ sum(sum(Azz_2d));

subplot(233)
pcolor(kn, kn, log10(Azz_2d)); shading interp
xlabel('$k_{x}~\mathrm{(m^{-1})}$', 'Interpreter', 'latex');
ylabel('$k_{y}~\mathrm{(m^{-1})}$', 'Interpreter', 'latex');
text(-.4,-.4,'c)','FontSize',36,'Color','w')
% Create the complex spectrum (amplitude * phase)
complex_spectrum = Azz_2d .* random_phases_2d;

% Perform inverse FFT to get the spatial domain surface
surface = ifft2(ifftshift(complex_spectrum));
surface = real(surface); % Take the real part

% Normalize the surface to have zero mean and unit variance
surface2 = surface ./ max(abs(surface(:)));

% Display the generated surface
subplot(224)
pcolor(x, y, surface2); shading interp
colormap(cmocean('balance', 'pivot', 0))
xlabel('$x~\mathrm{(m)}$', 'Interpreter', 'latex');
ylabel('$y~\mathrm{(m)}$', 'Interpreter', 'latex');
title('Model')
text(25,40,'e)','FontSize',36)
% Calculate RMSE for profiles
RMSE_profiles = sqrt(mean((profile1 - profile2).^2))

% Calculate RMSE for surfaces
RMSE_surfaces = sqrt(mean((surface1(:) - surface2(:)).^2))
