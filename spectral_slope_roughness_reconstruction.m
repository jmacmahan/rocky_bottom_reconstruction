clear all
close all
clc

%% Define constants and parameters
dk = 1/100;                   % Desired frequency resolution
kk = dk:dk:1;                 % Wavenumber vector
slope = fliplr(-3:0.25:-1);   % Slope values for different spectra

x = 1:0.25:100;               % Spatial domain
z = zeros(size(x));           % Initialize z

haxs = makeSubPlots(0.1, 0.1, 0.02, 0.1, 0.1, 0.06, 3, 3); % Create subplots

%% Loop over different slope values
for j = 1:length(slope)
    Gzz = 10.^(slope(j) * log10(kk));  % Compute Gzz based on slope
    Gzz = Gzz ./ trapz(kk, Gzz);       % Normalize Gzz
    var_z(j) = trapz(kk, Gzz);         % Compute variance
    
    % Generate random signal
    for i = 1:length(kk)
        varx = Gzz(i) * dk;
        A = sqrt(2 * varx);
        phase = rand(1) * 2 * pi;
        z = z + A * cos(2 * pi * kk(i) * x - phase);
    end
    z = z / std(z);  % Normalize z

    % Plotting
    figure(1)
    zz = diff(z) ./ diff(x);  % Compute slope
    percent_flat = length(find(abs(zz) <= 1/100)) / length(zz) * 100;  % Compute flatness
    tt = "$\alpha=~$" + num2str(-slope(j));  % Title with slope value
    tt2 = "$\sigma_{\beta}=" + num2str(round(std(zz) * 10) / 10) + "; " + num2str(round(percent_flat * 100) / 100) + "\% ~\mathrm{flat}$";  % Title with slope statistics

    subplot(haxs(j))
    if j == 1
        Anorm = max(abs(z));  % Normalize amplitude for first plot
    end
    plot(x / x(end), z / Anorm)
    title(tt, 'Interpreter', 'latex')
    ylim([-1 1])
    if ismember(j, [1, 4, 7])
        ylabel('$z{\prime}/z{\prime}_{max}$', 'Interpreter', 'latex')
    end
    if j >= 7
        xlabel('$L/L_{max}$', 'Interpreter', 'latex')
    end
    if j >= 1 && j <= 6
        set(gca, 'XTickLabel', [])
    end
    if ismember(j, [2, 3, 5, 6, 8, 9])
        set(gca, 'YTickLabel', [])
    end
    text(0.05, -0.85, tt2, 'FontSize', 25, 'Interpreter', 'latex')
end
