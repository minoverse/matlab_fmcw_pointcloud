%% User-defined parameters
% fc, B, Tchirp, fs, Nchirp
% target_range, target_velocity

%% Derived parameters
% S, Ns, tau, fb
% range_resolution, lambda, fd


%% Physical constant

c = 3e8;   % Speed of light [m/s]


%% FMCW parameters

fc = 77e9;          % Carrier frequency [Hz]
B = 1e9;            % Bandwidth [Hz]
Tchirp = 50e-6;     % Chirp duration [s]

fs = 20e6;          % ADC sampling rate [Hz]
Nchirp = 64;        % Number of chirps


%% Ground truth targets

target_range = [30, 55];          % Target ranges [m]
target_velocity = [10, -5];       % Target radial velocities [m/s]
target_amplitude = [1.0, 0.7];    % Relative amplitudes

Ntargets = length(target_range);


%% Derived FMCW parameters

S = B / Tchirp;

Ns = round(fs * Tchirp);

range_resolution = c / (2 * B);


%% Range-related parameters

tau = 2 * target_range / c;

fb = S .* tau;


%% Doppler-related parameters

lambda = c / fc;

fd = 2 * target_velocity / lambda;


%% Print configuration

fprintf("Ns = %d samples\n", Ns);
fprintf("Range resolution = %.3f m\n", range_resolution);

for k = 1:Ntargets

    fprintf("\nTarget %d\n", k);

    fprintf("Range = %.2f m\n", target_range(k));

    fprintf("Velocity = %.2f m/s\n", target_velocity(k));

    fprintf("Beat frequency = %.3f MHz\n", ...
        fb(k) / 1e6);

    fprintf("Doppler frequency = %.3f kHz\n", ...
        fd(k) / 1e3);

end