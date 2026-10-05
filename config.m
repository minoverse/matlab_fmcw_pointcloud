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


%% Ground truth target

target_range = 30;      % Target range [m]
target_velocity = 10;   % Target radial velocity [m/s]


%% Derived FMCW parameters

S = B / Tchirp;         % Chirp slope [Hz/s]

Ns = round(fs * Tchirp);   % Number of ADC samples per chirp

range_resolution = c / (2 * B);


%% Range-related parameters

tau = 2 * target_range / c;   % Round-trip delay [s]

fb = S * tau;                 % Beat frequency [Hz]


%% Doppler-related parameters

lambda = c / fc;              % Wavelength [m]

fd = 2 * target_velocity / lambda;   % Doppler frequency [Hz]


%% Print configuration

fprintf("Ns = %d samples\n", Ns);

fprintf("Range resolution = %.3f m\n", ...
    range_resolution);

fprintf("Round trip delay = %.3f ns\n", ...
    tau * 1e9);

fprintf("Expected beat frequency = %.3f MHz\n", ...
    fb / 1e6);

fprintf("Expected Doppler frequency = %.3f kHz\n", ...
    fd / 1e3);