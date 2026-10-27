%% Physical constant
c = 3e8;

%% FMCW parameters
fc = 77e9;
B = 1e9;
Tchirp = 50e-6;
fs = 20e6;
Nchirp = 64;

%% Ground truth targets
target_range = [30, 55];
target_velocity = [10, -5];
target_angle = [15, -20];
target_amplitude = [1.0, 0.7];

Ntargets = length(target_range);

%% Derived FMCW parameters
S = B / Tchirp;
Ns = round(fs * Tchirp);
range_resolution = c / (2 * B);

%% Range parameters
tau = 2 * target_range / c;
fb = S .* tau;

%% Doppler parameters
lambda = c / fc;
fd = 2 * target_velocity / lambda;

%% RX array
Nrx = 4;
d_rx = lambda / 2;
NangleFFT = 128;

%% Print configuration
fprintf("Ns = %d samples\n", Ns);
fprintf("Range resolution = %.3f m\n", range_resolution);

for k = 1:Ntargets
    fprintf("\nTarget %d\n", k);
    fprintf("Range = %.2f m\n", target_range(k));
    fprintf("Velocity = %.2f m/s\n", target_velocity(k));
    fprintf("Angle = %.2f deg\n", target_angle(k));
    fprintf("Beat frequency = %.3f MHz\n", fb(k) / 1e6);
    fprintf("Doppler frequency = %.3f kHz\n", fd(k) / 1e3);
end