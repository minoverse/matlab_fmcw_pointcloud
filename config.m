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

target_elevation = [10, -8];

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


%% 2 x 2 planar RX array

Nrx_az = 2;

Nrx_el = 2;

Nrx = Nrx_az * Nrx_el;

d_az = lambda / 2;

d_el = lambda / 2;


%% RX positions

rx_y = [0, d_az, 0, d_az];

rx_z = [0, 0, d_el, d_el];


%% Angle FFT

NangleFFT_az = 128;

NangleFFT_el = 128;