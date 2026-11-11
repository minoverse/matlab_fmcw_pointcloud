
function cfg = radar_config()

%% Physical constants
cfg.c = 3e8;

%% FMCW parameters
cfg.fc = 77e9;
cfg.B = 1e9;
cfg.Tchirp = 50e-6;
cfg.fs = 20e6;
cfg.Nchirp = 64;

%% Ground truth
cfg.target_range = [30, 55];
cfg.target_velocity = [10, -5];
cfg.target_angle = [15, -20];
cfg.target_elevation = [10, -8];
cfg.target_amplitude = [1.0, 0.7];

cfg.Ntargets = length(cfg.target_range);

%% Derived parameters
cfg.S = cfg.B / cfg.Tchirp;
cfg.Ns = round(cfg.fs * cfg.Tchirp);

cfg.range_resolution = cfg.c / (2 * cfg.B);

%% Range parameters
cfg.tau = 2 * cfg.target_range / cfg.c;
cfg.fb = cfg.S .* cfg.tau;

%% Doppler parameters
cfg.lambda = cfg.c / cfg.fc;
cfg.fd = 2 * cfg.target_velocity / cfg.lambda;

%% RX array
cfg.Nrx_az = 2;
cfg.Nrx_el = 2;
cfg.Nrx = cfg.Nrx_az * cfg.Nrx_el;

cfg.d_az = cfg.lambda / 2;
cfg.d_el = cfg.lambda / 2;

%% RX positions
cfg.rx_y = [0, cfg.d_az, 0, cfg.d_az];
cfg.rx_z = [0, 0, cfg.d_el, cfg.d_el];

%% Angle FFT
cfg.NangleFFT_az = 128;
cfg.NangleFFT_el = 128;

%% CFAR parameters
cfg.cfar.Tr = 8;
cfg.cfar.Td = 4;
cfg.cfar.Gr = 2;
cfg.cfar.Gd = 1;
cfg.cfar.Pfa = 1e-5;

end
