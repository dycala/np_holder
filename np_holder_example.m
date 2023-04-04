path = "C:\Users\dylan\NN MATLAB SCRIPTS\2022-05-17_17-38-45\experiment1\recording2";

% load data
data=np_holder(path,'quality_threshold',["good"]); % example of how you could load other quality cells

% bin data
data.bin_spikes('bin_size',0.01) % default is 10 ms bins and only good quality cells

% smooth data
data.smooth_binned % default is 20 ms kernel (sd), sampled at 10 ms, and only good quality cells
data.smooth_binned('kernel_size',0.01,'sample_time',0.01)

% get time of TTL
data.get_event_times(100)

% get cspk
data.collect_cspk()
data.select_cspk()
