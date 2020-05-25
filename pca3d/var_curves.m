% Plot cumulative variance curves for first 20 modes in each of the 13
% frequency bands

mode_dir = 'data/stage3b';
nfreq = 13;
S = zeros(600,nfreq);

% frequency bands in cpd
freqs = [ 6, 5.5, 5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1, 0.625, 0.3125 ]';
periods = 86400 ./ freqs;

for i = 1:nfreq
    Sfile = sprintf('%s/S_%d', mode_dir, i);
    tmp = load(Sfile);
    n = length(tmp);
    S(1:n,i) = tmp;
end

% compute cumulative sum of variances of each frequency band
cumsum = sum(S.^2);

modevar1 = S(1,:).^2 ./ cumsum;
modevar2 = sum(S(1:2,:).^2) ./ cumsum;
modevar5 = sum(S(1:5,:).^2) ./ cumsum;
modevar10 = sum(S(1:10,:).^2) ./ cumsum;
modevar15 = sum(S(1:15,:).^2) ./ cumsum;
modevar20 = sum(S(1:20,:).^2) ./ cumsum;
modevar30 = sum(S(1:30,:).^2) ./ cumsum;

modevar1 = modevar1(:);
modevar2 = modevar2(:);
modevar5 = modevar5(:);
modevar10 = modevar10(:);
modevar15 = modevar15(:);
modevar20 = modevar20(:);
modevar30 = modevar30(:);

A = [ freqs periods modevar1 modevar2 modevar5 modevar10 modevar15 modevar20 modevar30 ];
