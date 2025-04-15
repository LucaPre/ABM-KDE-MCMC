function [mode] = histmode(theta)
mode=zeros(1,5);
for i=1:5
num_bins = 200;

% Access histogram data
[counts, edges] = histcounts(theta(:,i), num_bins);

% Find the mode
[~, max_idx] = max(counts);
mode_bin_center = (edges(max_idx) + edges(max_idx + 1)) / 2;
mode(i)=mode_bin_center;
end
