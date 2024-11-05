% This script is used to analyse the data from the Werbellin lake experiment publicly available in the ASUNA database

% Copyright 2024 Nils Morozs, University of York (nils.morozs@york.ac.uk)
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
% associated documentation files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial
% portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
% NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
% OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Load the data
load('BerlinMultimodal-06-16.mat');
modem_ind = 1; % Evologics modem frequency band 
               % (1 - low frequency, 2 - medium frequency, 3 - high frequency)
               % The most data was available for the low frequency modems (18-34 kHz)
               % The data from the low frequency modems was used used for generating statistical models in the paper
num_tops = 5;  % number of topologies
num_pkts_per_top = 600; % number of transmitted packets per topology
num_nodes = size(FullTopMat, 2);

% Location data mapping (wrong node indexing in the dataset, N1 should be swapped with N6)
loc_ind = [6, 2, 3, 4, 5, 1];


%% Process the data to obtain the ON/OFF link state durations for every link
on_dur = cell(num_tops, num_nodes, num_nodes);
off_dur = cell(num_tops, num_nodes, num_nodes);
prob_link_fade = NaN(num_tops, num_nodes, num_nodes);
dist_mat = NaN(num_tops, num_nodes, num_nodes);
link_pair_corr = [];
link_pair_id_mat = []; % t, tx1, rx1, tx2, rx2

% Loop through evey topology and transmitter-receiver pair
for t = 1:num_tops
    time_ind = ((t-1)*num_pkts_per_top+1) : (t*num_pkts_per_top);
    timeline_mat = [];
    tx_rx_ids = [];
    for n = 1:num_nodes
        for k = 1:num_nodes
            % If there is data available for this transmitter-receiver link
            if (n ~= k) && any(squeeze(FullTopMat(time_ind, n, k, modem_ind)))
                
                % Carve out the link timeline for this link and topology
                link_timeline = squeeze(FullTopMat(time_ind, n, k, modem_ind));
                timeline_mat = [timeline_mat, link_timeline];
                tx_rx_ids = [tx_rx_ids; n, k];

                % Determine the link state change times and initial state
                state_change_times = find(diff(link_timeline) ~= 0);
                init_state = link_timeline(1);
                
                % Depending on the initial link state, store an array of link state durations
                if (init_state == 0)
                    on_start_times = state_change_times(1:2:end);
                    off_start_times = state_change_times(2:2:end);
                    on_dur{t, n, k} = off_start_times - on_start_times(1:numel(off_start_times));
                    off_dur{t, n, k} = on_start_times(2:end) - off_start_times(1:numel(on_start_times)-1);
                else
                    off_start_times = state_change_times(1:2:end);
                    on_start_times = state_change_times(2:2:end);
                    on_dur{t, n, k} = off_start_times(2:end) - on_start_times(1:numel(off_start_times)-1);
                    off_dur{t, n, k} = on_start_times - off_start_times(1:numel(on_start_times));
                end

                % Calculate the probability of link fading (same as packet error rate)
                prob_link_fade(t, n, k) = sum(link_timeline < 0.5) / num_pkts_per_top;

                % Note the distance between nodes
                earth_radius = 6371e3; % m;
                coord1 = squeeze(FullLocMat(time_ind(1), loc_ind(n), :));
                coord2 = squeeze(FullLocMat(time_ind(1), loc_ind(k), :));
                lat1 = coord1(2)*pi/180;
                lon1 = coord1(1)*pi/180;
                lat2 = coord2(2)*pi/180;
                lon2 = coord2(1)*pi/180;
                dist_mat(t, n, k) = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*earth_radius;
            end
        end
    end

    % Examine correlation between the timelines
    cc = corrcoef(timeline_mat);
    % Store all of them in the big matrix
    for i = 1:size(cc, 2)-1
        for j = i+1:size(cc, 2)
            link_pair_corr = [link_pair_corr; cc(i, j)];
            link_pair_id_mat = [link_pair_id_mat; t, tx_rx_ids(i, :), tx_rx_ids(j, :)];
        end
    end

end

%% Plot an example link timeline

% Specify the link to be plotted
top_ind = 2;
node_ind = [2, 5];
% The three links included in FIg. 2 in the paper
%   modem_ind = 1; top_ind = 2; node_ind = [2, 5];
%   modem_ind = 1; top_ind = 2; node_ind = [4, 2];
%   modem_ind = 2; top_ind = 1; node_ind = [3, 5];

figure; hold on;
link_timeline = squeeze(FullTopMat((1:num_pkts_per_top) + num_pkts_per_top*(top_ind-1),...
                                                       node_ind(1), node_ind(2), modem_ind));
plot(1:num_pkts_per_top, link_timeline, 'b-')
xlabel('Time, sec');
ylabel('Link state');
box on;
yticks([0 1])
xticks(0:100:num_pkts_per_top)
grid on;
yticklabels({'OFF', 'ON'})
ylim([-0.1, 1.1])

%% Plot probability of packet success vs link distance
figure; hold on;
markers = {'bo', 'r^', 'gx'};
plot(dist_mat(:), 1-prob_link_fade(:), markers{modem_ind});
xlabel('Link distance, m');
ylabel('Probability of packet success');
box on;
grid on;

%% Plot variance between pairs of link timelines using different grouping criteria

% Establish which cross correlation coeffs are for reverse links (same two nodes)
reverse_links = (link_pair_id_mat(:, 2) == link_pair_id_mat(:, 5)) & (link_pair_id_mat(:, 3) == link_pair_id_mat(:, 4));
same_rx_links = (link_pair_id_mat(:, 3) == link_pair_id_mat(:, 5));
same_tx_links = (link_pair_id_mat(:, 2) == link_pair_id_mat(:, 4));
other_links = ~ (reverse_links | same_rx_links | same_tx_links);
% Note distance between receiver for every pair of links
rx2rx_dist = NaN(size(link_pair_id_mat, 1), 1);
for n = 1:size(link_pair_id_mat, 1)
    rx2rx_dist(n) = dist_mat(link_pair_id_mat(n, 1), link_pair_id_mat(n, 3), link_pair_id_mat(n, 5));
end

figure; hold on;
plot(zeros(1, sum(same_rx_links)), link_pair_corr(same_rx_links).^2, 'bo');
plot(rx2rx_dist(same_tx_links), link_pair_corr(same_tx_links).^2, 'g^');
plot(rx2rx_dist(reverse_links), link_pair_corr(reverse_links).^2, 'mx');
plot(rx2rx_dist(other_links), link_pair_corr(other_links).^2, 'rv');
box on; grid on;
legend('Links with same Rx', 'Links with same Tx', 'Reverse links (same Tx <-> Rx)', 'Other pairs of links');
xlabel('Distance between receivers, m');
ylabel('Variance (R^2) between ON/OFF patterns of two links');

