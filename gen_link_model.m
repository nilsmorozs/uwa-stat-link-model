% This script is used to process the ASUNA Werbellin lake dataset 
% and generate statistical link models as described in the paper

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

% Specify the output MAT file to save the statistical model to
output_mat_file = 'on_off_dur_model_werbellin16.mat';

% Location data mapping (wrong node indexing in data, N1 swapped with N6)
loc_ind = [6, 2, 3, 4, 5, 1];


%% Process the data to obtain the ON/OFF link state durations for every link
on_dur = cell(num_tops, num_nodes, num_nodes);
off_dur = cell(num_tops, num_nodes, num_nodes);
prob_link_fade = NaN(num_tops, num_nodes, num_nodes);

for t = 1:num_tops
    time_ind = ((t-1)*num_pkts_per_top+1) : (t*num_pkts_per_top);
    timeline_mat = [];
    tx_rx_ids = [];
    for n = 1:num_nodes
        for k = 1:num_nodes
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

                % Calculate the probability of link fading
                prob_link_fade(t, n, k) = sum(link_timeline < 0.5) / num_pkts_per_top;
            end
        end
    end
end

%% Generate statistical link models from the data

% Divide links into "Good", "Average", "Bad" using their packet error rate (PER)
per_boundaries = [0.3, 0.6, 1.0];
on_dur_by_per = cell(1, 3);
off_dur_by_per = cell(1, 3);
for t = 1:num_tops
    for n = 1:num_nodes
        for k = 1:num_nodes
            if ~isnan(prob_link_fade(t, n, k))
                link_type = find((prob_link_fade(t, n, k) <= per_boundaries), 1);
                on_dur_by_per{link_type} = [on_dur_by_per{link_type}; on_dur{t, n, k}];
                off_dur_by_per{link_type} = [off_dur_by_per{link_type}; off_dur{t, n, k}];
            end
        end
    end
end

% For each link class plot its own CDF pair of ON/OFF state durations
on_off_dur_model = cell(1, numel(per_boundaries));
for n = 1:numel(per_boundaries)

    % Generate X and Y vectors for the discrete CDFs
    on_dur_sorted = sort(on_dur_by_per{n});
    off_dur_sorted = sort(off_dur_by_per{n});
    prc_vect_on = (1:numel(on_dur_sorted)) ./ numel(on_dur_sorted);
    prc_vect_off = (1:numel(off_dur_sorted)) ./ numel(off_dur_sorted);

    % Fit continuous CDFs to the discrete ON/OFF state duration data
    % ON state
    step_ind = find(diff(on_dur_sorted) > 0.5);
    on_dur_points = (on_dur_sorted(step_ind) + on_dur_sorted((step_ind)+1)) ./ 2;
    on_prc_points = prc_vect_on(step_ind)';
    on_prc_interp_range = linspace(0, 1, 1001);
    on_dur_interp = interp1([0; on_prc_points], [on_dur_sorted(1); on_dur_points], ...
                                  on_prc_interp_range, 'linear');
    on_dur_interp(isnan(on_dur_interp)) = max(on_dur_interp); % some CDF points at >99.9th percentile show NaN [hot fix]
    
    % OFF state
    step_ind = find(diff(off_dur_sorted) > 0.5);
    off_dur_points = (off_dur_sorted(step_ind) + off_dur_sorted((step_ind)+1)) ./ 2;
    off_prc_points = prc_vect_off(step_ind)';
    off_prc_interp_range = linspace(0, 1, 1001);
    off_dur_interp = interp1([0; off_prc_points], [off_dur_sorted(1); off_dur_points], ...
                                  off_prc_interp_range, 'linear');
    off_dur_interp(isnan(off_dur_interp)) = max(off_dur_interp); % some CDF points at >99.9th percentile show NaN [hot fix]

    % Save these CDFs in the structure as the model for this type link
    on_off_dur_model{n}.on_dur_vals = on_dur_interp;
    on_off_dur_model{n}.on_prc_vals = on_prc_interp_range;
    on_off_dur_model{n}.off_dur_vals = off_dur_interp;
    on_off_dur_model{n}.off_prc_vals = off_prc_interp_range;

    % Plot empirical CDFs and interpolated continuous CDFs for this type of link
    figure; hold on;
    cdfplot(on_dur_by_per{n});
    cdfplot(off_dur_by_per{n});
    plot(on_dur_interp, on_prc_interp_range, 'b-.', 'linewidth', 1.5);
    plot(off_dur_interp, off_prc_interp_range, 'r:', 'linewidth', 1.5);
    legend('ON state: empirical CDF (sea trial data)', 'OFF state: empirical CDF (sea trial data)', ...
            'ON state: interpolated continuous CDF',  'OFF state: interpolated continuous CDF', 'Location', 'SouthEast');
    xlabel('Link state duration, sec');
    legend('boxoff')
    ylabel('CDF');
    title('');
    box on;
    grid on;
    xlim([0 Inf]);
end

% Save the ON-OFF duration model as a MAT file
save(output_mat_file, 'on_off_dur_model');

