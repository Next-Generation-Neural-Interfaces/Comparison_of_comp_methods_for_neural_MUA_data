%% Testing encoding and decoding
% Where we encode the data according to the method, look at the BR, then
% decode it. Look at a sample of data and parameters.

% NOTE: be wary of MATLAB indexing, it starts at 1 whereas in the paper,
% i=0 is a valid index for MUA events.


    hh = tic;
    clearvars -except all_data hh index
    fprintf(['First clearvars time: ',num2str(toc(hh)),'\n'])
    rng(1);

    % Various BPs, channel counts, saturate at S
    timer = tic;

    % Look at different BP; and number of channels... add recordings.
    BP_vec = [1,5,10,20,50,100];
    channel_vec = [10, 100, 1000, 10000, 50000];% 10000 50000];
    S_vec = 2:15;
    hist_vec = [0 2:6];
    factor_vector = [1e-3 1e-2 1e-1 1 10 100];
    clip_length = 100; % in s


    %% Chosen index
    index = 5;

    %% Index desired channel nb, BP and S
    chan_index = mod(index,length(channel_vec));
    count_1 = floor(index/length(channel_vec));

    BP_index = mod(count_1,length(BP_vec));
    count_2 = floor(count_1/length(BP_vec));

    S_index = mod(count_2,length(S_vec));
    count_3 = floor(count_2/length(S_vec));
    
    hist_index = mod(count_3,length(hist_vec));

    nb_channels = channel_vec(chan_index+1);
    BP = BP_vec(BP_index+1);
    S = S_vec(S_index+1);
    hist_size = hist_vec(hist_index+1);
    % We overwrite indexing here, just for testing
    hist_size = 6;
    BP = 50;
    S = 2;
    nb_channels = 10;

    fprintf(['\nChannels: ',num2str(nb_channels),'; BP: ',num2str(BP),'; S: ',num2str(S),'; Hist size: ',num2str(hist_size),' bits \n \n'])
% 
    data_path = 'D:\Dropbox (Imperial NGNI)\NGNI Share\Workspace\Oscar\Work\Async MUA compression';
    save_file = [data_path,'\added_AH_and_E_for_delta_10\all_results_chan_', num2str(nb_channels),'_BP_',num2str(BP),'_S_',num2str(S),'.mat'];
%     data_path = '/rds/general/user/ows18/home/Async_MUA';
%     save_file = [data_path,'/added_AH_and_E_for_delta_10/all_results_chan_', num2str(nb_channels),'_BP_',num2str(BP),'_S_',num2str(S),'_hist_',num2str(hist_size),'.mat'];
    if isfile(save_file)
      done = [];
      fprintf('Results already in\n')
      rand(15) % makes file size over 1 kb, easy to see in HPC that the job succeeded.
      return
    end
    
%     % Load data
    hh = tic;
%     load([data_path,'/all_data_1ms.mat'])
    load([data_path,'\all_data_1ms.mat'])
    fprintf(['Load data time: ',num2str(toc(hh)),'\n'])

    signal_length_s = length(all_data(:,1)) / 1000;
    max_nb_channels = length(all_data(1,:));


    %% Bin data
    hist_channel_vector = 0.5:1:nb_channels+0.5;
    data_chans = randsample(length(all_data(1,:)),nb_channels);
    data = all_data(:,data_chans); % take a random sample of channels

    % Find where firing rates above 0 occur
    hh = tic;
    [a1,a2] = find(data~=0);
    fprintf(['Orig find events time: ',num2str(toc(hh)),'\n'])
    
    % Calibration histogram
    if hist_size == 0
        time_vector_post_cal = 0.5:BP:length(data(:,1))+0.5;
        mapping = 0;
        inv_mapping = 0;
        
    else % If we do the histogram for mapping
        
        % Calibration/mapping
        hh = tic;
        time_vector_cal = 0.5:BP:BP*2^hist_size+0.5; % 2^hist_size samples for the calibration histogram
        h = histogram2(a1,a2,time_vector_cal,hist_channel_vector);
        BP_MUA = h.Values;
        clear h time_vector_cal
        fprintf(['Sample hist time: ',num2str(toc(hh)),'\n'])
        
        % Saturate data
        max_val = S-1;
        BP_MUA(BP_MUA>max_val) = max_val;
         
        % Create mapping
        hh = tic;
        mapping = uint16(zeros(nb_channels,S));
        for channel = 1:nb_channels
            hist = histogram(BP_MUA(:,channel),-0.5:1:S-0.5);
            sample_hist = hist.Values;
            mapping(channel,:) = zheng_codewordAssign_flip(sample_hist); % used for mapping
        end
        peak_idx = mapping(:,1);
        
        % Inverse mapping
        inv_mapping = uint16(zeros(nb_channels,S));
        for j = 1:nb_channels
            for e = 1:length(mapping(j,:))
                inv_mapping(j,e) = find(mapping(j,:)==e);
            end
        end

        fprintf(['Create mapping time: ',num2str(toc(hh)),'\n'])

        clear BP_MUA hist
        time_vector_post_cal = BP*2^hist_size+0.5:BP:length(data(:,1))+0.5;
    end

    % Post-cal data
    hh = tic;
    h = histogram2(a1,a2,time_vector_post_cal,hist_channel_vector);
    BP_MUA = h.Values;
    clear h time_vector_post_cal a1 a2 data %all_data
    fprintf(['To-be-comp. data histogram time: ',num2str(toc(hh)),'\n'])

    % Saturate data
    max_val = S-1;
    BP_MUA(BP_MUA>max_val) = max_val;
    data_length = length(BP_MUA(:,1))*nb_channels;

    %% Basic implementation statistics
    % Sync codeword length:
    basic_sync_codeword_len = ceil(log2(S)); % m for sync

    % Sync pre-trianed static codeword lengths
    lossless_basic_sync_MUA_codeword_lengths = [1:S-1, S-1];
            
    % Explicit async binary codeword, nb of MUA events
    basic_exp_async_nb_MUA_events_codeword_len = ceil(log2(S-1)); % m for async
    lossless_basic_exp_async_MUA_codeword_lengths = [1:S-2, S-2];

    % Explicit async basic channel ID codeword
    basic_exp_async_chan_ID_codeword_len = ceil(log2(nb_channels)); % k1

    % Codeword representing channel ID (and stop symbol in basic
    % group async)
    basic_group_asycn_codeword_len = ceil(log2(nb_channels+S-2)); % k2

    %% Probs
    hh = tic;
    % Prob an event i will occur on channel j
    p = zeros(max_val+1,nb_channels);
    for j = 1:nb_channels  % for each channel
        for i = 0:max_val
            p(i+1,:) = sum(BP_MUA==i)/length(BP_MUA(:,1));
        end
    end
    
    % Average prob an i event will occur across channels
    p_channel_mean = mean(p,2)./sum(mean(p,2));
    

    % Async probs
    % Prob of sending out a channel ID (odds of it having more
    % than 0 events). Used in both async methods.
    % i.e. p_send
    p_send_chan = zeros(nb_channels,1);
    for j = 1:nb_channels  % for each channel
        p_send_chan(j,1) = 1-p(1,j);
    end
    % Relative prob of sending out a channel ID
    % i.e. p_send_norm. Used in both async methods.
    rel_p_send_chan = zeros(nb_channels,1);
    for j = 1:nb_channels  % for each channel
        rel_p_send_chan(j,1) = p_send_chan(j,1)/sum(p_send_chan);
    end

    % Prob an event i>0 will occur on a channel
    % i.e. p_norm. We use this in the entropy calculation for
    % explicit async.
    rel_p_positive_event = zeros(max_val+1,nb_channels);
    for i = 1:max_val
        rel_p_positive_event(i+1,:) = p(i+1,:)./sum(p(2:end,:),1);
    end
    
    % Average prob an event i>0 will occur across all channels
    % i.e. pe. Importantly, this isn't the channel-average of
    % hat(p). We have to normalise by each bin i > 0 before we take the
    % channel average. Used in Adaptive Huffman of exp async for firing
    % rate.
    rel_p_positive_event_mean_all_channels = zeros(max_val+1,1);
    temp = 0;
    for i = 1:max_val
        temp = temp + sum(p(i+1,:));
    end
    for i = 1:max_val
        rel_p_positive_event_mean_all_channels(i+1,1) = sum(p(i+1,:))/temp;
    end
    temp = [];


    % Prob of stop symbol, that atleast 1 channel will have i MUA
    % events occur on it. Used for stop symbols in group´async.
    % i.e. ps
    % NOTE: CAN FUCKING CHANGE. BUT, NOT RELATIVE TO P(0), BUT WHATEVER IS
    % THE MOST COMMON FIRING RATE FOR THE CHANNEL. SO STOP SYMBOL IS NOT
    % EQUAL TO A SPECIFIC FIRING RATE, BUT THE NTH MOST COMMON FIRING RATE
    % FOR EACH CHANNEL. BUUUT, no histogram for group async anyway...
    p_atleast_one_channel_has_i_events = zeros(max_val+1,1);
    for i = 2:max_val % start at 2 since no stop symbol for i=0,1
        p_atleast_one_channel_has_i_events(i+1,1) = 1-prod(1-p(i+1,:));
    end

    % NOTE: the exp- group Huffman dicts should be trained on the relative
    % prob of channel IDs AND stop symbols, since the pairing is
    % unknown: unlike with explicit async, you don't know that the
    % nb of events will follow the channel ID. In this case, you
    % may have x channel IDs, then y stop symbols, no idea what the
    % pairing is. The probs need to be relative to eachother, so
    % it's 1 encoder.            
    combined_prob_group = [p_send_chan; p_atleast_one_channel_has_i_events(3:end)]; % index of 3 since no codeword for i=0,1.
    rel_combined_prob_group = combined_prob_group./sum(combined_prob_group);

    % Mapped version: for working on mapping, firing rate stats only. We map the data
    % based on sample histogram.
    if hist_size ~= 0
        p_mapped = p; 
        for channel = 1:nb_channels
            p_mapped(:,channel) = p_mapped(mapping(channel,:),channel);
        end
        
        % Average prob an i event will occur across channels (mapped probs)
        p_channel_mean_mapped = mean(p_mapped,2)./sum(mean(p_mapped,2));
        
        % Async probs
        % Prob of sending out a channel ID (odds of it having more
        % than the most common FR on that channel). Used in both async methods.
        % i.e. p_send
        p_send_chan_mapped = zeros(nb_channels,1);
        for j = 1:nb_channels  % for each channel
            p_send_chan_mapped(j,1) = 1-p_mapped(1,j);
        end
        
        % Relative prob of sending out a channel ID
        % i.e. p_send_norm. Used in both async methods.
        rel_p_send_chan_mapped = zeros(nb_channels,1);
        for j = 1:nb_channels  % for each channel
            rel_p_send_chan_mapped(j,1) = p_send_chan_mapped(j,1)/sum(p_send_chan_mapped);
        end
        % Average prob an event i>0 will occur across all channels
        % i.e. pe. Importantly, this isn't the channel-average of
        % hat(p). We have to normalise by each bin i > 0 before we take the
        % channel average. USed in Adaptive Huffman. MAPPED VERSION.
        %
        % KEY POINT: we don't use p(i==0) as the base, but we use the most common firing
        % rate. It is included in the mapping information we send off-implant.
        rel_p_positive_event_mean_all_channels_mapped = zeros(max_val+1,1);
        temp = 0;
        for i = 1:max_val
            temp = temp + sum(p_mapped(i+1,:));
        end
        for i = 1:max_val % normalize to get a relative probility vector
            rel_p_positive_event_mean_all_channels_mapped(i+1,1) = sum(p_mapped(i+1,:))/temp;
        end
        temp = [];
        
        p_atleast_one_channel_has_i_events_mapped = zeros(max_val+1,1);
        for i = 2:max_val % start at 2 since no stop symbol for i=0,1
            p_atleast_one_channel_has_i_events_mapped(i+1,1) = 1-prod(1-p_mapped(i+1,:));
        end
        
        % NOTE: the exp- group Huffman dicts should be trained on the relative
        % prob of channel IDs AND stop symbols, since the pairing is
        % unknown: unlike with explicit async, you don't know that the
        % nb of events will follow the channel ID. In this case, you
        % may have x channel IDs, then y stop symbols, no idea what the
        % pairing is. The probs need to be relative to eachother, so
        % it's 1 encoder.            
        combined_prob_group_mapped = [p_send_chan_mapped; p_atleast_one_channel_has_i_events_mapped(3:end)]; % index of 3 since no codeword for i=0,1.
        rel_combined_prob_group_mapped = combined_prob_group_mapped./sum(combined_prob_group_mapped);

    end

    fprintf(['Calculate base probs time: ',num2str(toc(hh)),'\n'])
    
    %% Basic implementation BRs  
    hh = tic;

    % Encode and decode each, make sure they're decodable and get the BR
    fprintf('Checking basic implementations\n')
    [check_B_sync,B_sync] = check_sync_basic(BP_MUA,basic_sync_codeword_len);
    [check_B_exp_async,B_exp_async] = check_exp_async_basic(BP_MUA, basic_exp_async_chan_ID_codeword_len, basic_exp_async_nb_MUA_events_codeword_len);
    [check_B_group_async,B_group_async] = check_group_async_basic(BP_MUA, basic_group_asycn_codeword_len, S);

    fprintf(['Basic implementations time: ',num2str(toc(hh)),'\n'])
    
    %% Static Huffman versions
    hh = tic;
    [check_SH_sync,SH_sync] = check_sync_SH(BP_MUA, mapping, inv_mapping, S);
    [check_SH_exp_async,SH_exp_async] = check_exp_async_SH(BP_MUA, mapping, inv_mapping, basic_exp_async_chan_ID_codeword_len,S);
    [check_SH_delta_async,SH_delta_async] = check_delta_async_SH(BP_MUA, mapping, inv_mapping, basic_exp_async_chan_ID_codeword_len, S);
    
    fprintf(['SH implementations time: ',num2str(toc(hh)),'\n'])

    
    %% Adaptive Huffman encoder versions
    hh = tic;
    if hist_size == 0
        [check_AH_sync,AH_sync] = check_sync_AH(BP_MUA, mapping, inv_mapping, S, p_channel_mean);
        [check_AH_exp_async,AH_exp_async] = check_exp_async_AH(BP_MUA, mapping, inv_mapping, S, rel_p_send_chan, rel_p_positive_event_mean_all_channels);
        [check_AH_delta_async,AH_delta_async] = check_delta_async_AH(BP_MUA, mapping, inv_mapping, S);
        [check_AH_group_async,AH_group_async] = check_group_async_AH(BP_MUA, mapping, inv_mapping, S,rel_combined_prob_group); % we use idx instead if idx_fr_above_0 since we map the BP_MUA data directly
    else
        [check_AH_sync,AH_sync] = check_sync_AH(BP_MUA, mapping, inv_mapping, S, p_channel_mean_mapped);
        [check_AH_exp_async,AH_exp_async] = check_exp_async_AH(BP_MUA, mapping, inv_mapping, S, rel_p_send_chan_mapped,rel_p_positive_event_mean_all_channels_mapped);
        [check_AH_delta_async,AH_delta_async] = check_delta_async_AH(BP_MUA, mapping, inv_mapping, S);
        [check_AH_group_async,AH_group_async] = check_group_async_AH(BP_MUA, mapping, inv_mapping, S, rel_combined_prob_group_mapped); % we use idx instead if idx_fr_above_0 since we map the BP_MUA data directly
    end

    fprintf(['AH implementations time: ',num2str(toc(hh)),'\n'])

    % Turn average codeword lengths into BR by dividing by BP
    B_sync = B_sync./(BP/1000);
    B_exp_async = B_exp_async./(BP/1000);
    B_group_async = B_group_async./(BP/1000);
    
    SH_sync = SH_sync./(BP/1000);
    SH_exp_async = SH_exp_async./(BP/1000);
    SH_delta_async = SH_delta_async./(BP/1000);
    
    AH_sync = AH_sync./(BP/1000);
    AH_exp_async = AH_exp_async./(BP/1000);
    AH_delta_async = AH_delta_async./(BP/1000);
    AH_group_async = AH_group_async./(BP/1000);

    all_results = [B_sync; B_exp_async; B_group_async;...
                   SH_sync; SH_exp_async; SH_delta_async;...
                   AH_sync; AH_exp_async; AH_delta_async; AH_group_async];

    
    fprintf('\nBRs given in bits/s/channel:\n')
    
    fprintf(['Basic sync BR: ',num2str(B_sync),'\n'])
    fprintf(['Basic exp async BR: ',num2str(B_exp_async),'\n'])
    fprintf(['Basic group async BR: ',num2str(B_group_async),'\n \n'])
    
    fprintf(['SH sync BR: ',num2str(SH_sync),'\n'])
    fprintf(['SH exp async BR: ',num2str(SH_exp_async),'\n'])
    fprintf(['SH delta exp async BR: ',num2str(SH_delta_async),'\n \n'])
    
    fprintf(['AH sync BR: ',num2str(AH_sync),'\n'])
    fprintf(['AH exp async BR: ',num2str(AH_exp_async),'\n'])
    fprintf(['AH delta exp async BR: ',num2str(AH_delta_async),'\n'])
    fprintf(['AH group async BR: ',num2str(AH_group_async),'\n \n'])
    
        %% Verify relative BRs
    if (AH_sync - buffer > SH_sync) || (SH_sync - buffer > B_sync)
        error(['Check sync encodings, something wrong'])
    end
    if (AH_exp_async - buffer > SH_exp_async) || (SH_exp_async - buffer > B_exp_async)
        error(['Check exp async encodings, something wrong'])
    end
    if (AH_delta_async - buffer > SH_delta_async)
        error(['Check delta exp async encodings, something wrong'])
    end
    if (AH_group_async - buffer > B_group_async)
        fprintf(['Check group async encodings, something wrong \n'])
    end

%     save(save_file,'all_results','SH_delta_async')

    done = toc(timer);
    fprintf(['Whole script took ',num2str(done),' s\n'])
