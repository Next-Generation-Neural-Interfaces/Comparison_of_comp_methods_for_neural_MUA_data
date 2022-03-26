function [check_AH_delta_async,AH_delta_async] = check_delta_async_AH(BP_MUA, mapping, inv_mapping, S)

    nb_channels = length(BP_MUA(1,:));
    
    % Map BP_MUA directly
    BP_MUA_copy = BP_MUA; % for storage
    if mapping ~= 0
        for t = 1:length(BP_MUA_copy(:,1))
            for j = 1:nb_channels
                FR = BP_MUA_copy(t,j);
                BP_MUA_copy(t,j) = mapping(j,FR+1)-1; % + and - 1 for MATLAB indexing
            end
        end
    end

    % Get delta-sample channel ID frequencies to train the encoder. While
    % we're at, we get the firing rates too (with mapping option). It's a
    % bit slow, since we iterate through the data to get the frequencies,
    % then again to encode it.

    % These arrays are larger than necessary, but saves time in the for
    % loop. The values == 0 are ignored anyway in the histogram below.
    stored_d_channel_IDs = zeros(nb_channels*length(BP_MUA_copy(:,1)),1);
    stored_firing_rates = zeros(nb_channels*length(BP_MUA_copy(:,1)),1);
    stored_d_channel_IDs_counter = 1; % faster to index an existing array than increase its size +1 for each added eMUA event above 0
    stored_firing_rates_counter = 1;
    for timestep = 1:length(BP_MUA_copy(:,1))
        samp = BP_MUA_copy(timestep,:); % our multi-channel timestep of the MUA data
        a1 = find(samp>0); % IDs of channels with firing rates above 0
        for event = 1:length(a1) % iterate through each MUA event in timestep
            if event == 1 % If it's the first event in the timestep, we assume it's delta-sampled relative to 0.
                d_chan_ID = a1(event);
            else % If it's the 2nd or more, we delta-sample the channel IDs
                d_chan_ID = a1(event) - a1(event-1); % delta sampled channel ID
            end
            % Firing rate.
            fr = samp(a1(event));

            % Store all channel ID lengths in this variable
            stored_d_channel_IDs(stored_d_channel_IDs_counter) = d_chan_ID; % add channel ID length
            stored_d_channel_IDs_counter = stored_d_channel_IDs_counter + 1;
            stored_firing_rates(stored_firing_rates_counter) = fr;
            stored_firing_rates_counter = stored_firing_rates_counter + 1;
        end
    end

    rel_prob_delta_channel_IDs = histogram(stored_d_channel_IDs,0.5:1:nb_channels+0.5);
    rel_prob_delta_channel_IDs = rel_prob_delta_channel_IDs.Values./sum(rel_prob_delta_channel_IDs.Values);

    rel_firing_rates_above_1 = histogram(stored_firing_rates,0.5:1:S-0.5);
    rel_firing_rates_above_1 = rel_firing_rates_above_1.Values./sum(rel_firing_rates_above_1.Values);
    
    % Channel IDs encoder
    [dict_AH_delta_channel_ID,~] = huffmandict(1:nb_channels,rel_prob_delta_channel_IDs);
    % Firing rate encoder
    [dict_AH_MUA_events_above_1,~] = huffmandict(1:S-1,rel_firing_rates_above_1);
    

    % Encode data
    AH_delta_async_enc = [];
    AH_delta_async = 0;
    for t = 1:length(BP_MUA_copy(:,1))
        events = find(BP_MUA_copy(t,:));
        collated_encoding = [];
        previous_channel = 0; % reset delta-sampling
        
        for i = 1:length(events)
            % Channel ID
            delta_chan_ID = events(i) - previous_channel;
            previous_channel = events(i);
            delta_chan_ID_enc = huffmanenco(delta_chan_ID,dict_AH_delta_channel_ID);
            
            % FR
            FR = BP_MUA_copy(t,events(i));
            FR_enc = huffmanenco(FR,dict_AH_MUA_events_above_1);
            
            collated_encoding = [collated_encoding, delta_chan_ID_enc', FR_enc'];
            
            AH_delta_async = AH_delta_async + length(delta_chan_ID_enc) + length(FR_enc);
        end
        collated_encoding = num2str(collated_encoding);
        AH_delta_async_enc = [AH_delta_async_enc, collated_encoding, 'stop'];
    end
    
    % Remove blank spaces added by FR codeword
    AH_delta_async_enc = strrep(AH_delta_async_enc,' ',"");
    AH_delta_async_enc = AH_delta_async_enc{1}; % formatting
    
    % Average bits per encoded symbol
    AH_delta_async = AH_delta_async/(nb_channels * length(BP_MUA_copy(:,1)));

    % Decode data
    D = int16(zeros(length(BP_MUA_copy(:,1)),nb_channels));
    % IMPORTANT: We prepare the decoded matrix (D): we give each channel, for each
    % timestep, its most common FR since this is the assumed default,
    % overwritten by any communicated data for that timestep.
    if mapping ~= 0
        for j = 1:nb_channels
            D(:,j) = find(mapping(j,:)==1)-1; % '-1' because of the MATLAB indexing
        end
    end
    
    counter = 1; row = 1;
    previous_channel = 0;
    while counter < length(AH_delta_async_enc)
        
        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(AH_delta_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            
            % Last stop symbol, we end the while loop
            if counter > length(AH_delta_async_enc)
                break
            end
            previous_channel = 0; % reset delta
            continue % important so that the algorithm doesn't fail if two stop symbols in a row
        end
        
        %% First, we look at the channel ID.
        flag_delta_chan_ID = 0; % flag == 1 when we find a correct Huffman codeword for a chan ID.
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_delta_chan_ID == 0
            temp_delta_chan_ID_enc = AH_delta_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for chan IDs.

            % Convert string to numeric array, because Huffman decoding formatting is a pain
            numeric_temp_delta_chan_ID_enc = zeros(1,length(temp_delta_chan_ID_enc));
            for k = 1:length(temp_delta_chan_ID_enc)
                numeric_temp_delta_chan_ID_enc(k) = str2num(temp_delta_chan_ID_enc(k));
            end

            try
                % If a real codeword, we found the correct one (a benefit of
                % prefix coding: there is only one valid codeword and it
                % can be found by increasing the amount of considered bits
                % and checking if the codeword matches one of the Huffman
                % codewords)
                decoded_delta_chan_ID = huffmandeco(numeric_temp_delta_chan_ID_enc,dict_AH_delta_channel_ID);
                if ~isempty(decoded_delta_chan_ID)
                    flag_delta_chan_ID = 1; % we're done, we found the codeword
                end
            catch
            end
            increment = increment + 1;

            if increment > nb_channels
                error('We have a problem, we have not found the correct codeword in the AH explicit async encoded data during decoding')
            end
        end
        channel_ID = previous_channel + decoded_delta_chan_ID;
        previous_channel = channel_ID;
        counter = counter + increment; % increase counter by however many bits were in the firing rate codeword
        
        %% Then we decode the following firing rate
        flag_1 = 0; % flag == 1 when we find a correct Huffman codeword for a FR.
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_1 == 0
           
            temp_fr_enc = AH_delta_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for FRs.

            % Convert to numeric array, because Huffman decoding formatting is a pain
            numeric_temp_fr_enc = zeros(1,length(temp_fr_enc));
            for k = 1:length(temp_fr_enc)
                numeric_temp_fr_enc(k) = str2num(temp_fr_enc(k));
            end
            
            try
                % If a real codeword, we found the correct one (a benefit of
                % prefix coding: there is only one valid codeword and it
                % can be found by increasing the amount of considered bits
                % and checking if the codeword matches one of the Huffman
                % codewords)
                decoded_FR = huffmandeco(numeric_temp_fr_enc,dict_AH_MUA_events_above_1);
                if ~isempty(decoded_FR)
                    flag_1 = 1; % we're done, we found the codeword
                end
            catch
            end
            increment = increment + 1;
            
            if increment > S
                error('We have a problem, we have not found the correct codeword in the SH explicit async encoded data during decoding')
            end
        end
        counter = counter + increment; % increase counter by however many bits were in the firing rate codeword
        
        % Store decoded data
        if mapping == 0
            D(row,channel_ID) = decoded_FR;
        else
            D(row,channel_ID) = inv_mapping(channel_ID,decoded_FR+1)-1; % inverse mapping
        end
        
        
    end
    check_AH_delta_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
    if check_AH_delta_async == 0
        error('Check the AH Delta Async encoding, the decoded data does not match the original')
    end
end
