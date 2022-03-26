function [check_AH_exp_async,AH_exp_async] = check_exp_async_AH(BP_MUA, mapping,inv_mapping, S, rel_p_send_chan, rel_p_positive_event_mean_all_channels)

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

    % Channel IDs encoder
    [dict_AH_channel_ID,~] = huffmandict(1:nb_channels,rel_p_send_chan);
    
    % Firing rate encoder
    [dict_AH_MUA_events_above_1,~] = huffmandict(1:S-1,rel_p_positive_event_mean_all_channels(2:end));
    
    % Encode data
    AH_exp_async_enc = [];
    AH_exp_async = 0;
    for t = 1:length(BP_MUA_copy(:,1))
        events = find(BP_MUA_copy(t,:)); % channels with firing rates above 0 (or if mapped, different to most common FR in that channel)
        collated_encoding = [];
        
        for i = 1:length(events)
            channel = events(i);
            chan_ID_enc = huffmanenco(channel,dict_AH_channel_ID);
                           
            FR = BP_MUA_copy(t,channel);
                
            FR_enc = huffmanenco(FR,dict_AH_MUA_events_above_1);
            collated_encoding = [collated_encoding, chan_ID_enc', FR_enc'];
            
            AH_exp_async = AH_exp_async + length(chan_ID_enc) + length(FR_enc);
        end
        collated_encoding = num2str(collated_encoding);
        AH_exp_async_enc = [AH_exp_async_enc, collated_encoding, 'stop'];
    end
    % Remove blank spaces added by FR codeword
    AH_exp_async_enc = strrep(AH_exp_async_enc,' ',"");
    AH_exp_async_enc = AH_exp_async_enc{1}; % formatting
    
    % Average bits per encoded symbol
    AH_exp_async = AH_exp_async/(nb_channels * length(BP_MUA_copy(:,1)));

    
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
    while counter < length(AH_exp_async_enc)
        
        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(AH_exp_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            
            % Last stop symbol, we end the while loop
            if counter > length(AH_exp_async_enc)
                break
            end
            continue % important so that the algorithm doesn't fail if two stop symbols in a row
        end
        
        %% First, we look at the channel ID.
        flag_chan_ID = 0; % flag == 1 when we find a correct Huffman codeword for a chan ID.
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_chan_ID == 0
            temp_chan_ID_enc = AH_exp_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for chan IDs.

            % Convert string to numeric array, because Huffman decoding formatting is a pain
            numeric_temp_chan_ID_enc = zeros(1,length(temp_chan_ID_enc));
            for k = 1:length(temp_chan_ID_enc)
                numeric_temp_chan_ID_enc(k) = str2num(temp_chan_ID_enc(k));
            end

            try
                % If a real codeword, we found the correct one (a benefit of
                % prefix coding: there is only one valid codeword and it
                % can be found by increasing the amount of considered bits
                % and checking if the codeword matches one of the Huffman
                % codewords)
                decoded_chan_ID = huffmandeco(numeric_temp_chan_ID_enc,dict_AH_channel_ID);
                if ~isempty(decoded_chan_ID)
                    flag_chan_ID = 1; % we're done, we found the codeword
                end
            catch
            end
            increment = increment + 1;

            if increment > nb_channels
                error('We have a problem, we have not found the correct codeword in the AH explicit async encoded data during decoding')
            end
        end
        counter = counter + increment; % increase counter by however many bits were in the firing rate codeword

        %% Then we decode the following firing rate
        flag_1 = 0; % flag == 1 when we find a correct Huffman codeword for a FR.
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_1 == 0
           
            temp_fr_enc = AH_exp_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for FRs.

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
            D(row,decoded_chan_ID) = decoded_FR;
        else
            D(row,decoded_chan_ID) = inv_mapping(decoded_chan_ID,decoded_FR+1)-1; % inverse mapping
        end
        
    end
    check_AH_exp_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
%     if check_AH_exp_async == 0
%         error('Check the AH Explicit Async encoding, the decoded data does not match the original')
%     end
end
    