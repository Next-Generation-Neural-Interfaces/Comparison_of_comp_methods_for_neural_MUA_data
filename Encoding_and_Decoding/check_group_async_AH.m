     
function [check_AH_group_async,AH_group_async] = check_group_async_AH(BP_MUA,mapping,inv_mapping,S,rel_combined_prob_group)

    nb_channels = length(BP_MUA(1,:));
    
    % NOTE: the exp- group Huffman dicts should be trained on the relative
    % prob of channel IDs AND stop symbols, since the pairing is
    % unknown: unlike with explicit async, you don't know that the
    % nb of events will follow the channel ID. In this case, you
    % may have x channel IDs, then y stop symbols, with no idea what the
    % sequence is. Therefore, you need to have the codewoirds be deodable
    % by the same encoder, and do the probs need to be relative to eachother, so
    % it's 1 encoder.            

    % Get probs, train AH on it. Combined encoder for channel IDs and stop
    % symbols.
    [dict_group_async_AH,~] = huffmandict(1:S+nb_channels-2,rel_combined_prob_group);

    % Map BP_MUA directly
    BP_MUA_copy = BP_MUA;
    if mapping ~= 0
        for t = 1:length(BP_MUA(:,1))
            for j = 1:nb_channels
                FR = BP_MUA_copy(t,j);
                BP_MUA_copy(t,j) = mapping(j,FR+1)-1; % + and - 1 for MATLAB indexing
            end
        end
    end
    
    % Encode data
    AH_group_async_enc = [];
    AH_group_async = 0;
    for t = 1:length(BP_MUA_copy(:,1))
        collated_encoding = [];
        for i = 1:S-1
            events = find(BP_MUA_copy(t,:)==i); % channel IDs that have firing rate == i (or if mapping, channels that have their i^{th} most common fr)
            if i == 1 % no stop symbol
                FR_stop_enc = [];
            else
                FR_stop_enc = huffmanenco(nb_channels+i-1,dict_group_async_AH)';
            end
            % Really weird, MATLAB code transposes the output depending on
            % the number of inputs to the encoder
            if length(events) == 1
                channel_IDs_enc = huffmanenco(events,dict_group_async_AH)';
            else
                channel_IDs_enc = huffmanenco(events,dict_group_async_AH);
            end
            
            collated_encoding = [collated_encoding, FR_stop_enc, channel_IDs_enc];
            
            AH_group_async = AH_group_async + length(channel_IDs_enc) + length(FR_stop_enc);
        end

        collated_encoding = num2str(collated_encoding);
        AH_group_async_enc = [AH_group_async_enc, collated_encoding, 'stop'];
    end
    % Remove blank spaces added by FR codeword
    AH_group_async_enc = strrep(AH_group_async_enc,' ',"");
    AH_group_async_enc = AH_group_async_enc{1}; % formatting
    
    % Average bits per encoded symbol
    AH_group_async = AH_group_async/(nb_channels * length(BP_MUA(:,1)));

    % Decode
    D = int16(zeros(length(BP_MUA(:,1)),nb_channels));
    
    % IMPORTANT: We prepare the decoded matrix (D): we give each channel, for each
    % timestep, its most common FR since this is the assumed default,
    % overwritten by any communicated data for that timestep.
    if mapping ~= 0
        for j = 1:nb_channels
            D(:,j) = find(mapping(j,:)==1)-1; % '-1' because of the MATLAB indexing
        end
    end
    
    counter = 1; row = 1; firing_rate = 1;
    while counter < length(AH_group_async_enc)

        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(AH_group_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            firing_rate = 1; % we reset the firing rate to 1

            % Last stop symbol, we end the while loop
            if counter > length(AH_group_async_enc)
                break
            end
            continue % important so that the algorithm doesn't fail if two stop symbols in a row
        end

        %% Decode channel IDs and stop symbols
        % We progressively increase the number of bits we're looking at
        % untiul we fing a legitimate codeword from out encoder.
        flag_1 = 0; % flag == 1 when we find a correct Huffman codeword
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_1 == 0
           
            temp_enc = AH_group_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for FRs.

            % Convert to numeric array, because Huffman decoding formatting is a pain
            numeric_temp_enc = zeros(1,length(temp_enc));
            for k = 1:length(temp_enc)
                numeric_temp_enc(k) = str2num(temp_enc(k));
            end
            
            try
                % If a real codeword, we found the correct one (a benefit of
                % prefix coding: there is only one valid codeword and it
                % can be found by increasing the amount of considered bits
                % and checking if the codeword matches one of the Huffman
                % codewords)
                decoded_temp = huffmandeco(numeric_temp_enc,dict_group_async_AH);
                if ~isempty(decoded_temp)
                    flag_1 = 1; % we're done, we found the codeword
                end
            catch
            end
            increment = increment + 1;
            
            if increment > S + nb_channels - 2
                error('We have a problem, we have not found the correct codeword in the SH explicit async encoded data during decoding')
            end
        end
        counter = counter + increment; % increase counter by however many bits were in the firing rate codeword
        

        % It's a stop symbol (group async firing rate stop symbol, not between timesteps)
        if decoded_temp > nb_channels
            firing_rate = decoded_temp-nb_channels+1;
        else % if it's a channel ID
            
            if mapping == 0
                D(row,decoded_temp) = firing_rate;
            else
                D(row,decoded_temp) = inv_mapping(decoded_temp,firing_rate+1)-1;
            end
        end

    end
    check_AH_group_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
    if check_AH_group_async == 0
        error('Check the AH Group Async encoding, the decoded data does not match the original')
    end
end
