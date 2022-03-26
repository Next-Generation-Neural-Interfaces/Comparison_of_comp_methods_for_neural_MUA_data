function [check_B_group_async,B_group_async] = check_group_async_basic(BP_MUA,basic_group_asycn_codeword_len,S)

    % Our strategy is to do all the collation in decimal, and then convert the
    % entire sequence into binary.
    
    nb_channels = length(BP_MUA(1,:));

    % Get stop symbols
    stop = {[]}; % first stop symbol (between 0 and 1) is empty
    for i = 1:S-2
        stop{i+1} = nb_channels+i-1;
    end

    % Encode
    B_group_async = 0; % data length
    B_group_async_enc = ''; % encoded data
    for t = 1:length(BP_MUA(:,1)) % iterate through timesteps
        collated_events  = []; % per timestep
        for i = 1:S-1
            events = find(BP_MUA(t,:)==i);
            collated_events = [collated_events, stop{i}, events-1];
        end
        temp = convertCharsToStrings(dec2bin(collated_events,basic_group_asycn_codeword_len)'); % collated char array into 1 binary string
        B_group_async = B_group_async + length(temp{1}); % how many bits

        % We add a stop symbol between timesteps so we can tell them apart.
        % This doesn't contribute to the encoded data length
        B_group_async_enc = strcat(B_group_async_enc,temp{1},'stop'); 
    end
    B_group_async = B_group_async/(length(BP_MUA(:,1))*nb_channels); % average number of bits to encode each symbol in the original data
    % B_group_async_enc = B_group_async_enc{1};

    % Decode
    D = int16(zeros(length(BP_MUA(:,1)),nb_channels));
    counter = 1; row = 1; firing_rate = 1;
    while counter < length(B_group_async_enc)

        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(B_group_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            firing_rate = 1; % we reset the firing rate to 1

            % Last stop symbol, we end the while loop
            if counter > length(B_group_async_enc)
                break
            end
        end

        % We're decoding channel IDs and stop symbols
        temp = bin2dec(B_group_async_enc(counter:counter+basic_group_asycn_codeword_len-1));
        % Since we're not encoding a fr of 0, this makes everything a bit more
        % intuitive. e.g. 0000 -> 1 = chan. 1, since there is no channel 0.
        temp = temp + 1; 
        counter = counter + basic_group_asycn_codeword_len;

        % It's a stop symbol (group async firing rate stop symbol, not between timesteps)
        if temp > nb_channels
            firing_rate = temp-nb_channels+1;
        else
            D(row,temp) = firing_rate;
        end

    end
    check_B_group_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
    if check_B_group_async == 0
        error('Check the Basic Group Async encoding, the decoded data does not match the original')
    end
end

