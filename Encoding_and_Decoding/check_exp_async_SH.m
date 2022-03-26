function [check_SH_exp_async,SH_exp_async] = check_exp_async_SH(BP_MUA, mapping, inv_mapping, basic_exp_async_chan_ID_codeword_len, S)

    % Explicit async BR SH
    
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
    
    % Train SH decoder on a decaying exponential, for the 0 =< firing
    % rates =< S-1.
    decaying_prob_vector = 1./exp(4*(1:S-1)); % This should guarantee codeword lengths of 1,2,...,S-2,S-2.
    decaying_prob_vector = decaying_prob_vector ./ sum(decaying_prob_vector);
    [dict_exp_async_FR_SH,~] = huffmandict(1:S-1,decaying_prob_vector);
    
    % Initialise
    SH_exp_async_enc = []; % encoded data
    SH_exp_async = 0; % data length
    
    % Iterate through timesteps, find firing rates above 0
    for t = 1:length(BP_MUA_copy(:,1))
        for j = 1:nb_channels
            FR = BP_MUA_copy(t,j);
            if FR > 0 % different to most common FR for that channel
                % Encode channel ID. fixed length codewor. 0 = chan 1, 1 =
                % chan 2, etc.
                % Maybe use a same-length Huffman encoder?
                temp_ID = dec2bin(j-1,basic_exp_async_chan_ID_codeword_len);
                
                % Encode firing rate, SH codeword
                temp_FR_codeword = dict_exp_async_FR_SH(FR,2); 

                SH_exp_async_enc = [strcat(SH_exp_async_enc, temp_ID), num2str(temp_FR_codeword{1})];
                SH_exp_async  = SH_exp_async + length(temp_ID) + length(temp_FR_codeword{1});
            end
        end
        
        % Add stop symbol between timesteps so we can tell them apart. THis
        % doesn't contribute to the encoded data length
        SH_exp_async_enc = strcat(SH_exp_async_enc, 'stop');
    end
    % Remove blank spaces added by FR codeword
    SH_exp_async_enc = strrep(SH_exp_async_enc,' ',"");
    SH_exp_async_enc = SH_exp_async_enc{1}; % formatting
    
    % Average bits per encoded symbol
    SH_exp_async = SH_exp_async/(nb_channels * length(BP_MUA_copy(:,1)));

    
    % Decode
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
    while counter < length(SH_exp_async_enc)
        
        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(SH_exp_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            
            % Last stop symbol, we end the while loop
            if counter > length(SH_exp_async_enc)
                break
            end
            continue % important so that the algorithm doesn't fail if two stop symbols in a row
        end
        
        % First, we look at the channel ID.
        temp_ID = bin2dec(SH_exp_async_enc(counter:counter+basic_exp_async_chan_ID_codeword_len-1)) + 1;
        counter = counter + basic_exp_async_chan_ID_codeword_len;
        
        % Then at firing rates
        flag_1 = 0; % flag == 1 when we find a correct Huffman codeword for a FR.
        increment = 0; % how many bits we're currently looking at for the codeword
        while flag_1 == 0
           
            temp_fr_enc = SH_exp_async_enc(counter:counter+increment); % our sample, we check if it's eqaul to a SH codeword for FRs.

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
                decoded_FR = huffmandeco(numeric_temp_fr_enc,dict_exp_async_FR_SH);
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
            D(row,temp_ID) = decoded_FR;
        else
            D(row,temp_ID) = inv_mapping(temp_ID,decoded_FR+1)-1; % inverse mapping
        end
        
    end
    check_SH_exp_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
    if check_SH_exp_async == 0
        error('Check the SH Explicit Async encoding, the decoded data does not match the original')
    end
end