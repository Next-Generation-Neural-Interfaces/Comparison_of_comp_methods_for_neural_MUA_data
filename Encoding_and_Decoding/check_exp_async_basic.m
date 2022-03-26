function [check_B_exp_async,B_exp_async] = check_exp_async_basic(   BP_MUA, basic_exp_async_chan_ID_codeword_len, basic_exp_async_nb_MUA_events_codeword_len)

    % Explicit async BR basic
    
    nb_channels = length(BP_MUA(1,:));
    
    % Initialise
    B_exp_async_enc = []; % encoded data
    B_exp_async = 0; % data length
    
    % Iterate through timesteps, find firing rates above 0
    for t = 1:length(BP_MUA(:,1))
        for j = 1:nb_channels
            if BP_MUA(t,j) > 0
                % Encode channel ID
                temp = dec2bin(j,basic_exp_async_chan_ID_codeword_len);
                % Encode firing rate
                if basic_exp_async_nb_MUA_events_codeword_len > 0 % if S > 2
                    temp = [temp,dec2bin(BP_MUA(t,j),basic_exp_async_nb_MUA_events_codeword_len)];
                end
                B_exp_async_enc = [B_exp_async_enc, temp];
                B_exp_async  = B_exp_async + length(temp);
            end
        end
        
        % Add stop symbol between timesteps so we can tell them apart. THis
        % doesn't contribute to the encoded data length
        B_exp_async_enc = [B_exp_async_enc, 'stop'];
    end
    
    % Average bits per encoded symbol
    B_exp_async = B_exp_async /(nb_channels * length(BP_MUA(:,1)));
    
    % Decode
    D = int16(zeros(length(BP_MUA(:,1)),nb_channels));
    counter = 1; row = 1;
    while counter < length(B_exp_async_enc)
        
        % Stop symbol, means we're at the next timestep. In practice, this
        % won't be included in the communicated data since we assume the
        % communication rate is sufficient to have a buffer between
        % timesteps.
        if strcmp(B_exp_async_enc(counter:counter+3),'stop')
            row = row + 1;
            counter = counter + 4;
            
            % Last stop symbol, we end the while loop
            if counter > length(B_exp_async_enc)
                break
            end
        end
        
        % We're looking at a channel ID.
        temp_ID = bin2dec(B_exp_async_enc(counter:counter+basic_exp_async_chan_ID_codeword_len-1));
        counter = counter + basic_exp_async_chan_ID_codeword_len;
        
        % We're looking at firing rates
        if basic_exp_async_nb_MUA_events_codeword_len > 0
            temp_fr = bin2dec(B_exp_async_enc(counter:counter+basic_exp_async_nb_MUA_events_codeword_len-1));
            counter = counter + basic_exp_async_nb_MUA_events_codeword_len;
        else
            temp_fr = 1;
        end

        D(row,temp_ID) = temp_fr;
    end
    check_B_exp_async = sum(sum(abs(int16(BP_MUA)-D))) == 0;
    
    if check_B_exp_async == 0
        error('Check the Basic Explicit Async encoding, the decoded data does not match the original')
    end
end