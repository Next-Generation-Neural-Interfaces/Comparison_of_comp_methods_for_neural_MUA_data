    
function [check_B_sync,B_sync] = check_sync_basic(BP_MUA,basic_sync_codeword_len)

    nb_channels = length(BP_MUA(1,:));
    data_length = length(BP_MUA(:,1)) * nb_channels;
    
    % Used for comparing to the decoded data
    r_BP_MUA = reshape(BP_MUA',length(BP_MUA(:,1))*nb_channels,1); % reshaped BP_MUA as a 1d vector
            
    B_sync_enc = [];
    for t = 1:length(BP_MUA(:,1))
        for j = 1:nb_channels
            B_sync_enc = [B_sync_enc, dec2bin(BP_MUA(t,j),basic_sync_codeword_len)];
        end
    end
    B_sync = length(B_sync_enc)/data_length; % average symbol length in encoding

    % Decode
    D = zeros(data_length,1);
    k = 0;
    for counter = 1:basic_sync_codeword_len:length(B_sync_enc)
        k = k + 1;
        D(k) = bin2dec(B_sync_enc(counter:counter+basic_sync_codeword_len-1));
    end
    check_B_sync = sum(abs(r_BP_MUA-D)) == 0;
    
    if check_B_sync == 0
        error('Check the Basic Sync encoding, the decoded data does not match the original \n')
    end

end