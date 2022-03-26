     
function [check_AH_sync,AH_sync] = check_sync_AH(BP_MUA,mapping,inv_mapping,S,p_channel_mean)

    % This is highly trivial: it's just training a Huffman encoder,
    % encoding the data, and then decoding, all using the standard Huffman
    % functions.
    
    nb_channels = length(BP_MUA(1,:));
    data_length = length(BP_MUA(:,1)) * nb_channels;
    
    % Map data
    r_BP_MUA = BP_MUA;
    if mapping ~= 0
        for t = 1:length(BP_MUA(:,1))
            for j = 1:nb_channels
                r_BP_MUA(t,j) = mapping(j,r_BP_MUA(t,j)+1) -1; % + and - 1 for MATLAB indexing
            end
        end
    end
    
    % Used for comparing to the decoded data
    r_BP_MUA = reshape(r_BP_MUA',length(BP_MUA(:,1))*nb_channels,1); % reshaped BP_MUA as a 1d vector
    
    % Get probs, train AH on it
    [dict_sync_AH,avglen_sync] = huffmandict(0:S-1, p_channel_mean);

    AH_sync_enc = huffmanenco(r_BP_MUA,dict_sync_AH);
    AH_sync = length(AH_sync_enc)/data_length; % average symbol length in encoding
    % Rounding errors can make the below throw up an error sometimes, not
    % trustworthy as a test
%     if AH_sync ~= avglen_sync
%         error('Something wrong, encoded lengths do not match for AH Sync encoding. "avglen_sync" given by the Huffmman dict function should be equal to the average symbol length derived analytically. See code for details.')
%     end
    
    
    % Decode
    D = huffmandeco(AH_sync_enc,dict_sync_AH);
   
    % Inverse mapping
    D = reshape(D,nb_channels,length(BP_MUA(:,1)))'; % we reshape to 2D matrix
    if mapping ~= 0
        for t = 1:length(BP_MUA(:,1))
            for j = 1:nb_channels
                D(t,j) = inv_mapping(j,D(t,j)+1)-1; % inverse mapping
            end
        end
    end

    check_AH_sync = sum(sum(abs(BP_MUA-D))) == 0;
    
    if check_AH_sync == 0
        error('Check the AH Sync encoding, the decoded data does not match the original')
    end

end