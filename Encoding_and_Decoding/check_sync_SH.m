     
function [check_SH_sync,SH_sync] = check_sync_SH(BP_MUA,mapping,inv_mapping,S)

    nb_channels = length(BP_MUA(1,:));
    data_length = length(BP_MUA(:,1)) * nb_channels;
    
    % Train SH decoder on a decaying exponential
    decaying_prob_vector = 1./exp(4*(0:S-1)); % This should guarantee codeword lengths of 1,2,3,...,S-1,S-1.
    decaying_prob_vector = decaying_prob_vector ./ sum(decaying_prob_vector);
    [dict_sync_SH,~] = huffmandict(0:S-1,decaying_prob_vector);
    
    SH_sync_enc = [];
    for t = 1:length(BP_MUA(:,1))
        for j = 1:nb_channels
            FR = BP_MUA(t,j);
            if mapping ~= 0
                FR = mapping(j,FR+1);
            else
                FR = FR + 1; % for indexing Huffman (MATLAB indexing starts at 1)
            end
            temp = dict_sync_SH(FR,2);
            SH_sync_enc = [SH_sync_enc, temp{1}];
        end
    end
    SH_sync = length(SH_sync_enc)/data_length; % average symbol length in encoding

    % Decode
    D = huffmandeco(SH_sync_enc,dict_sync_SH);
    
    % Inverse mapping
    D = reshape(D,nb_channels,length(BP_MUA(:,1)))'; % we reshape to 2D matrix
    if mapping ~= 0
        for t = 1:length(BP_MUA(:,1))
            for j = 1:nb_channels
                D(t,j) = inv_mapping(j,D(t,j)+1)-1; % inverse mapping
            end
        end
    end
    
    check_SH_sync = sum(sum(abs(BP_MUA-D))) == 0;
    
    if check_SH_sync == 0
        error('Check the SH Sync encoding, the decoded data does not match the original')
    end

end