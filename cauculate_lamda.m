function lamda=cauculate_lamda(data,nodes)
    data_temp_i=data;
    data_temp_j=data;
%     Setting up temporary variables
    for i=1:nodes
        for j=1:nodes
            data(:,1)=data_temp_i(:,i);
            data(:,2)=data_temp_j(:,j);
            data(:,i)=data_temp_i(:,1);
            data(:,j)=data_temp_j(:,2);
%             Use temporary variable assignments to avoid logical errors
            if i~=j
                SCIM((i-1)*nodes+j)=SCI(data(:,1),data(:,2),0);
            else 
                SCIM((i-1)*nodes+j)=0;
            end
        end
    end
    lamda=min(SCIM)+var(SCIM)*(mean(SCIM)-min(SCIM));
end