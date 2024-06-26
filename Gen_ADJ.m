function ADJT=Gen_ADJ(lamda,data,nodes,T) 
    U=1:nodes;
    CADJ=U(find(U~=T));
%     Step1: Remove the target node T from all node sets U
    data_temp_i=data;
    data_temp_j=data;
%     Setting up temporary variables
    for X=1:length(CADJ)
        data(:,1)=data_temp_i(:,CADJ(X));
        data(:,2)=data_temp_j(:,T);
        data(:,CADJ(X))=data_temp_i(:,1);
        data(:,T)=data_temp_j(:,2);
 %     Use temporary variable assignments to avoid logical errors
        node_score(X)=SCI(data(:,1),data(:,2),0);
        if node_score(X)<lamda
            ADJ(X)=CADJ(X);
        else node_score(X)=[];
        end
%     Step 2:Judge the strength of the dependency between each node X in CADJ(T) and the target node T
    end
    [node_score, idx] = sort(node_score);
    ADJT = ADJ(idx);
%     Step 3: Ascend order of the variables ADJ(T) by SCI score, result in an adjacency node set of final target nodes T.
end