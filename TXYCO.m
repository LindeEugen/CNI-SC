function NEW_CAM=TXYCO(T,X,Y,CAM)
    if CAM==0
        inter_dir_sp=[];
        inter_sp=[];
    end
    [index]=find(Y(:,1)==0);
    Y(index,:)=[];
    for j=1:size(Y,1)  
        vdir(T,Y(j,1))=1;     
        %v structural orientation
        %stru(sp(j,1),i)=0;
        inter_dir_sp=union(inter_dir_sp,Y(j,:));   
        %y oriented parent-child nodes and their spouses
        inter_sp=union(inter_sp,Y(j,2:end)); 
        %Spouse Node
    end
    inter_sp(inter_sp==0)=[];
    inter_dir_sp(inter_dir_sp==0)=[];
    ture_pc=mysetdiff(MB,inter_sp);   
    %Real parent-child nodes
    stru(T,ture_pc)=1;
    %% Processing the structure matrix and outputting pairs of undetermined relations
    [index1,index2]=find(vdir==1);
    stru(index2,index1)=0;       
    %Setting the other direction with the v-structure
    nodirstru=stru-vdir; 
    %Remove undirected structures
    nodirstru(nodirstru<0)=0; 
    %Deletion of negative elements (i.e., directed structures)
    flag=1;
    for i=1:n
        for j=i:n
            if nodirstru(i,j)==1||nodirstru(j,i)==1
                dir_index(flag,1:2)=[i,j];
                flag=flag+1;
            end
        end
    end
end 