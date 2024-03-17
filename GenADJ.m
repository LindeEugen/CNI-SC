function ADJT=GenADJ(T,Data,lamda,SCIM) 
    [CADJ,SepT]=PCsimple_Z(data,T,lamda,SCIM);   
    for i=1:length(CADJ)
        X = CADJ(i);
        ADJ = PCsimple_Z(data,X,lamda,CADJ);
        if ~ismember(T,ADJ)
            continue;
        end
        for j=1:length(ADJ)
            Y=PCX(j);
            if Y==T
                continue;
            end
            ADJ=union(SepT{1,Y},X);
        end
    end
end