function [CPCT,sepset]=RePC(T,ADJT,nodes,data)
    maxK=3;
    sepset = cell(1,nodes);
    cutSetSize = 0;
    tmp_ADJT=ADJT(find(ADJT~=T));
    ADJT_length = length(tmp_ADJT);
    temp=tmp_ADJT;
% Step1:initialize
    while ADJT_length > cutSetSize&&cutSetSize<=maxK
        for i=1:length(tmp_ADJT)
            X = temp(i);
            nbrs = setdiff(tmp_ADJT, X);
            SS = subsets(nbrs, cutSetSize);
% Step2: Perform the GenADJ algorithm to generate ADJ(T), and generate a series of candidate cut sets S of size cutSetSize for each node X in ADJ(T)            
            for si=1:length(SS)
                S = SS{si};
                CI=SCI(data(:,X),data(:,T),data(:,S));     
                if(CI<=0)            
                    tmp_ADJT = setdiff(tmp_ADJT,X);
                    ADJT_length = ADJT_length-1;
                    sepset{1,X} = S;
                    break;
% Step3: Treat each candidate cut set S as a conditional set separately, and test each node X in ADJ(T) for conditional independence from the target node T                    
                end
            end
        end
    ADJT=tmp_ADJT;
    cutSetSize = cutSetSize + 1;
    end
    CPCT=ADJT;
% Step 4: Remove NonPC from ADJ(T), the set after deletion the paternal child node set of the target node T, noted asPC(T).
end