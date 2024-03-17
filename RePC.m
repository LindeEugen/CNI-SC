
function [PC,sepset] = RePc(Data,target,lamda,ADJT)     
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       lamda is the significance level
% OUTPUT:
%       PC is the parents and children of the target
%       time is the runtime of the algorithm
%       sepset: the condition sets for each node to make it independt to the target

[samples,p]=size(Data);
maxK=3;
start=tic;
test=0;
sepset = cell(1,p);
ADJT_length = length(ADJT);
cutSetSize = 0;
%/* Initialise cutset size to 0 */
tmp_ADJT = ADJT;
while ADJT_length > cutSetSize&&cutSetSize<=maxK
    %/* Generate all sets of size cutSetSize as candidate cut sets */
    for i=1:length(ADJT)
        X = ADJT(i);
        nbrs = setdiff(tmp_ADJT, X);
        SS = subsets(nbrs, cutSetSize);    
        %/* Conduct conditional independence tests to remove non-causal nodes */
        for si=1:length(SS)
            S = SS{si};
            test=test+1;
            CI=SCI(Data(:,X),Data(:,target),Data(:,S));     
            if(CI==0)            
                tmp_ADJT = setdiff(tmp_ADJT,X);
                ADJT_length = ADJT_length-1;
                sepset{1,X} = S;
                break;          
            end
        end
        
    end

    ADJT=tmp_ADJT;
    
    cutSetSize = cutSetSize + 1;
end
PC=ADJT;
time=toc(start);


