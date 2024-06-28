data = importdata('alarm_500.txt');
[samples,nodes]=size(data)
Vstru=[];
nodir=[];
CanT=1:nodes;
length_CanT=length(CanT)
lamda=cauculate_lamda(data,nodes)
%  Initialize
while length_CanT>0
    for t=1:length_CanT
%   structure learning and direction learning are performed on the target node T to construct a local causal network about T
        T=CanT(t);
        ADJT=Gen_ADJ(lamda,data,nodes,T);
        [CPCT,SepT]=RePC(T,ADJT,nodes,data);
        PCT=CPCT
        for i=1:length(PCT)
% completing the local causal network inference of T, and can be removed from the candidate target node set
            X = PCT(i);
            ADJX=Gen_ADJ(lamda,data,nodes,X) ;
            [CPCX,SepT]=RePC(X,ADJT,nodes,data);
            if find(CPCX==T)==0
                CPCT=CPCT(find(CPCT~=X));
                continue;
            end
            dir=[T,X];
            nodir=[nodir;dir]
            for j=1:length(CPCX)
                Y=CPCX(j);
                if Y==T
                    continue;
                end
                S=union(SepT{1,Y},X);
                [SC]=SCI(data(:,X),data(:,Y),data(:,S));
                if(SC>0)
                    v=[Y,X,T];
                    Vstru=[Vstru;v]
                end
            end
            v=[Y,X,T];
            Vstru=[Vstru;v]
% save  nodir and  Vstru       
        end
    length_CanT=length_CanT-1
    end
    CAM=TXYCO(Vstru,data);
end
xlswrite('Vstru.xls',V);
xlswrite('nodir.xls',nodir);
