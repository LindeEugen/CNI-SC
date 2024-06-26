function CAM=TXYCO(Vstru,data)
    CAM=[]
    for i=1:length(Vstru)
        X=Vstru(i,1);
        Y=Vstru(i,2);
        S=Vstru(i,3);
% Step1: Judge how T- Xi -Yj is connected;
        [SC]=SCI(data(:,X),data(:,Y),data(:,S));
% Step2: Use SCI judge      
        if SC>0
            cam=[Vstru(i,:),1];
            CAM=[cam;CAM]
        end
% Step3: note the V-structure
    end
end