function [Fout,Vout]=cleanMesh(F,V)

if ~isempty(F)&&~isempty(V)&&length(F)>1
    % Calculate mean edge length
    E1=F';
    E2=F(:,[2:end 1])';
    E=[E1(:) E2(:)];
    V_E1=V(E(:,1),:);
    V_E2=V(E(:,2),:);
    VD=(V_E1-V_E2);
    D=sqrt(sum(VD.^2,2));
    mean_D=mean(D); %Mean across array
    numDigitsMerge=6-round(log10(abs(mean_D))); %base number of digits on mean

    % Merge Close Vertices
    [~,indKeep,indFix]=unique(round(V.*(10^numDigitsMerge))./(10^numDigitsMerge),'rows');
    V=V(indKeep,:);
    F=indFix(F); %Fix indices in F

    % Remove Duplicate Faces
    [~,iunique,~]=unique(sort(F,2),'rows');
    F=F(iunique,:);

    % Remove collapsed faces
    F_sort=sort(F,2); %Sort faces in 2nd direction so 2 1 2 -> 1 2 2
    d=diff(F_sort,[],2); %Difference in 2nd direction so 1 2 2 -> 1 0
    logicKeep=~any(d==0,2); %Logic for faces without zeros in difference measure of sorted faces

    % Create a new face set
    F=F(logicKeep,:); %Selecting faces without repeated indice

    %Remove unused points
    indUni=unique(F(:));
    indUni=unique(indUni(~isnan(indUni)));

    numPointsOriginal=size(V,1); %Number of original points
    numPointsNew=numel(indUni); %Number of points in new set
    Vout=V(indUni,:); %Select relevant points

    %Fix indices in faces matrix
    indFix1=1:1:numPointsNew;
    indFix2=zeros(numPointsOriginal,1);
    indFix2(indUni)=indFix1;

    Fout=F;
    logicValid=~isnan(F);
    Fout(logicValid)=indFix2(F(logicValid));
else
    Fout = F;
    Vout = V;
end
end