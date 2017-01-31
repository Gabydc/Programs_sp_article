ic=[0 0 ];
dic=[0 0 ];
icits=[0 0 ];
dicits=[0 0 ];
avicits=[0 0];
avdicits=[0 0];
for i=1:2

        for j=1:10
           
        if miter(j,i)~=0
            icits(i)= icits(i)+miter(j,i);
            ic(i)=ic(i)+1;

        end
    end
    avicits(i)=icits(i)/ic(i);
 
    
    for j=11:52
        if miter(j,i)~=0
            dicits(i)= dicits(i)+miter(j,i);
            dic(i)=dic(i)+1;
        end
    end
    avdicits(i)=dicits(i)/dic(i);
end

itsd(1)=ic(1);
itsd(2)=avicits(1);
itsd=round(itsd);
itsd(3)=itsd(1)*itsd(2);
itsd(4)=ic(2);
itsd(5)=avicits(2);
itsd=round(itsd);
itsd(6)=itsd(4)*itsd(5);
itsd(7)=dic(1);
itsd(8)=avdicits(1);
itsd=round(itsd);
itsd(9)=itsd(7)*itsd(8);
itsd(10)=dic(2);
itsd(11)=avdicits(2);
itsd=round(itsd);
itsd(12)=itsd(10)*itsd(11);

fprintf('IC_it_1 & Av_IC_it_1 & Tot_IC_it_1 & IC_it_2 & Av_IC_it_2 & Tot_IC_it_2 & DIC_it_1 & Av_DIC_it_1 & Tot_DIC_it_1 & DIC_it_2 & Av_DIC_it_2 & Tot_DIC_it_2 \\\\ \n');
fprintf(' %3d    &   %3d      &  %3d        & %3d     &   %3d      &  %3d        & %3d      &   %3d       &  %3d         & %3d      &   %3d       &  %3d\\\\ \n', itsd(1),itsd(2),itsd(3),itsd(4),itsd(5),itsd(6),itsd(7),itsd(8),itsd(9),itsd(10),itsd(11),itsd(12));
text1 = [dir 'DItanalysis.txt'];
fileID = fopen(text1,'w');
fprintf(fileID,'IC_it_1 & Av_IC_it_1 & Tot_IC_it_1 & IC_it_2 & Av_IC_it_2 & Tot_IC_it_2 & DIC_it_1 & Av_DIC_it_1 & Tot_DIC_it_1 & DIC_it_2 & Av_DIC_it_2 & Tot_DIC_it_2 \\\\ \n');
fprintf(fileID,' %3d    &   %3d      &  %3d        & %3d     &   %3d      &  %3d        & %3d      &   %3d       &  %3d         & %3d      &   %3d       &  %3d\\\\ \n', itsd(1),itsd(2),itsd(3),itsd(4),itsd(5),itsd(6),itsd(7),itsd(8),itsd(9),itsd(10),itsd(11),itsd(12));
fclose(fileID);

file=['itsd'];
filename=[dir file];
save(filename,file)       