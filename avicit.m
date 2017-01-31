ic=[0 0 ];
dic=[0 0 ];
icits=[0 0 ];
dicits=[0 0 ];
avicits=[0 0];
avdicits=[0 0];
for i=1:2

        for j=1:52
           
        if miter(j,i)~=0
            icits(i)= icits(i)+miter(j,i);
            ic(i)=ic(i)+1;

        end
    end
    avicits(i)=icits(i)/ic(i);
 
   
end

itsd(1)=ic(1);
itsd(2)=avicits(1);
itsd(3)=round(avicits(1)*ic(1));
itsd(4)=ic(2);
itsd(5)=avicits(2);
itsd(6)=round(avicits(2)*ic(2));
itsd=round(itsd);
fprintf('IC_it_1 & Av_IC_it_1 & Tot_IC_it_1 & IC_it_2 & Av_IC_it_2 & Tot_IC_it_2  \\\\ \n');
fprintf(' %3d      &   %3d       &  %3d         & %3d      &   %3d       &  %3d\\\\ \n', itsd(1),itsd(2),itsd(3),itsd(4),itsd(5),itsd(6));

text1 = [dir 'Itanalysis.txt'];
fileID = fopen(text1,'w');
fprintf(fileID,'IC_it_1 & Av_IC_it_1 & Tot_IC_it_1 & IC_it_2 & Av_IC_it_2 & Tot_IC_it_2  \\\\ \n');
fprintf(fileID,' %3d      &   %3d       &  %3d         & %3d      &   %3d       &  %3d\\\\ \n', itsd(1),itsd(2),itsd(3),itsd(4),itsd(5),itsd(6));
fclose(fileID);

file=['itsd'];
filename=[dir file];
save(filename,file) 