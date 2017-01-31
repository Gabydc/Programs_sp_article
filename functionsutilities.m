%subplot
        
        plotStyle = {'b-','k:','r.','g--','m*','yo','c+', 'g^','rv'};
        colors = [1 0 0;0 1 0; 0 1 1; 1 0 1; 0 0 1;1 0 1;0.2 0.5 0.5; 0.5 0.2 1; 1 0.2 0.5; 0.2 1 0.6];
         subplot(2,3,step) 
         h=semilogy(iter,ee,'o','Color',colors(nit+1,:));
         hold on
       title( ['Timestep ' num2str(t )] )
 %variable legend
        leginf{nit+1} = ['N-R ' num2str(nit+1)];
             legend(leginf{1:nit+1})
             
 %save text      
 def=['tol-'];
text = [dir 'iter' def '.txt'];
          fileID = fopen(text,'a');
 %fprintf(fileID,'Time step', 'Newton iteration', '# Iterations solver \\');
fprintf(fileID,' %6d  &   %6d  &  %6d \\\\ \n', step, nit+1, max(iter));
 fclose(fileID);
 
 %save matrix, vector 
 file=['x'];
 filename=[dir file];
 save(filename,file) 
 
 %save graph
  file='sol';
   B=[dir  file def '.fig'];
   saveas(h1(1),B)
    B=[dir  file def '.jpg'];
      saveas(h1(1),B)
      
      %load matrix, vector

file=[dir 'xi'];
load(file)
file=[dir 'A' num];
      %write text on screen
        fprintf('\n DICCG_pod         %8d           %1.0d\n',iter6, max(e6));
       