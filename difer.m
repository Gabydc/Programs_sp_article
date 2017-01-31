close all
for i=2:3
    figure(1)
    p1(1,:)=pff(1,i,:);
    subplot(1,2,i-1),  %caxis([120, 205])
    plotCellData(G, convertTo(p1', barsa), 'EdgeColor', 'w');
    title(['step 10, NR  ' num2str(i)],'FontSize',14)
    axis equal tight
    colorbar
    hold on
    figure(2)
    p2(1,:)=pff(2,i,:);
    subplot(1,2,i-1),  %caxis([120, 205])
    plotCellData(G, convertTo(p2', barsa), 'EdgeColor', 'w');
    title(['step 11, NR ' num2str(i)],'FontSize',14)
    axis equal tight
    colorbar
    hold on
    figure(3)
    subplot(1,2,i-1),  %caxis([120, 205])
    diff(i,:)=pff(1,i,:)-pff(2,i,:);
    plotCellData(G, convertTo(diff(i,:)', barsa), 'EdgeColor', 'w');
    title(['Diff  NR ' num2str(i)],'FontSize',14)
    axis equal tight
    colorbar
    hold on
end