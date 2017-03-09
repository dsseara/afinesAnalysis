%enter number of polymers
npoly=500;

nmotors=mnum(1);
pmotors=pnum(1);
nt=size(adata,3);
h_old=[];

%enter the x and y ranges below
xmin=-25;
xmax=25;
ymin=-25;
ymax=25;

timestep=25;
count=1;  %added by mm

for tm=1:timestep:nt
    

    adat=adata(:,:,tm);
    mdat=mdata(:,:,tm);
    pdat=pdata(:,:,tm);

    for k=1:1:npoly
        idx=find(adat(:,4)==k-1);
        X1=adat(idx(:),1);Y1=adat(idx(:),2);
        ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
        h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
        hold on
    end
    if nmotors>1
        for k=1:1:nmotors-1
            xh=mdat(k,1);yh=mdat(k,2);dx=mdat(k,3);dy=mdat(k,4);
            h=plot([xh xh+dx],[yh yh+dy],'-ok','MarkerSize',3,'MarkerFaceColor','k','LineWidth',1);
            hold on
        end
    end
    
    if pmotors>1
        for k=1:1:pmotors-1
            xh=pdat(k,1);yh=pdat(k,2);dx=pdat(k,3);dy=pdat(k,4);
            h=plot([xh xh+dx],[yh yh+dy],'-pb','MarkerSize',3,'MarkerFaceColor','b','LineWidth',1);
            hold on
        end
    end

    count=count+1;   %added by mm
    xlim([xmin xmax])
    ylim([ymin ymax])
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    axis tight  %added by MM
    saveas(gcf, num2str([count]), 'fig') %added by MM
    
    close all;
%     hold off
%     clear h_old  %added by MM
% %    delete(h_old); %commented by MM
%     h_old=h;
%     drawnow;
end