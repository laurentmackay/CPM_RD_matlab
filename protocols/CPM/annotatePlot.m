function h=annotatePlot(l,s,d)
if nargin <3
    dx=0.02;% these deltas are the only arbitrary part of the code
    dy=0.02;% might need to be adjusted for specific resolution/fontsize 
    %         => global variables might be of use but unsure how to scope
    %         properly
else
    dx=d(1);
    dy=d(2);
end
    if nargin < 2
        s=12;
    end
    pos=get(gca,'Position');
    
    h=annotation('textbox',[0 0 0 0],'String',l,'FontSize',s,'FontName','Arial','LineStyle','none','Margin',0,'HorizontalAlignment','Right');
    set(h,'Position',[pos(1)-pos(3)-dx pos(2)+dy pos(3) pos(4) ]);
end

