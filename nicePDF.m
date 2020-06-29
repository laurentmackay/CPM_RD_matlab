function [] = nicePDF(name)
%obtained from lawrence oprea
%laurentiu.oprea at mail dot mcgill dot ca

if nargin == 0
    name = 'default';
end
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,name,'-dpdf','-r0')
end