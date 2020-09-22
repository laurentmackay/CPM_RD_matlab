function [outputArg1,outputArg2] = plotCellIm(ax,u,v,cell_mask,i0,j0)
    im=ones(size(u));
    im(~cell_mask)=nan;
    cla(ax)
    u(~cell_mask)=0;
    v(~cell_mask)=0;
    surface(ax,j0,flip(i0),-ones(size(i0)),nan(size(u,1),size(v,2)),'EdgeColor',[1 1 1]*0.95);
    surface(ax,j0,flip(i0),zeros(size(i0)),im,'Linestyle','none');
    hold(ax,"on")
    quiver3(ax,j0+0.5,flip(i0)-0.5,ones(size(i0)),u,v,zeros(size(u)));
     view(ax,0,90)
         hold(ax,"off")
    axis(ax,'tight');
end

