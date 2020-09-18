function [outputArg1,outputArg2] = plotCellIm(ax,im,cell_mask,i0,j0)
    im(~cell_mask)=nan;
    cla(ax)
    surface(ax,j0,i0,zeros(size(i0)),nan(size(im,1),size(im,2)),'EdgeColor',[1 1 1]*0.95);
    surface(ax,j0,flip(i0),ones(size(i0)),im,'LineStyle','none');
    axis(ax,'tight');
end

