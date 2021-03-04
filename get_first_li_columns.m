function [M_li, p] = get_first_li_columns(M)

           
     [~, R, p] = qr(M,0); 
     

     diagr = abs(diag(R));
     tol=1e-10;
     r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
     p=p(1:r);
     p=sort(p);
     M_li=M(:,p);             
end

