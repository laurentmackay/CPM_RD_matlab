if rx==1
    x(vox+[0  1]*sz)=x(vox+[0  1]*sz)+[-1  1];
elseif rx==2
    x(vox+[0  1]*sz)=x(vox+[0  1]*sz)+[1 -1];
elseif rx==3
    x(vox+[2  3]*sz)=x(vox+[0  1]*sz)+[-1  1];
elseif rx==4
    x(vox+[2  3]*sz)=x(vox+[0  1]*sz)+[1 -1];
end
