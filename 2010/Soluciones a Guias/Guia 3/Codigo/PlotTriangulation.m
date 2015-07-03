function []=PlotTriangulation(malla)
    x = malla.n(:,1);
    y = malla.n(:,2);
    
    scatter(x,y);
    
    hold on;
    
    %tri = delaunay(x,y);
    tri = malla.e;
    [triCount, triElem] = size(tri);
    
    tri(:,4) = tri(:,1);
    for i=1:triCount
       line( x(tri(i,:)) , y(tri(i,:)) );
    end
    
end
