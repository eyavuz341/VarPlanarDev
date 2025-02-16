function [x,y]=My_midpoint2D(V,point1_inx,point2_inx)

    x=(V(point1_inx,1)+V(point2_inx,1))/2;
    y=(V(point1_inx,2)+V(point2_inx,2))/2;
end