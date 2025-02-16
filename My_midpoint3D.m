function [x,y,z]=My_midpoint3D(V,point1_inx,point2_inx)

    x=(V(point1_inx,1)+V(point2_inx,1))/2;
    y=(V(point1_inx,2)+V(point2_inx,2))/2;
    z=(V(point1_inx,3)+V(point2_inx,3))/2;
end