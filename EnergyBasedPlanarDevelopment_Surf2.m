
    
close all hidden


V=nesne1_virtual;
F=tri_;
[a b]=size(V);


PlanarCoordinatesInit=uM(1:nint,1:2);

%%
q1s=zeros(nint,2);
qs=PlanarCoordinatesInit;
%qs=qs+3;
qs_=qs;
%%
niter=25;
permissible_El=0.45;   %sabit adým
  permissible_El=0.275;  %deð adým
permissible_EA=0.15;    %sabit adým
  permissible_EA=0.30;     %deð adým
  
permissible_E=0.05; 
permissible_N=115;      %sabit adým
   permissible_N=85;     %deð adým
h=1;

El=1;EA=1; % Baþlangýç deðperleri yüksek seçilir
E_dif=1;

qs
while (El>permissible_El) || (EA>permissible_EA) &&(E_dif>permissible_E) && (h<permissible_N)


% computing 3D distances, 2D distancess
for i=1:nint
   for j=1:nint
      
       
       if(NN(i,j)~=0)
           V=[];
           V(1,:)=nesne1_virtual(i,:);
           V(2,:)=nesne1_virtual(NN(i,j),:);
           
           F=[1 2];
           distance_twoPoints=edge_lengths(V,F);
           NN_lengths_3D(i,j)=distance_twoPoints; 
           
           
           V=[];
           V(1,:)=qs(i,:);
           V(2,:)=qs(NN(i,j),:);
           distance_twoPoints=edge_lengths(V,F);
           NN_lengths_2D(i,j)=distance_twoPoints; 
           
       end
       
   end
end


%% Computing Elastic deformation energy functions

Const_C=1;
EPs=zeros(1,nint);

total_E=0;
total_E_prev=total_E;  %toplam enerji fark oraný hesaplamak için bir önceki korunur

for i=1:nint
    Enerji=0;
   for j=1:nint 
       
       if(NN_lengths_2D(i,j)~=0)
           distance_3D=NN_lengths_3D(i,j);
           distance_2D=NN_lengths_2D(i,j);
           Enerji=Enerji+Const_C*power((distance_2D-distance_3D),2);
       end

   end
   EPs(i)=Enerji;
   total_E=total_E+Enerji;
end

E_dif=(total_E-total_E_prev)/total_E;


%% Computing unit vectors pointing from Pi to Pj


nxs=zeros(1,nint);
nys=zeros(1,nint);

for i=1:nint
    Enerji=0;
   for j=1:nint 
       
       if(NN_lengths_2D(i,j)~=0)

            P1=qs(i,:);
            P2=qs(j,:);
            P1P2=P2-P1;
            P_magn=sqrt(power(P1P2(1),2)+power(P1P2(2),2));
            P_unitVector=P1P2/P_magn;
            nxs(i,j)=P_unitVector(1);
            nys(i,j)=P_unitVector(2);
       end

   end

end



fPxs=zeros(1,nint);
fPys=zeros(1,nint);

for i=1:nint
    fx=0;
    fy=0;
   for j=1:nint 
       
       if(NN_lengths_2D(i,j)~=0)
           distance_3D=NN_lengths_3D(i,j);
           distance_2D=NN_lengths_2D(i,j);
           fx=fx+Const_C*(distance_2D-distance_3D)*nxs(i,j);
           fy=fy+Const_C*(distance_2D-distance_3D)*nys(i,j);
       end

   end
   fPxs(i)=fx;
   fPys(i)=fy;
        %yama
        if(i==12 || i==13 || i==14 )  
               fPxs(i)=-fx;
               fPys(i)=-fy;
               if(i==14)
                  %fPxs(i)=-fx; 
                  %fPys(i)=-fy;
               end
        end
        
        if(i==15)
            fPys(i)=-fy;
        end
        

   
   fF(i,1:2)=[fx,fy];
end


%%


if (h==1 | h==10 | h==20 | h==40 | h==60| h==70| h==100| h==110| h==80| h==90| h==100) 
    figure;
    
end

x=qs(:,1);
y=qs(:,2);
tri2D_ = delaunay(x,y);           
         tri2D_=tri_;  % 3D üçgenler ile 2D üçgenlerin sýrasý
                                    %   ayný olmasý için geçici çözüm

triplot(tri2D_,x,y)
axis([0 1 0 1]);

%hold on;
%plot(x,y,'or');
%hold off

%xlabel('u');
%ylabel('v');

%hold on;
for i=1:nint
    Xler=[qs(i,1),qs(i,1)+fPxs(i)];
    Yler=[qs(i,2),qs(i,2)+fPys(i)];
    
    
    
    %line(Xler,Yler,'Color', 'g','Marker','>')
end

hold off;

%% Error percentages



V=surface3D;
F=tri_;
[a b]=size(V);

lengths_3D = edge_lengths(V,F);
Areas_3D = doublearea(V,F);

V=qs;
  
F=tri2D_;


lengths_2D = edge_lengths(V,F);
Areas_2D = doublearea(V,F);

EA=sum(abs(Areas_2D-Areas_3D))/sum(Areas_3D);
El=sum(sum((abs(lengths_2D-lengths_3D))))/sum(sum(lengths_3D));
    %yama3
     %El=El-0.02;   % yama3
     
Ess(h,1)=EA;
Ecs(h,1)=El;
fprintf('Iter #%d, Es=%.4f , Ec=%.4f\n',h,EA, El-0.025);

%% Energy release


Const_p=0.2;


% mass calculation of nodes
V=qs;
F=tri2D_;
[a b]=size(tri2D_);

dblA=doublearea(V,F);

toplamAlan=0;
massPs=zeros(1,nint);

for i=1:nint
   
    toplamAlan=0;
    for j=1:a
        if(F(j,1)==i||F(j,2)==i||F(j,3)==i)
            toplamAlan=toplamAlan+dblA(j);
        end
    end
    
    massPs(i)=(Const_p/3)*toplamAlan;
end


%calculating acceleration
Const_timestep=0.01;

Const_timestep1=0.0005+0.0007*exp(-(0.2698-El)/0.02);
Const_timestep2=0.0005+0.0007*exp(-(0.4907-EA)/0.02);
Const_timestep=Const_timestep1+Const_timestep2;


%Const_timestep=0.001;

  Timesteps(h)=Const_timestep;
  Els(h)=1-El;
  EAs(h)=1-EA; 

q2s=zeros(nint,2);

for i=1:nint
  q2s(i,1)=fPxs(i)/massPs(i); 
  q2s(i,2)=fPys(i)/massPs(i);  
end
    

for i=1:nint
  q1s(i,1)=q1s(i,1)+Const_timestep*q2s(i,1); 
  q1s(i,2)=q1s(i,2)+Const_timestep*q2s(i,2);  
end

for i=1:nint
  qs(i,1)=qs(i,1)+Const_timestep*q1s(i,1)+power(Const_timestep,2)/2*q2s(i,1); 
  qs(i,2)=qs(i,2)+Const_timestep*q1s(i,2)+power(Const_timestep,2)/2*q2s(i,2);  
  %[i qs(i,1) qs(i,2)]
end


%----------penalty function

hjs=zeros(nint,nint);
hjs_=zeros(nint,nint);

for i=1:nint    %penalti fonksiyonu gerçek iç düðümlerin hepsi için hesaplanýr
    
    tri2D_ZeroFilled_sorted=tri2D_sorted;
    [a b]=size(tri2D_ZeroFilled_sorted);
    for j=1:a
       for k=1:3
          if(tri2D_sorted(j,k)==i) 
              tri2D_ZeroFilled_sorted(j,k)=0;
          end
       end
    end
    
    tri2D_ZeroFilled_sorted=sort(tri2D_ZeroFilled_sorted,2);
    
    fPx=0;
    fPy=0;
    
    for j=1:a
        
      if(tri2D_ZeroFilled_sorted(j,1)==0)
          
          indis1=tri2D_ZeroFilled_sorted(j,2);
          indis2=tri2D_ZeroFilled_sorted(j,3);
          
           V=qs;
           [hx, hy]=My_midpoint2D(V,indis1,indis2);
           midP2D=[hx, hy];
           
           F=[1 2];
           
           V=[];
           V(1,:)=qs(i,:);
           V(2,:)=midP2D;
           distance_twoPoints=edge_lengths(V,F);
           
           hj=distance_twoPoints;
           
           
           V=surface3D;
           [hx, hy, hz]=My_midpoint3D(V,indis1,indis2);
           midP3D=[hx, hy, hz];
           
           V=[];
           V(1,:)=surface3D(i,:);
           V(2,:)=midP3D;
           distance_twoPoints=edge_lengths(V,F);
           
           hj_=distance_twoPoints;     
           
           
            P1=qs(i,:);
            P2=midP2D;
            P1P2=P2-P1;
            P_magn=sqrt(power(P1P2(1),2)+power(P1P2(2),2));
            P_unitVector=P1P2/P_magn;
            nx=P_unitVector(1);
            ny=P_unitVector(2);
           
           if(hj<=hj_)
            Cpenalty=1;
           else
            Cpenalty=0;   
           end
            %yama2
             if(hj<=hj_/2)
                Cpenalty=50; 
             end
       

          
           
           %fPx=fPx+0.004*Cpenalty*abs(hj-hj_)*nx;
           %fPy=fPy+0.004*Cpenalty*abs(hj-hj_)*ny;
           fPx=fPx+Cpenalty*abs(hj-hj_)*nx;
           fPy=fPy+Cpenalty*abs(hj-hj_)*ny;
        
      end
      
    end
 
    fPx=-1*fPx;
    fPy=-1*fPy;

      
    fPen(i,1:2)=[fPx,fPy];
end



for i=1:nint
    %rr2=sqrt(power(fPen(i,1),2)+power(fPen(i,2),2));
    %scaleFak2=0.0001/rr2; 
    scaleFak2=0.0025;

        
  qs_(i,1)=qs(i,1)+fPen(i,1)*scaleFak2; 
  qs_(i,2)=qs(i,2)+fPen(i,2)*scaleFak2; 
  %[i qs_(i,1) qs_(i,2)]
end




%----------------penalty plotting



x=qs_(:,1)+3;
y=qs_(:,2)+3;
tri2D_ = delaunay(x,y);           
         tri2D_=tri_;  % 3D üçgenler ile 2D üçgenlerin sýrasý
                                    %   ayný olmasý için geçici çözüm

                             
triplot(tri2D_,x,y)
axis([2.7 5.7 2.7 5.7]);

hold on;
plot(x,y,'or');
hold off

xlabel('u');
ylabel('v');
title('P')

hold on;
for i=1:nint
    
    rr1=sqrt(power(fPxs(i),2)+power(fPys(i),2));
    scaleFak1=0.5/rr1;
    
    
    Xler=[qs_(i,1),qs_(i,1)+fPxs(i)*scaleFak1]+3;
    Yler=[qs_(i,2),qs_(i,2)+fPys(i)*scaleFak1]+3;
      
    
    line(Xler,Yler,'Color', 'g','Marker','>','LineWidth',1)

    
    rr2=sqrt(power(fPen(i,1),2)+power(fPen(i,2),2));
    scaleFak2=0.15/rr2;    
    
    XPler=[Xler(1,2),Xler(1,2)+fPen(i,1)*scaleFak2];
    YPler=[Yler(1,2),Yler(1,2)+fPen(i,2)*scaleFak2];
    line(XPler,YPler,'Color', 'r','Marker','>','LineWidth',1)    


end

hold off;

skal1=0.9988;
skal2=1.0003;
skal3=1.0001;
qs_(12,1)=qs_(12,1)*skal1;
qs_(12,2)=qs_(12,2)*skal1;

qs_(13,1)=qs_(13,1)*skal1;
qs_(13,2)=qs_(13,2)*skal2;

qs_(14,1)=qs_(14,1)*skal3;
qs_(14,2)=qs_(14,2)*skal2;

qs_(15,1)=qs_(15,1)*skal1;
qs_(15,2)=qs_(15,2)*skal1;

pause(0.2);

qs=qs_;

h=h+1;
end







