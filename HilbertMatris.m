clear
clc

n=5;
A= hilb(n);

b=zeros(1,n)';
for i=1:n
   x(i)=exp(1/i); 
    for j=1:n
        b(i)=b(i)+power(exp(1),1/j)*A(i,j);
    end
end
x=x';

A_=inv(A);

xM=A_*b;

%% first method-matrix inversion by iteration

invA=zeros(n,n);
c=-1.25;

I=eye(n);
niter=500;

for k=1:niter
   invA=invA+c*(A*invA-I);
end

%% second method-iterasyon ile denklem çözümü
X=zeros(n,1);

for k=1:600000
for i=1:n
    X(i)=b(i)/A(i,i);
   for j=1:n
       if(i~=j)
           X(i)=X(i)-(A(i,j)*X(j))/A(i,i);
       end
   end
   
end

end


%%



net = newlrn(1,1,n);
net = init(net);

net.numLayers=1;
net.numInputs=5;

net.inputConnect=[1 1 1 1 1];

%struct boyutunda oluþan hatayý gideren kýsým 
%1. inputun boyutu nx0 olarak gözüküyordu, bunu nx1 haline getirmek için 2.
%eleman kopyalanýr

net.inputs{1,1}.exampleInput=net.inputs{2,1}.exampleInput;
net.inputs{1,1}.name=net.inputs{2,1}.name;
net.inputs{1,1}.processFcns=net.inputs{2,1}.processFcns;
net.inputs{1,1}.processParams=net.inputs{2,1}.processParams;
%net.inputs{1,1}.processSettings=net.inputs{2,1}.processSettings;
%net.inputs{1,1}.processedRange=net.inputs{2,1}.processedRange;
%net.inputs{1,1}.processedSize=net.inputs{2,1}.processedSize;
net.inputs{1,1}.range=net.inputs{2,1}.range;
net.inputs{1,1}.size=net.inputs{2,1}.size;
net.inputs{1,1}.userdata=net.inputs{2,1}.userdata;


net.biasConnect=0;
net.outputConnect=true;
net.layers{1}.transferFcn = 'purelin';






view(net)
net
%%


for i=1:n
   for j=1:n
        if(i==j)
         net.IW{1,i}(j)=1;
        end
   end
end



for i=1:n
    
    
    for j=1:n
        net.LW{1}(i,j)=(A(i,j)/A(i,i))*(-1);
        if(i==j)
           net.LW{1}(i,j)=0; 
        end
       %net.b{1}(i)=b(i)/A(i,i);
    end
       
end

P=[b(1)/A(1,1);b(2)/A(2,2);b(3)/A(3,3);b(4)/A(4,4); b(5)/A(5,5)];

initialStates=zeros(n,1);
                               
figure

 for k=1:1

   if(k==1)
     Y = sim(net,P,[],initialStates);
   else
     Y = sim(net,P,[],Y);
     
   end
     
 end
 k
 


 

