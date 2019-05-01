
%--X: X={x1,x2,...,xm}Ϊ�ҹ������Ӽ�.����Matrix
%--y: y={y(0),y(1),...,y(n)}�ο�����
%--xi:xi={xi(0),xi(1),...,xi(n)}�Ƚ�����
function R =GreyIncidenceDegree(y,X)
    m = size(X,1);
    n = length(y);
    D = zeros(m,n);%�������
    z = zeros(m,n);%����ϵ������
    R = zeros(1,m);%����������
    p=0.5;%�ֱ�ϵ��
    %������еĳ�ֵ��
    %y = y/y(1);
    %for k = 1:m
    %    a = X(k,1);
    %	X(k,:) = X(k,:)/a;
    %end
    %���������
    for k = 1:m
        D(k,:)=abs(y-X(k,:));
    end
    dmin = min(min(D));
    dmax = max(max(D));
    
    for i = 1:m
        for j = 1:n
            z(i,j) = (dmin + p*dmax)/(D(i,j)+p*dmax);
        end
    end
    
    for i = 1:m
          R(i) = sum(z(i,:))/n;
    end
