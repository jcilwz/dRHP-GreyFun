%greyPsePssm(pssmCell,n)��pssm��������αpssm����
%pssmCell��1*L��cell�����ݣ����L�����������е�PSSM����
%pssmCell{i}��ʾ��i�������ʵ�PSSM����
%����һ��L*M�ľ���M=60��80)�������ÿһ��������ʾһ�����������е�αpssm������
%�������������֣�ǰ20����PSSM����İ��з���ľ�ֵ,
%��40��60���ɻ�ɫģ�Ͳ�����
%IF n=1,����GM(1,1)ģ�Ͳ����Ļ�ɫϵ����Ϊα������ɷ�
%IF n=2,����GM(2,1)ģ�Ͳ����Ļ�ɫϵ����Ϊα������ɷ�
%���巽������PSSM�����ÿһ��Ϊ���й�����ɫģ�ͣ�����PSSM��Ԫ�ص�ֵ-7~7��
%���Խ�ģǰ�Ժ���1/(1+exp(-x))����
%���õ�αPSSM����û�б�׼����Ҳû�жԻ�ɫϵ����������Ȩֵ����Щ��ʹ�����ں��������
function psepssm=greyPsePssm(pssmCell,n)

row = length(pssmCell);

if n == 1 %GM(1,1) generate pseudo 
    psepssm = zeros(row,60);
    for i = 1 : row
        pssm = 1./(1+exp(-pssmCell{i}));
        psepssm(i,1:20) = mean(pssm(:,1:20));
        for j = 1 : 20
           p = GMParam(pssm(:,j));
           psepssm(i,20+2*j-1:20+2*j) = [abs(p(1)) abs(p(2))];
        end
%         psepssm(i,:) = psepssm(i,:)/sum(psepssm(i,:));
    end
elseif n == 2 %GM(2,1) generate pseudo
    psepssm = zeros(row,80);
    for i = 1 : row
        pssm = 1./(1+exp(-pssmCell{i}));
        psepssm(i,1:20) = mean(pssm(:,1:20));
        for j = 1 : 20
           p = GM21Param(pssm(:,j));
           psepssm(i,20+3*j-2:20+3*j) = [abs(p(1)) abs(p(2)) abs(p(3))];
        end
%         psepssm(i,:) = psepssm(i,:)/sum(psepssm(i,:));
    end
end