% pf1: �����򼯺�1
% pf2: �����򼯺�2
% clanDict: ��������������ֵ䣬cell���ͣ����Ϊ���������š�
%           ÿһ����һ���ַ���cell,��������Ϊ��clan_id,clan_name,clan_abb,
%           clan_annotation
function s = pfamcmp(pf1, pf2, clanDict)
s = 0;
if isempty(pf1) || isempty(pf2)
    %��������������һ��Ϊ�ռ�������0
    return
end
%�Ƚ���������������Ƴ̶�
intset = intersect(pf1,pf2);
uniset = union(pf1,pf2);
if ~isempty(intset)
    %�����������н���
    s = size(intset,1)/size(uniset,1);
    return;
else
    %����������û�н������ж��Ƿ�����ͬ������
    for p1 = pf1
        m = extractBetween(p1,3,7);
        m = str2num(char(m));
        clan1 = clanDict{m};
        for p2 = pf2
            m = extractBetween(p2,3,7);
            m = str2num(char(m));
            clan2 = clanDict{m};
            if strcmp(clan1{1},'NONE') || strcmp(clan2{1}, 'NONE')
                s = 0;
            else
                %����������������ע��
                if strcmp(clan1{1},clan2{1})
                    s = 0.2;
                end
            end
        end
    end
    return
end
            