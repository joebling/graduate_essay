function re = match(a,b)
% a ��1*2��������b��N*2������
% ��b���ҳ���aƥ���������������λ��
% ���û�У��򷵻ؿ�ֵ
re = [];
[m,n]=size(b);
for i=1:m
    if a(1)==b(i,1) && a(2)==b(i,2)
        re = [re;i];
    end
end
