function re = match(a,b)
% a 是1*2的向量，b是N*2的向量
% 在b中找出与a匹配的向量，返回其位置
% 如果没有，则返回空值
re = [];
[m,n]=size(b);
for i=1:m
    if a(1)==b(i,1) && a(2)==b(i,2)
        re = [re;i];
    end
end
