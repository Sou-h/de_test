// csvファイルを読み込む
M1 = csvRead('C:\Users\itolab\Desktop\homoto\冬休み\jade_improbement\de_intermediate_copletion\DE_gBestHistory201901081528_DE_NO1_FUNC_NO1_NP100_D30_C0.100000.csv')


x=size(M1)
sum_data=zeros(1,x(1,1));

for i=1:x(1,1)
    for j=1:x(1,2)-1
    sum_data(1,i)=sum_data(1,i)+M1(i,j)
    end
end

