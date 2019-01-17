import math
from matplotlib.colors import LogNorm

NP=30
D=20
EXT=50
MAXG=500
data_No=1

x_rabel=['0.1', '0.2', '0.3', '0.4','0.5','0.6','0.7','0.8','0.9']
y_rabel=['0.1', '0.2', '0.3', '0.4','0.5','0.6','0.7','0.8','0.9']


def log_hoge(data):
    """
    for i in range(9):
        for j in range(9):
            if(data[i,j]==0):
                data[i,j]=-11
            else:
                data[i,j]=math.log10(data[i,j])
    """
    return data


    
def label_hi_to1(data,name):
    data=log_hoge(data)
    max_data=pow(10,-7)
    min_data=0
    data = pd.DataFrame(data,index=x_rabel,columns=y_rabel)
    matplotlib.colors.LogNorm()
    ax=sns.heatmap(data,cmap="jet",vmin = min_data ,vmax = max_data,center = (max_data+min_data)/2)

#    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.xlabel("MRATE",fontsize=12)
    plt.ylabel("CRATE")
    plt.savefig('{}'.format(name))
    plt.figure()
    
def label_hi_to2(data,name):
    data=log_hoge(data)
    max_data=100
    min_data=10
    data = pd.DataFrame(data,index=x_rabel,columns=y_rabel)
    sns.heatmap(data,cmap="jet",vmin = min_data ,vmax = max_data,center = (max_data+min_data)/2)
    plt.xlabel("MRATE",fontsize=12)
    plt.ylabel("CRATE")
    plt.savefig('{}'.format(name))
    plt.figure()
    
def label_hi_to3(data,name):
    data=log_hoge(data)
    max_data=75
    min_data=0
    data = pd.DataFrame(data,index=x_rabel,columns=y_rabel)
    sns.heatmap(data,cmap="jet",vmin = min_data ,vmax = max_data,center = (max_data+min_data)/2)
    plt.xlabel("MRATE",fontsize=12)
    plt.ylabel("CRATE")
    plt.savefig('{}'.format(name))
    plt.figure()

def label_hi_to4(data,name):
    data=log_hoge(data)
    max_data=1/50
    min_data=0
    data = pd.DataFrame(data,index=x_rabel,columns=y_rabel)
    sns.heatmap(data,cmap="jet",vmin = min_data ,vmax = max_data,center = (max_data+min_data)/2)
    plt.xlabel("MRATE",fontsize=12)
    plt.ylabel("CRATE")
    plt.savefig('{}'.format(name))
    plt.figure()

    
def v1_1():
    name='de1_func1.pdf'
    data = np.loadtxt('data{}/DE_NO1_FUNC_NO1_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    #sns.heatmap(data,robust=True)
    label_hi_to1(data,name)
    
def v2_1():
    name='de2_func1.pdf'
    data = np.loadtxt('data{}/DE_NO2_FUNC_NO1_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to1(data,name)

def v3_1():
    name='de3_func1.pdf'
    data = np.loadtxt('data{}/DE_NO3_FUNC_NO1_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to1(data,name) 
    
def v4_1():
    name='de4_func1.pdf'
    data = np.loadtxt('data{}/DE_NO4_FUNC_NO1_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to1(data,name)    
#------------------------------------------------------------------------------    
    
def v1_2():
    name='de1_func2.pdf'
    data = np.loadtxt('data{}/DE_NO1_FUNC_NO2_NP30_D20_extime50_MaxG500_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    #sns.heatmap(data,robust=True)
    label_hi_to2(data,name)
    
def v2_2():
    name='de2_func2.pdf'
    data = np.loadtxt('data{}/DE_NO2_FUNC_NO2_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to2(data,name)



def v3_2():
    name='de3_func2.pdf'
    data = np.loadtxt('data{}/DE_NO3_FUNC_NO2_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to2(data,name)
 
    
def v4_2():
    name='de4_func2.pdf'
    data = np.loadtxt('data{}/DE_NO4_FUNC_NO2_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to2(data,name)
    
    
    #------------------------------------------------------------------------------    
    
def v1_3():
    name='de1_func3.pdf'
    data = np.loadtxt('data{}/DE_NO1_FUNC_NO3_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    #sns.heatmap(data,robust=True)
    label_hi_to3(data,name)
    
def v2_3():
    name='de2_func3.pdf'
    data = np.loadtxt('data{}/DE_NO2_FUNC_NO3_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to3(data,name)



def v3_3():
    name='de3_func3.pdf'
    data = np.loadtxt('data{}/DE_NO3_FUNC_NO3_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to3(data,name)
 
    
def v4_3():
    name='de4_func3.pdf'
    data = np.loadtxt('data{}/DE_NO4_FUNC_NO3_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    label_hi_to3(data,name)
    
    
#-------------------------------------------------------------------------------    
    
def v1_4():
    data = np.loadtxt('data{}/DE_NO1_FUNC_NO4_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    name='de1_func4.pdf'
    label_hi_to4(data,name)  
    
def v2_4():
    data = np.loadtxt('data{}/DE_NO2_FUNC_NO4_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    name='de2_func4.pdf'
    label_hi_to4(data,name)  



def v3_4():
    data = np.loadtxt('data{}/DE_NO3_FUNC_NO4_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    name='de3_func4.pdf'
    label_hi_to4(data,name)  
    
def v4_4():
    data = np.loadtxt('data{}/DE_NO4_FUNC_NO4_NP{}_D{}_extime{}_MaxG{}_ave.txt'.format(data_No,NP,D,EXT,MAXG))
    name='de4_func4.pdf'
    label_hi_to4(data,name)
    
    
    
        