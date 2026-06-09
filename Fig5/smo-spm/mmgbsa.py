import os,re,math
import numpy as np,pandas as pd,seaborn as sns,matplotlib.pyplot as plt
from itertools import product
from matplotlib.pyplot import MultipleLocator
from matplotlib.pylab import mpl
import matplotlib
import matplotlib.patches as patches

class mmgbsa(object):
    def __init__(self):


        # delTo()和decEn()
        self.path = "E:/2_smo/score/apao-spm/gbsa/2_tegong"
        self.file="final_decomp.dat"
        self.resi=[1,490]
        self.sori= 9 # 7-525

       #  self.sori=0 # dm

    def decEn(self):    
        '''生成能量分解文件decEn.csv，并修改残基序号，其中最后一个序号1028是配体小分子'''
        d=open(os.path.join(self.path,self.file),"r").readlines()
        kw=["Total Energy Decomposition:\n","Sidechain Energy Decomposition:\n","Backbone Energy Decomposition:\n"]
        dfall,txt=pd.DataFrame(),[]  
        for k in kw:
            ln=d.index(k)+3
            ori,end=self.resi[0],self.resi[1]
            ln1=ln
            ln2=ln1+(end-ori+1)
            txt=d[ln1:ln2]                  
            df=pd.DataFrame(np.genfromtxt(txt,delimiter=',',dtype=str))   
            # 生成氨基酸名称，序号，标签（序号+氨基酸），total（sidechain、backbone）                        
            if k=="Total Energy Decomposition:\n":
                aa,num=[],list(range(ori,end+1))
                for j in range(len(df)):
                    aa.append(re.findall("[A-Z]{3}",df[0][j])[0])
                dfEn=pd.DataFrame({"num":num,"Residue":aa,"van der Waals":df[5],"Electrostatic":df[8],"Polar Solvation":df[11],"Non-Polar Solv.":df[14],k.split()[0]:df[17]})
            else:
                dfEn=pd.DataFrame({k.split()[0]:df[17]})
            dfall=pd.concat([dfall,dfEn],axis=1)
        dfall.to_csv(os.path.join(self.path,"decEnj.csv"))

    def TolBar(self):  
            # 导入数据
            fn = "decEnj"
            d=pd.read_csv(os.path.join(self.path,"%s.csv"%(fn)))[:-1]
            num,tol=d["num"].tolist(),d["Total"].tolist()
            ds=d.sort_values("Total")
            
            # 绘图区域参数设置
            high,width,pad=20,8,10
            colhline,colbar,colhot="black",'black','#F21B42',
            size1,size2,size3,size4=2.5,4,45,2.5

            # 绘图
            plt.rcParams["axes.labelweight"] ="bold"
            plt.rcParams["font.family"]="arial"
            plt.rcParams["font.weight"]="bold"

            fig= plt.figure(figsize=(high,width), dpi=120)
            ax = fig.add_subplot()
            fig.subplots_adjust(0.05,0.05,0.95,0.95,wspace=0.1,hspace=0.35)

            # 设置坐标范围
            x_min,x_max,y_min,y_max=min(num),max(num),math.floor(min(tol)),math.ceil(max(tol)) 
            # ax.axis([x_min-30,x_max+30,y_min-1,y_max]) 原
            ax.axis([x_min - 30, x_max + 30, -10, 2])
            ax.axhline(0, color=colhline,linewidth=size4)
             
            # 将colors中热点残基的序号改为红色
            hotns,colors=ds.num.tolist(),[colbar]*len(tol) 
            for item in hotns[:self.sori]:
                colors[item-314]=colhot
            print(colors)
            
            # 绘图
            ax.bar(num,tol,color=colors,width=size1)  
    
            ax.tick_params(axis="both",direction="out",length=size2*1.5,width=size2,pad=pad,labelsize=size3)
        
            font={"family":"arial","weight":"bold","size":size3}
            ax.set_ylabel("Energy (kcal/mol)",font=font)
            ax.set_xlabel("Residues Sequence",font=font)

            ax.spines['bottom'].set_linewidth(size2)
            ax.spines['top'].set_linewidth(size2)
            ax.spines['left'].set_linewidth(size2)
            ax.spines['right'].set_linewidth(size2)

            plt.tight_layout()
            plt.subplots_adjust(0.125,0.20,0.96,0.92)
            plt.savefig(os.path.join(self.path,"Tol_all%s.tif"%(fn)))
            plt.show()
            plt.close()

    def decomp(self):
        # 导入数据
        fn = "decEnj"
        d=pd.read_csv(os.path.join(self.path,"%s.csv"%(fn)))[:-1]
        # d=pd.read_csv(os.path.join(self.path,"%s.csv"%(fn)))
        ds=d.sort_values("Total")[:self.sori]
        # ds=d.sort_values("Total")
        # print(ds)
        
        lab=[]
        for item in ds.num.tolist():
            res=str(d[d["num"]==item]["Residue"].iloc[0])
            if res =="GLY":res ="G"
            elif res =="ALA":res ="A"               
            elif res =="VAL":res ="V"               
            elif res =="LEU":res ="L"               
            elif res =="ILE":res ="I"               
            elif res =="PRO":res ="P" 
            elif res =="SER":res ="S"               
            elif res =="CYS":res ="C"               
            elif res =="MET":res ="M"               
            elif res =="ASN":res ="N"               
            elif res =="GLN":res ="Q"               
            elif res =="THR":res ="T" 
            elif res =="PHE":res ="F"               
            elif res =="TYR":res ="Y"               
            elif res =="TRP":res ="W" 
            elif res =="ASP":res ="D"               
            elif res =="GLU":res ="E"   
            elif res =="ARG":res ="R"               
            elif res =="LYS":res ="K"               
            elif res =="HIS":res ="H"               
            elif res =="HIE":res ="H"               
            elif res =="HID":res ="H"  
            lab.append(res+str(item))
        print(lab)

        # 6组数据
        vdw,nonpol,ele,pol,sidtol,bactol=ds["van der Waals"],ds["Non-Polar Solv."],ds["Electrostatic"],ds["Polar Solvation"],ds["Sidechain"],ds["Backbone"]
        # 4组数据
        yNon,yPol=np.sum([vdw,nonpol],axis=0).tolist(),np.sum([ele,pol],axis=0).tolist() # 将两个列表中的元素一一对应进行相加
     
        # 参数设置
        plt.rcParams["axes.labelweight"] ="bold"
        plt.rcParams["font.family"]="arial"
        plt.rcParams["font.weight"]="bold"
        
        x_ind,x_width= np.arange(self.sori),0.15 
        # x_ind,x_width= np.arange(len(ds)),0.15 
        high,width,pad=20,8,10
        # high,width,pad=40,8,10
        size1,size2,size3,size4=5,45,35,1
        fig=plt.figure(figsize=(high,width), dpi=120)

        fig.subplots_adjust(0.08,0.08,0.95,0.94,wspace=0.1,hspace=0.35) 
        sub = plt.subplot() 

        color1,color2,color3,color4,="#BF1B28","#1730BF","#F2A20C","#1C8C18"

        # 绘图
        sub.set_axisbelow(True) # 使网格出现在图像的下方
        sub.axhline(0, color='grey',linestyle='--',linewidth=size4)

        sub.bar(x_ind,yNon,x_width, bottom=0,label="vdw+nonpol",color=color1)
        sub.bar(x_ind+x_width,yPol,x_width,bottom=0,label="ele+pol",color=color2)
        sub.bar(x_ind+x_width*2,sidtol,x_width,bottom=0,label="side chain",color=color3)
        sub.bar(x_ind+x_width*3,bactol,x_width, bottom=0,label="backbone",color=color4) 
 
        sub.set_xticklabels(lab)
        sub.set_xticks(x_ind)
        sub.tick_params(axis="y",direction="out",length=size1*1.5,width=size1,pad=pad,labelsize=size2)
        sub.tick_params(axis="x",direction="out",length=size1*1.5,width=size1,pad=pad,labelsize=size2,rotation=45)

        sub.set_ylabel("Energy (kcal/mol)",fontdict=dict(fontsize=size2))

        sub.spines['bottom'].set_linewidth(size1)
        sub.spines['top'].set_linewidth(size1)
        sub.spines['left'].set_linewidth(size1)
        sub.spines['right'].set_linewidth(size1)
  
        sub.legend(loc=4,handlelength=2,fontsize=size3,ncol=2) 

        plt.tight_layout()
        plt.ylim(-10,2)
        plt.subplots_adjust(0.11,0.3,0.99,0.95) 
        # plt.subplots_adjust(0.12,0.3,0.98,0.92) 
        plt.savefig(os.path.join(self.path,"Dec%s.tif"%(fn)))
        plt.show()
        

m = mmgbsa()
# m.decEn()

#m.TolBar()

m.decomp()

