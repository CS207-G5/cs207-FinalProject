try:
    # for Python2
    from Tkinter import *   ## notice capitalized T in Tkinter 
except ImportError:
    # for Python3
    from tkinter import * 
    
root=Tk(className="G5 CRRS")
title = Label(root, text="Welcome to Chemical Reaction Rate Simulation (CRRS) by G5",font = "Helvetica 24 bold")
title.pack()

### First Reaction ###
l=Label(root, text="Reactants and coefficients",font = "Helvetica 16 bold")
l.pack()
r={}
c_r={}
con_r={}
for count in ["1st","2nd","3rd"]:
    row=Frame(root)
    row.pack()
    Label(row,text=count+" reactant is: ",anchor='e',width=18).pack(side=LEFT)
    r[count]=Entry(row, width=6)
    r[count].insert(0,"None")
    r[count].pack(side=LEFT)
    Label(row,text=count+" coeffient is: ",anchor='e',width=20).pack(side=LEFT)
    c_r[count]=Entry(row, width=2)
    c_r[count].insert(0,"0")
    c_r[count].pack(side=LEFT)
    Label(row,text=count+" concentration is: ",anchor='e',width=22).pack(side=LEFT)
    con_r[count]=Entry(row,width=4)
    con_r[count].insert(0,"0")
    con_r[count].pack(side=LEFT)

l=Label(root, text="Products and coefficients",font = "Helvetica 16 bold")
l.pack()
p={}
c_p={}
con_p={}
for count in ["1st","2nd","3rd"]:
    row=Frame(root)
    row.pack()
    Label(row,text=count+" product is: ",anchor='e',width=18).pack(side=LEFT)
    p[count]=Entry(row, width=6)
    p[count].insert(0,"None")
    p[count].pack(side=LEFT)    
    Label(row,text=count+" coefficient is: ",anchor='e',width=20).pack(side=LEFT)
    c_p[count]=Entry(row, width=2)
    c_p[count].insert(0,"0")
    c_p[count].pack(side=LEFT)
    Label(row,text=count+" concentration is: ",anchor='e',width=22).pack(side=LEFT)
    con_p[count]=Entry(row, width=4)
    con_p[count].insert(0,"0")
    con_p[count].pack(side=LEFT)



### Choose K ###
k=[("Constant",0),("Arrhenius",1),("Modified Arrhenius",2)]

l=Label(root,text="Choose the type of k",font = "Helvetica 16 bold")
l.pack()
k1=IntVar()
row=Frame(root)
row.pack()
for t,v in k:
    rb=Radiobutton(row, text=t,variable=k1, value=v)
    rb.pack(side=LEFT)


### Calculate reaction coefficient k ###
lst=["A","E","T","b"]
param={}

row=Frame(root)
row.pack()
for item in lst:
    Label(row,text=item+": ").pack(side=LEFT)
    param[item]=Entry(row,width=5)
    param[item].insert(0,"1")
    param[item].pack(side=LEFT)

def cal_k():
    p={}
    for item in lst+["k"]:
        content=param[item].get()
        p[item]=float(content)
    import math
    type_k=k1.get()
    if type_k==0:
        k=p["k"]
    elif type_k==1:
        k=p["A"]*math.exp(-p["E"]/(8.341*p["T"]))
    elif type_k==2:
        k=p["A"]*p["T"]**p["b"]*math.exp(-(p["E"]/(8.341*p["T"])))
    return k

def show_k():
    k=cal_k()
    t=("The react coeifficient k is: {}".format(k))
    param["k"].delete(0,END)
    param["k"].insert(0,k)
    
row=Frame(root)
row.pack()
Label(row,text="k: ").pack(side=LEFT)
param["k"]=Entry(row)
param["k"].insert(0,"1")
param["k"].pack(side=LEFT)

b1=Button(row,text="Calculate k",highlightbackground = "gray",command=show_k)
b1.pack(side=LEFT)


### Progress rate ###
Label(root,text="Progress rate",font = "Helvetica 16 bold").pack()
row=Frame(root)
row.pack()
Label(row,text="omega: ").pack(side=LEFT)
value=Entry(row)
value.insert(0,"0")
value.pack(side=LEFT)


def cal_p():
    p=cal_k()
    for count in ["1st","2nd","3rd"]:
        if r[count].get()=="None" or c_r[count].get()=="0" or con_r[count].get()=="0":
            continue
        p=p*float(con_r[count].get())**int(c_r[count].get())  
    return p
        
def update_p():
    p=cal_p() 
    value.delete(0,END)
    value.insert(0,p)

b2=Button(row, text="Calculate p",highlightbackground = "gray", command=update_p)
b2.pack(side=LEFT)
    
### Reaction rate ###
Label(root,text="The Reaction Rate",font = "Helvetica 16 bold").pack()
def cal_rates():
    for count in ["1st","2nd","3rd"]:
        if r[count].get()=="None" or c_r[count].get()=="0" or con_r[count].get()=="0":
            continue
        name=r[count].get()
        pro=cal_p()
        rate=-1*pro*int(c_r[count].get())
        t="The reaction rate of {} is: {}".format(name,rate)
        Label(root,text=t).pack()
        
    for count in ["1st","2nd","3rd"]:
        if p[count].get()=="None" or c_p[count].get()=="0" or con_p[count].get()=="0":
            continue
        name=p[count].get()
        pro=cal_p()
        rate=pro*int(c_p[count].get())
        t="The reaction rate of {} is: {}".format(name,rate)
        Label(root,text=t).pack()
b3=Button(root, text="Calculate reaction rates",highlightbackground = "gray", command=cal_rates)
b3.pack()

root.mainloop()