import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

start_time = time.time()    #record the start time of the programe
object_fx = lambda x : x[0]**2 + 1000*x[1]**2   #input the objective function
x =[1,1]    #input the initial seaching point

epsilon =10**(-4)    # set the accuracy
#seting the wolfe parameters
alpha = 1
rho = 0.5
c = 10**(-4)
y = []
x1 = []


#Newton's method
x_k = x

fx = object_fx

def G(x_k):
    a1=np.array([2*x_k[0],2000*x_k[1]])
    a2=np.array(a1)
    return a2
   #calculate the gradient of the objective function at point x_k

def H(x_k):
    alpha = 1000
    hess = np.diag([2*(alpha**((i-1+1)/(len(x_k)-1))) for i in range(len(x_k))])
   
    return hess
#H = nd.Hessian(fx)
iteration =0

   #start the iteration times


while 1>0 :    
    #to test the stopping criteion

    G_val = G(x_k)
    H_val = H(x_k)
    H_val_inv = np.linalg.inv(H_val)


    np.dot(G_val.T,G_val)
    norm = np.log(np.linalg.norm(G_val))
    y.append(norm)
    

    
    stop_test = np.dot(G_val.T,G_val)<=epsilon
    
    iteration =iteration + 1
    x1.append(iteration)

    if stop_test: 
        break
    else :
        p = -np.dot(H_val_inv,G_val)
        while 1>0:
            
            one = fx(x_k+alpha*p) > fx(x_k)+c*alpha*np.dot(G_val.T,p)
            
            if not one :
                x_k = x_k+alpha*p
                print(x_k)
                break
            else:
                alpha = rho*alpha

#plt.xlim(0,8000)
plt.plot(x1,y,color='r',markerfacecolor='blue',marker='o')
plt.legend(labels = ['f1'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||g||")  
plt.show()


while 1>0 :    
    #to test the stopping criteion

    G_val = G(x_k)
    norm = np.log(np.linalg.norm(G_val))
    y.append(norm)
   
    stop_test = np.dot(G_val.T,G_val)<=epsilon
    
    iteration =iteration + 1
    x1.append(iteration)
#    print(stop_test)
    if stop_test: 
        break
    else :

        while 1>0:
            
            one = fx(x_k-alpha*G_val) > fx(x_k)-c*alpha*np.dot(G_val.T,G_val)
            
            if not one:
                x_k = x_k-alpha*G_val
                break
            else:
                alpha = rho*alpha
                
                
x_major_locator=MultipleLocator(2000)
#把x轴的刻度间隔设置为1，并存在变量里

ax=plt.gca()
#ax为两条坐标轴的实例
ax.xaxis.set_major_locator(x_major_locator)
#把x轴的主刻度设置为1的倍数


                
#plt.xlim(-500,8000)
#plt.ylim(-15,10)#%plt.plot(1,7.6,markerfacecolor='blue',marker='o')              
plt.plot(x1,y)
plt.legend(labels = ['f1'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||f||")  
plt.show()
                
print ('iteration times : ',iteration)                
                
