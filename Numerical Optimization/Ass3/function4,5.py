import numpy as np
import time
import matplotlib.pyplot as plt


start_time = time.time() #to record the start time of the programe
q = 10**8   
object_fx = lambda x : (np.log(1+np.exp(q*x[0]))/q)**2+100*(np.log(1+np.exp(-q*x[1]))/q)**2  #objective function
x =[10**(-6),10**(-6)] #input the initial guess

epsilon =10**(-15)    #set the accuracy
#set Backtracking parameters
alpha = 1
rho = 0.8
c = 0.1
y = []
x1 = []


#Newton's method

x_k = x

fx = object_fx

def G(x_k):
#    a1=np.array([2*x_k[0],2000*x_k[1]])
#    a2=np.array(a1)
    q = 10**(8)

    gra = np.array([np.exp(q*x_k[i])/(1+np.exp(q*x_k[i]))-100*np.exp(-q*x_k[i])/(1+np.exp(-q*x_k[i])) for i in range(len(x_k))])

    return gra

  #to calculate the gradient of the objective function at point x_k

def H(x_k):
    
    q = 10**(8)
    hess = np.diag([101*np.exp(q*x_k[i])*q/(1+np.exp(q*x_k[i]))**2 for i in range(len(x_k))])
    
    return hess

iteration =0     #to start the iteration times


while 1>0 :    
    #to test the stopping criteion
    G_val = G(x_k)

    H_val = H(x_k)
    H_val_inv = np.linalg.inv(H_val)
    np.dot(G_val.T,G_val)
    
 
    stop_test = np.dot(G_val.T,G_val)<=epsilon
#    print(stop_test)
    if stop_test: 
        break
    else :
        p = -abs(np.dot(H_val_inv,G_val))
        print(p)
        while 1>0:
            
            one = fx(x_k+alpha*p) > fx(x_k)+c*alpha*np.dot(G_val.T,p)
#            print('yes')
#            print(fx(x_k)+c*alpha*np.dot(G_val.T,p))
            
            if not one:
                x_k = x_k+alpha*p
                print(alpha*p)
                print(x_k)

                break
            else:
                alpha = rho*alpha
                print('a',alpha)

while 1>0 :    
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
#            print('yes')
#            print(fx(x_k)+c*alpha*np.dot(G_val.T,p))
            
            if not one:
                x_k = x_k-alpha*G_val
                break
            else:
                alpha = rho*alpha  
                
plt.plot(x1,y)
plt.legend(labels = ['GD'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||f||")  
plt.show()