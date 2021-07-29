import numpy as np
import time
import matplotlib.pyplot as plt

start_time = time.time()    #to record the start time of the programe
object_fx = lambda x : np.log(x[0]**2 + 1000*x[1]**2 + e)  #to input the objective function
x =[1,1]    #to input the initial seaching point
e = 10**(-16)
epsilon =10**(-4)    #to set the accuracy

alpha = 1
rho = 0.8
c = 0.0004
y = []
x1 = []





#Newton's  method
x_k = x

fx = object_fx

def G(x_k):
#    a1=np.array([2*x_k[0],2000*x_k[1]])
#    a2=np.array(a1)
    alpha1 = 1000
    epsilon = 10**(-16)
    
    gra = np.array([2*x_k[0]/(epsilon+x_k[0]**2+alpha1*x_k[1]**2),2*alpha1*x_k[1]/(epsilon+x_k[0]**2+alpha1*x_k[1]**2)])

    return gra
   #calculate the gradient of the objective function at point x_k

def H(x_k):
    alpha1 = 1000
    epsilon = 10**(-16)
    
    a = (2*(epsilon+x_k[0]**2+alpha1*x_k[1]**2)-4*x_k[0]**2)/((epsilon+x_k[0]**2+alpha1*x_k[1]**2)**2)
    b = (-2*alpha1*x_k[1]*2*x_k[0])/((epsilon+x_k[0]**2+alpha1*x_k[1]**2)**2)
    c = (-2*alpha1*x_k[1]*2*x_k[0])/((epsilon+x_k[0]**2+alpha1*x_k[1]**2)**2)
    d = (2*alpha1*(epsilon+x_k[0]**2+alpha1*x_k[1]**2)-2*x_k[1]*2*x_k[1]*alpha1*alpha1)/((epsilon+x_k[0]**2+alpha1*x_k[1]**2)**2)
    
    hess = np.array([[a,b],[c,d]])
    return hess

iteration =0     #to start the iteration times


while 1>0 :    
    #to test the stopping criteion

    G_val = G(x_k)
  
    H_val = H(x_k)
    H_val_inv = np.linalg.inv(H_val)
    np.dot(G_val.T,G_val)
    
    norm = np.log(np.linalg.norm(G_val))
    y.append(norm)
   
    stop_test = np.dot(G_val.T,G_val)<=10**(-8)
    
    iteration =iteration + 1
    x1.append(iteration)

    if stop_test or iteration == 10: 
        break
    else :
        p = -np.dot(H_val_inv,G_val)
        iteration1 = 0


        while 1>0:
        
            
            one = fx(x_k+alpha*p) > fx(x_k)+c*alpha*np.dot(G_val.T,p)
            iteration1 = iteration1 + 1
            
            if not one or iteration1 == 50:
                x_k = x_k+alpha*p
#                print(x_k)

                break
            else:
                alpha = rho*alpha
                
plt.plot(x1,y)
plt.legend(labels = ['newton'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||f||")  
plt.show()

   


             
while 1>0 :    
    #to test the stopping criteion
    G_val = G(x_k)
    norm = np.log(np.linalg.norm(G_val))
    y.append(norm)
   
    stop_test = np.dot(G_val.T,G_val)<=epsilon
    
    iteration =iteration + 1
    x1.append(iteration)

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
                
plt.plot(x1,y)
plt.legend(labels = ['f3'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||g||")  
plt.show()