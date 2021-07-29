import numpy as np
import time
import matplotlib.pyplot as plt

start_time = time.time()    #record the start time of the programe
object_fx = lambda x : (1-x[0])**2 + 100*(x[0]-x[1]**2)**2   #input the objective function
x =np.array([2,2])    #input the initial seaching point

epsilon =10**(-10)    # set the accuracy
#seting the wolfe parameters
alpha = 1
rho = 0.8
c = 0.1
y = []
x1 = []


#Newton's method
x_k = x

fx = object_fx

def G(x_k):
    a1=np.array([-2*(1-x_k[0])-400*x_k[0]*(x_k[1]-x_k[0]**2),200*(x_k[1]-x_k[0]**2)])
    a2=np.array(a1)
    return a2
   #calculate the gradient of the objective function at point x_k

def H(x_k):
    
    a = 2-400*x_k[1]+1200*x_k[0]**2
    b = -400*x_k[0]
    c = -400*x_k[0]
    d = 200
    
    hess = np.array([[a,b],[c,d]])
   
    return hess
iteration =0    #start the iteration times


while 1>0 :    

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
#        print('stop_test',stop_test)
        break
    else :
        p = -abs(np.dot(H_val_inv,G_val))
#        print(p)
        while 1>0:
            
            one = fx(x_k+alpha*p) > fx(x_k)+c*alpha*np.dot(G_val.T,p)
#            print('one',one)
#            two = abs(alpha*p[0]*rho) >= 1.2*(10**(-12))
#            three = abs(alpha*p[1]*rho) >= 1.2*(10**(-12))
            
            if not one:
                if not two or not three:
                    one = True
                else:
#                print('x_k',x_k)
                    x_k = x_k+alpha*p
#                print('x_k',x_k)
#                print('alpha*p',alpha*p)

#                print(type(x_k+alpha*p))
                break
            else:
#                alpha = rho*alpha

#                if abs(alpha*p[0]*rho) >= 1.2*(10**(-12)) or abs(alpha*p[1]*rho) >= 1.2*(10**(-12)):
#                    print(type(alpha*p*rho[0]))
                if two or three:
                
                    alpha = rho*alpha
#                    print('alpha1',alpha)
                    break
                else:
                    alpha = alpha
#                    print('alpha2',alpha)
                    break
                break
#                alpha = rho*alpha
#                print('alpha',alpha)

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
   
    stop_test = np.dot(G_val.T,G_val)<=10**(-4)
    
    iteration =iteration + 1
    x1.append(iteration)
    print(x_k)
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
                
                
plt.plot(x1,y)
plt.legend(labels = ['newton'],loc = 'upper right')
plt.xlabel("Number of iterations")
plt.ylabel("log||f||")  
plt.show()
                