import numpy as np
import pprint
import scipy
import scipy.linalg

object_fx = lambda x : (1-x[0])**2 + 100*(x[0]-x[1]**2)**2   #input the objective function
x =np.array([2,2])    #input the initial seaching point

   # set the accuracy
#seting the wolfe parameters


delta = 10
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

if np.all(np.linalg.eigvals(H(x_k))) > 0:
    
    p_0 = np.dot(-np.linalg.inv(H(x_k)),G(x_k))
    p = p_0
    
    if np.linalg.norm(p_0) <= delta**2:
        
        
        lamda_l = 1
        iteration = 0
        
        while 1>0:
            
            RR = H(x_k) + lamda_l * np.diag([1,1])
            RT = scipy.linalg.cholesky(RR, lower=True)
            
            
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))

            
            p = np.linalg.norm(p_l)
            q = np.linalg.norm(q_l)
            
            p_lamda = -np.dot(G(x_k),np.linalg.inv(H(x_k) + lamda_l * np.diag([1,1])))
          
            
            diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*(np.linalg.norm(p_l)-delta)/delta
            
            iteration = iteration + 1
            if iteration == 2 or not np.all(H(x_k) + (lamda_l + diff)* np.diag([1,1])) < 0:
                lamda_0 = 1
                lamda_l = lamda_l
#                print(lamda_l)
                break
            else:
                lamda_l = lamda_l + diff
                
            
                

                
                
       
    else:
        
        lamda_l = 0
        while 1>0:
            
            RR = H(x_k)+lamda_l*np.diag([1,1])
            RT = scipy.linalg.cholesky(RR, lower=True)
            
            p_lamda = -np.dot(G(x_k),np.linalg.inv(H(x_k) + lamda_l * np.diag([1,1])))
            
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))
            
            if np.linalg.norm(p_lamda) <= delta:
                
                lamda_0 = 0
                lamda_l = lamda_l
                
                break
            else:
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta

                lamda_l = lamda_l + diff
                
#            print(lamda_l)
                
                
                


else:
    
    lamda_l = 0
    while 1>0:
        RR = H(x_k)+lamda_l*np.diag([1,1])
#            print(RR)
        RT = scipy.linalg.cholesky(RR, lower=True)
        
        p_l = -np.dot(G(x_k),np.linalg.inv(RR))
        q_l = np.dot(p_l,np.linalg.inv(RT))
        
        equa = H(x_k) + lamda_l * np.diag([1,1])
        p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
            
        if np.linalg.norm(p_lamda) > delta and np.all(np.linalg.eigvals(equa)) > 0:
            
            c = lamda_l
            lamda_0 = lamda_l
            
#            print("lamda_0",lamda_0)
               
            break
        else:
            
            diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
            lamda_l = lamda_l + diff
            
            
            
    while 1>0:
        
        RR = H(x_k)+lamda_l*np.diag([1,1])
#            print(RR)
        RT = scipy.linalg.cholesky(RR, lower=True)
        
        p_l = -np.dot(G(x_k),np.linalg.inv(RR))
        q_l = np.dot(p_l,np.linalg.inv(RT))
        
        equa = H(x_k) + lamda_l * np.diag([1,1])
        p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
            
        if np.linalg.norm(p_lamda) < delta and lamda_l > c:
            
#            print("lamda_l",lamda_l)
               
            break
        else:
            
            diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
            lamda_l = lamda_l + diff

while 1>0:
    if lamda_l - lamda_0 <=10**(-2):
        
        lamda = (lamda_l + lamda_0)/2
        break
    
    else:
        
        lamda = (lamda_0 + lamda_l)/2
        equa = H(x_k) + lamda * np.diag([1,1])
        p_lamda = -np.dot(G(x_k),np.linalg.inv(equa))
        
        if np.linalg.norm(p-lamda) > delta:
            
            lamda_0 = lamda
        else:
            lamda_1 = lamda
            
print(lamda)