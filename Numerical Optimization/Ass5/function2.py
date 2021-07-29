import numpy as np
import scipy
import scipy.linalg
import matplotlib.pyplot as plt


x =[2,2]
x1 = []
y = []
deltap = []
iteration = 0
r = 10**(-8)

while 1 > 0:
    
    object_fx = lambda x : (1-x[0])**2 + 100*(x[0]-x[1]**2)**2
    

    epsilon =10**(-10) 

    delta = 0.5



    x_k = x

    fx = object_fx

    def G(x_k):
    
        a1=np.array([-2*(1-x_k[0])-400*x_k[0]*(x_k[1]-x_k[0]**2),200*(x_k[1]-x_k[0]**2)])
        a2=np.array(a1)
        print(a2)
    
        return a2
    
#    hess = np.array([[6000,200],[200,1000]])
#    hess = np.array([[6000,0],[0,1000]])
    hess = np.array([[6000,-200],[-200,1000]])



    if np.all(np.linalg.eigvals(hess)) > 0:
        

        p_0 = np.dot(-np.linalg.inv(hess),G(x_k))
        p = p_0
#        print('0000000000')
    
        if np.linalg.norm(p_0) <= delta**2:
        
        
            lamda_l = 1
            iteration = 0
#            print('aaaaaaaaaa')
        
            while 1>0:
#                print('bbbbbbbbb')
            
                RR = hess + lamda_l * np.diag([1,1])
                RT = scipy.linalg.cholesky(RR, lower=True)
            
            
                p_l = -np.dot(hess,np.linalg.inv(RR))
                q_l = np.dot(p_l,np.linalg.inv(RT))

            
                p = np.linalg.norm(p_l)
                q = np.linalg.norm(q_l)
              
                p_lamda = -np.dot(np.linalg.inv(hess + lamda_l * np.diag([1,1])),np.mat(G(x_k)).T)
#                print(p_lamda)
          
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*(np.linalg.norm(p_l)-delta)/delta
            
                iteration = iteration + 1
                if iteration == 2 or not np.all(hess + (lamda_l + diff)* np.diag([1,1])) > 0:
                    lamda_0 = 1
                    lamda_l = lamda_l
#                print(lamda_l)
                    break
                else:
                    lamda_l = lamda_l + diff
                
        else:
            
#            print('ccccccccccc')
            lamda_l = 0
            while 1>0:
            
                RR = hess+lamda_l*np.diag([1,1])
                RT = scipy.linalg.cholesky(RR, lower=True)
            
                p_lamda = -np.dot(np.linalg.inv(hess + lamda_l * np.diag([1,1])),np.mat(G(x_k)).T)
            
                p_l = -np.dot(G(x_k),np.linalg.inv(RR))
                q_l = np.dot(p_l,np.linalg.inv(RT))
              
                if np.linalg.norm(p_lamda) <= delta:
                    
#                    print('lambda0000',lamda_l)
                    lamda_0 = 0
                    lamda_l = lamda_l
#                    print('lambda1111',lamda_l)
                    print('lamda',lamda_l)
                
                    break
                else:
#                    print('dddddddddddddd')
                    diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta

                    lamda_l = lamda_l + diff
                
#            print(lamda_l)
                
                
                


    else:
    
        lamda_l = 0
        while 1>0:
            RR = hess+lamda_l*np.diag([1,1])
#            print(RR)
            RT = scipy.linalg.cholesky(RR, lower=True)
        
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))
        
        
            p_lamda = -np.dot(np.linalg.inv(hess + lamda_l * np.diag([1,1])),np.mat(G(x_k)).T)
            
            if np.linalg.norm(p_lamda) > delta :
            
                c = lamda_l
                lamda_0 = lamda_l
            
#            print("lamda_0",lamda_0)
               
                break
            else:
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
                lamda_l = lamda_l + diff
            
              
        while 1>0:
            
        
            RR = hess+lamda_l*np.diag([1,1])
#            print(RR)
            RT = scipy.linalg.cholesky(RR, lower=True)
        
            p_l = -np.dot(G(x_k),np.linalg.inv(RR))
            q_l = np.dot(p_l,np.linalg.inv(RT))
        
            equa = hess + lamda_l * np.diag([1,1])
            p_lamda = -np.dot(np.linalg.inv(hess + lamda_l * np.diag([1,1])),np.mat(G(x_k)).T)
            
            if np.linalg.norm(p_lamda) < delta and lamda_l > c:
            
#            print("lamda_l",lamda_l)
               
                break
            else:
            
                diff = (np.linalg.norm(p_l)/np.linalg.norm(q_l))**2*2*(np.linalg.norm(p_l)-delta)/delta
            
                lamda_l = lamda_l + diff
                print(lamda_l)
                
    while 1>0:
        
#        print(lamda_l)
                
        iteration =iteration + 1
        x1.append(iteration)
        
        norm = np.log(np.linalg.norm(G(x_k)))
        y.append(norm)
        
        deltap.append(delta)
        
                      
        if np.dot(G(x_k),np.mat(G(x_k)).T) <= epsilon or iteration == 20000:
            
            print('result',np.dot(G(x_k),np.mat(G(x_k)).T))
            
            break
        else:
            
            
            p_lamda = -np.dot(scipy.linalg.inv(hess + lamda_l * np.diag([1,1])),G(x_k))
#            print(p_lamda)
            mx = fx(x_k) + scipy.dot(G(x_k),p_lamda) +0.5*scipy.dot(np.dot(p_lamda.T,hess),p_lamda)
#            print('mx',mx)
            
#            p_lamda = np.array(p_lamda.T)
#            print(p_lamda)
            
        
            rho_l = (fx(x_k)-fx((x_k+p_lamda.T).tolist()))/(fx(x_k)-mx)
#            print('rho',rho_l)
            
        

    
            if rho_l < 0.25:
                delta = 0.25*delta
             
            else:
            
                if rho_l >0.75 and np.linalg.norm(p_lamda) == delta:
                    delta = min(2*delta,10)
                
                else:
                    delta = delta
                    
            
            if 1 > 0:
            
                x_k = x_k + p_lamda
                
                y_k = G(x_k) - G(x_k-p_lamda)
#                print("lamda",lamda_l)
#                print("p_lamda",p_lamda)
        
            else:
                x_k = x_k
                y_k = 0
                
               
            sec = np.mat(y_k).T-np.dot(hess,np.mat(p_lamda).T)
#            print(sec)
#            print(type(sec))
            
            if abs(np.dot(np.mat(p_lamda),sec)) >= r*np.linalg.norm(p_lamda)*np.linalg.norm(sec):
                
                
                hess = hess + np.dot(sec,sec.T)/np.dot(sec.T,p_lamda)
                hess = np.array(hess)
#                print('hess',hess)
#                print(type(hess))
                
            else:
                    
                hess = hess
                
    print(iteration)
    print("x_k",x_k) 
    plt.plot(x1,y)
    plt.legend(labels = ['f2'],loc = 'upper right')
    plt.xlabel("Number of iterations")
    plt.ylabel("log||g||")  
    plt.show()

        
    break