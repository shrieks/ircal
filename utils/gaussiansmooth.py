### This is the Gaussian data smoothing function I wrote ###  

 def smoothListGaussian(list,degree=5):  

     window=degree*2-1  

     weight=numpy.array([1.0]*window)  

     weightGauss=[]  

     for i in range(window):  

         i=i-degree+1  

         frac=i/float(window)  

         gauss=1/(numpy.exp((4*(frac))**2))  

         weightGauss.append(gauss)  

     weight=numpy.array(weightGauss)*weight  

     smoothed=[0.0]*(len(list)-window)  

     for i in range(len(smoothed)):  

         smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  

     return smoothed  


 pylab.figure(figsize=(550/80,700/80))  

 pylab.suptitle('1D Data Smoothing', fontsize=16)  

   

 pylab.subplot(4,1,1)  

 p1=pylab.plot(data,".k")  

 p1=pylab.plot(data,"-k")  

 a=pylab.axis()  

 pylab.axis([a[0],a[1],-.1,1.1])  

 pylab.text(2,.8,"raw data",fontsize=14)  

   

 pylab.subplot(4,1,2)  

 p1=pylab.plot(smoothList(data),".k")  

 p1=pylab.plot(smoothList(data),"-k")  

 a=pylab.axis()  

 pylab.axis([a[0],a[1],-.1,.4])  

 pylab.text(2,.3,"moving window average",fontsize=14)  

   

 pylab.subplot(4,1,3)  

 p1=pylab.plot(smoothListTriangle(data),".k")  

 p1=pylab.plot(smoothListTriangle(data),"-k")  

 pylab.axis([a[0],a[1],-.1,.4])  

 pylab.text(2,.3,"moving triangle",fontsize=14)  

   

 pylab.subplot(4,1,4)  

 p1=pylab.plot(smoothListGaussian(data),".k")  

 p1=pylab.plot(smoothListGaussian(data),"-k")  

 pylab.axis([a[0],a[1],-.1,.4])  

 pylab.text(2,.3,"moving gaussian",fontsize=14)  

   

 #pylab.show()  

 pylab.savefig("smooth.png",dpi=80)  
