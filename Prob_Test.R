library(plotly)

N=20;

num=2:101;
P_1 = array(NA,c(N+1,100));
P_2 = array(NA,c(N+1,100));
for (x in num){ 
for (k in 0:N){
  P_1[k+1,x-1]=dbeta(1/x,k+1,N-k+1)/(N+1);
  P_2[k+1,x-1]=(beta(x,N+1)*(x+N))/(beta(x-1,N-k+1)*(N-k+x-1));
}
}

choice=c(1,2,3,5,10,20,40,50,70,100);
for (x in choice){ 
  # for (k in 0:N){
    # win.graph()
    # plot(P_1[,x])
    # }
}

y <- num;
x <- 0:N;
z <- P_1;
win.graph()
persp(x,y,z)
z <- P_2;
win.graph()
persp(x,y,z)

# P_1
# sum(P_1)
# 
# P_2
# sum(P_2)

# plot(P_1)