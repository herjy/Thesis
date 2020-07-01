import numpy as np
import matplotlib.pyplot as plt

a1 = np.random.randn(1)*10
a2 = a1*(np.random.rand(1)*2-4)

x1_true,x2_true = np.random.randn(2)*10


print(a1, a2, x1_true, x2_true)

y = a1*x1_true+a2*x2_true

x1 = np.linspace(-20,20,10000)
x2 = y/a2-a1/a2*x1
x2_1 = y/a2
x1_1 = y/a1

#L2
r2 = np.sqrt(x1**2+x2**2)
rL2 = np.min(r2)
x1_L2 = x1[r2==rL2]
x2_L2 = x2[r2==rL2]
epsilon = r2/20.

#L1, L0
r1 = np.abs(x1)+np.abs(x2)
rL1 = np.min(r1)
x1_L1 = x1[r1==rL1]
x2_L1 = x2[r1==rL1]
x1_L0 = x1_L1
x2_L0 = x2_L1




#Linf
rinf = np.max([np.abs(x1),np.abs(x2)], axis = 0)
rLinf = np.min(rinf)

xLinf_ball = [-rLinf, rLinf, rLinf, -rLinf, -rLinf]
yLinf_ball = [rLinf, rLinf, -rLinf, -rLinf, rLinf]

x1_Linf = x1[rinf==rLinf]
x2_Linf = x2[rinf==rLinf]

#L0
xL1_ball = [0,rL1,0,-rL1,0]
yL1_ball = [rL1,0,-rL1,0,rL1]
xL0_ball = [0,0 ,0, 0 ,0,15,0,-15]
yL0_ball = [0,15,0,-15,0,0 ,0,0]
theta = np.pi*np.linspace(-180,180,10000)/180.
xL2_ball = rL2*np.cos(theta)
yL2_ball = rL2*np.sin(theta)

with plt.xkcd():
    plt.figure(figsize = (20,20))
    ax = plt.gca()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_label_coords(1,0.48)
    ax.yaxis.set_label_coords(0.48,1)
    plt.axis([-15,15,-15,15])
    ax.set_yticks([])
    ax.set_xticks([])

    ax.plot(x1,x2, 'k', linewidth = 10, label = '$x_2 = y/a_2-x_1 a_1/a_2$')
    ax.plot(xLinf_ball, yLinf_ball, 'y', linewidth=8, label='$\ell_\infty sphere$')
    ax.plot(xL2_ball,yL2_ball ,'r', linewidth = 8, label = '$\ell_2 sphere$')
    ax.plot(xL1_ball, yL1_ball,'b', linewidth = 8, label = '$\ell_1 sphere$')
    ax.plot(xL0_ball, yL0_ball,'g', linewidth = 9,label = '$\ell_0 sphere$')

    #ax.fill([np.min(x1), np.min(x1), np.max(x1),np.max(x1)], [np.min(x2-epsilon), np.min(x2+epsilon), np.max(x2+epsilon), np.max(x2-epsilon)], 'grey')

    plt.legend(loc=2, fontsize=35)
    plt.plot(x1_L2, x2_L2,'*r', markersize = 35)
    plt.plot(x1_L0, x2_L0, '*g', markersize=55)
    plt.plot(x1_L1, x2_L1, '*b', markersize=25)
    plt.plot(x1_Linf, x2_Linf, '*y', markersize=25)
    plt.xlabel('$x_1$', fontsize = 50)
    plt.ylabel('$x_2$', fontsize = 50)

    plt.savefig('xkcd.png')
    plt.close()

    plt.figure(figsize=(20, 20))
    ax = plt.gca()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    # ax.spines['left'].set_smart_bounds(True)
    # ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_label_coords(1, 0.48)
    ax.yaxis.set_label_coords(0.48, 1)
    plt.axis([-15, 15, -15, 15])
    ax.set_yticks([])
    ax.set_xticks([])
    plt.xlabel('$x_1$', fontsize = 50)
    plt.ylabel('$x_2$', fontsize = 50)

   # ax.fill([np.min(x1), np.min(x1), np.max(x1), np.max(x1)],
   #         [np.min(x2 - epsilon), np.min(x2 + epsilon), np.max(x2 + epsilon), np.max(x2 - epsilon)], 'grey')

    ax.plot(x1, x2, 'k', linewidth=10, label='$x_2 = y/a_2-x_1 a_1/a_2$')


    #  ax.fill([np.min(x1), np.min(x1), np.max(x1),np.max(x1)], [np.min(x2-epsilon), np.min(x2+epsilon), np.max(x2+epsilon), np.max(x2-epsilon)], 'grey')

    plt.legend(loc=2, fontsize=35)


    plt.savefig('xkcd_problem.png')

    plt.close()

    plt.figure(figsize=(20, 20))
    ax = plt.gca()

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    # ax.spines['left'].set_smart_bounds(True)
    # ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_label_coords(1, 0.48)
    ax.yaxis.set_label_coords(0.48, 1)
    plt.axis([-15, 15, -15, 15])
    ax.set_yticks([])
    ax.set_xticks([])
    plt.xlabel('$x_1$', fontsize = 50)
    plt.ylabel('$x_2$', fontsize = 50)

  #  ax.fill([np.min(x1), np.min(x1), np.max(x1), np.max(x1)],
   #         [np.min(x2 - epsilon), np.min(x2 + epsilon), np.max(x2 + epsilon), np.max(x2 - epsilon)], 'grey')

    ax.plot(x1, x2, 'k', linewidth=10, label='$x_2 = y/a_2-x_1 a_1/a_2$')
    plt.plot(x1_L2, x2_L2, '*r', markersize=35)
    ax.plot(xL2_ball, yL2_ball, 'r', linewidth=8, label='$\ell_2 sphere$')

    #  ax.fill([np.min(x1), np.min(x1), np.max(x1),np.max(x1)], [np.min(x2-epsilon), np.min(x2+epsilon), np.max(x2+epsilon), np.max(x2-epsilon)], 'grey')

    plt.legend(loc=2, fontsize=35)


    plt.savefig('xkcd_L2.png')

