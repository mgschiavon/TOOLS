clear;
DTT = [0.0,0.66,1.0,1.5,2.2,3.3,5.0];

%% Fitted simulations:
x = importdata('TEST_DYN.txt');
t = x.data(1,:);
Y = x.data(2:end,:);
Xe = [ 0  0   0.0   0.0    0.0    0.0   0.0   0.0    0
 0  0   0.5   1.0    1.0    1.0   1.2   1.6    2
 0  1   2.0   2.75   3.25   3.5   4.0   4.0    4
 0  3   5.0   7.5   10.0   10.0  11.0  11.0   11
 0  3  10.0  16.0   23.0   28.0  33.0  35.5   38
 0  3  13.0  25.0   41.0   51.0  61.0  65.5   70
 0  3  13.0  26.0   47.0   61.0  76.0  84.0   92
 0  0   0.0   0.0    0.0    1.5   3.0   3.0    3
 1  3  12.0  27.0   43.0   58.0  71.0  73.5   76
 1  3  19.0  38.0   55.0   70.0  81.0  83.0   85
 1  6  22.0  43.0   63.0   78.0  84.0  88.0   92
 1  8  26.0  46.0   63.0   85.0  88.0  90.0   92
 1  8  26.0  43.0   58.0   85.0  89.0  94.5  100
 1  8  26.0  44.0   57.0   79.0  92.0  93.0   94
 0  1   1.0   1.0    1.0    0.0   0.0   0.0    0
 1  3   8.0  14.0   19.0   15.0  11.0  11.5   12
 1  3  11.0  19.0   26.0   30.0  31.0  32.0   28
 1  5  15.0  35.0   53.0   70.0  82.0  84.0   88
 1  5  20.0  39.0   59.0   75.0  85.0  94.0  100
 1  6  23.0  39.0   56.0   72.0  81.0  92.0   88
 1  7  20.0  34.0   48.0   69.0  81.0  89.0   88];

%% Model dynamics - Plot
figure();
subplot(3,1,1)
hold on;
    plot(t,Y(1:7,:),'LineWidth',2)
        %'DisplayName',cat(2,'DTT=',num2str(DTT(i)),'mM'))
    plot(t,Xe(1:7,:),'LineStyle',':','Marker','o')
        xlabel('time (minutes)')
        ylabel('SR %')
        title('WT')
        xlim([0 240])
        box

subplot(3,1,2)
hold on;
    plot(t,Y(8:14,:),'LineWidth',2)
        %'DisplayName',cat(2,'DTT=',num2str(DTT(i)),'mM'))
    plot(t,Xe(8:14,:),'LineStyle',':','Marker','o')
        xlabel('time (minutes)')
        ylabel('SR %')
        title('\Delta hac1')
        xlim([0 240])
        box

subplot(3,1,3)
hold on;
    plot(t,Y(15:21,:),'LineWidth',2)
        %'DisplayName',cat(2,'DTT=',num2str(DTT(i)),'mM'))
    plot(t,Xe(15:21,:),'LineStyle',':','Marker','o')
        xlabel('time (minutes)')
        ylabel('SR %')
        title('Ire1-bipless')
        xlim([0 240])
        box

%% Longer simulations:
x = importdata('TEST_DYN_LONG.txt');
t = x.data(1,:);
Y = x.data(2:end,:);

%% Plot
figure();
subplot(3,1,1)
hold on;
    plot(t,Y(1:7,:),'LineWidth',2)
        xlabel('time (minutes)')
        ylabel('SR %')
        title('WT')
        xlim([0 240]*10)
        box

subplot(3,1,2)
hold on;
    plot(t,Y(8:14,:),'LineWidth',2)
        xlabel('time (minutes)')
        ylabel('SR %')
        title('\Delta hac1')
        xlim([0 240]*10)
        box

subplot(3,1,3)
hold on;
    plot(t,Y(15:21,:),'LineWidth',2)
        xlabel('time (minutes)')
        ylabel('SR %')
        title('Ire1-bipless')
        xlim([0 240]*10)
        box