P = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5 ];  
R1 = -1.1; R2 = -2.6;

t1=P(1) ; d1=P(2);
t2=P(3);  d2=P(4);
t3=P(5); d3 = P(6);
t4=P(7); d4=P(8);


t_win= 0:0.1:60;

a1= sqrt(t1)/d1; a2= sqrt(t1)*d1 ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);  
a1= sqrt(t2)/d2; a2= sqrt(t2)*d2 ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);     

a1= sqrt(t3)/d3; a2= sqrt(t3)*d3 ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   
a1= sqrt(t4)/d4; a2= sqrt(t4)*d4 ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);  


CRF_sc = IR_cardF + R1 * IR_cardS;     CRF_stan = CRF_sc/max(abs(CRF_sc));     
RRF_sc = IR_respF + R2 * IR_respS;     RRF_stan = RRF_sc/max(abs(RRF_sc));


save curve_stand.mat CRF_stan RRF_stan;