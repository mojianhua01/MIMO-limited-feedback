clc;
close all;

Sampling_rate = [10
15
20
20
20
20
20
20
20
20
20
20
20
25
25
25
25
25
30
30
30
30
32
35
35
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
40
42
45
50
50
50
50
50
50
50
50
50
50
50
50
53
60
60
60
60
60
62
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
65
66
66
70
70
75
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
80
100
100
100
100
105
105
105
105
105
105
105
105
105
105
105
105
105
105
105
105
105
105
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
125
130
135
155
155
155
160
160
160
160
160
160
160
160
160
160
160
170
170
170
170
170
190
200
200
200
200
200
200
200
200
200
200
200
210
210
210
210
210
210
210
210
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
250
370
400
400
500
500
500
500
500
500
500
500
500
500
550
800
800
900
1000
1000
1000
1000
1000
1000
1500
1600
1600
2000
2000
2000
2000
2000
2700
3000
3000
3000
3000
3000
3000
3200
3200
3200
3600
3600
4000
5000
];

Power = [33
41
49
53
54
55.5
55.5
60
66
68.4
69
78.6
85
88
95
95
98
103
106
106
120
120
125
126
127.5
136
136
148
148
150
150
150
150
160
160
177
177
178
183
183
183
183
200
200
200
200
200
210
210
212
215
228
228
232
232
235
265
265
267
270
275
277
277
277
281
281
285
285
285
285
287
288
288
295
295
300
300
300
310
310
310
318
318
320
329
329
330
332
332
335
335
340
340
350
350
351
354
357
360
360
370
374
374
376
385
390
390
391
391
400
400
400
401
401
417
417
425
430
440
447
450
454
454
464
470
470
473
491
491
500
505
512
518
518
543
543
560
583
584
584
587
587
599
600
600
607
608
616
616
628
628
630
630
640
656
660
674
674
686
687
687
700
700
710
715
715
715
720
720
725
733
736
740
740
740
755
760
760
780
780
780
781
790
790
792
792
794
800
800
800
800
812
812
828
840
845
865
900
900
900
900
905
907
915
951
957
967
967
983
1000
1000
1000
1003
1040
1050
1050
1081
1100
1100
1130
1140
1140
1180
1180
1200
1200
1200
1230
1230
1230
1250
1250
1250
1250
1340
1350
1350
1350
1350
1350
1360
1400
1400
1460
1460
1470
1600
1600
1600
1600
1607
1640
1640
1650
1680
1680
1700
1700
1800
1850
1900
1900
1900
1900
1900
2000
2000
2000
2000
2020
2100
2100
2100
2100
2160
2180
2200
2200
2250
2250
2250
2250
2250
2250
2250
2250
2500
2500
2500
2660
2700
2700
2770
2900
3000
3140
3380
3500
3510
3880
];

figure, semilogx(Sampling_rate, Power, 'linewidth',2)
xlabel('Sampling Rate (MHz)');
ylabel('Power Consumption (mW)')
title('ADC of Texas Instruments')
grid on;
set(gca, 'fontsize', 18)