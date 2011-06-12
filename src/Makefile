CC = g++ 
FASTOPT = -O3 -march=native
DEBUGOPT = -pg -Wall

tlda:tldavb tldacgs tldacvb tldairem tldasrem

tldavb:ldavb_g.o test_ldavb.o
	$(CC) -o $@ $^ $(DEBUGOPT)
ldavb_g.o:lda_vb.cc
	$(CC) -o $@ -c $^ $(DEBUGOPT)

tldacgs:ldacgs_g.o test_ldacgs.o
	$(CC) -o $@ $^ $(DEBUGOPT)
ldacgs_g.o:lda_cgs.cc
	$(CC) -o $@ -c $^ $(DEBUGOPT)

tldacvb:ldacvb_g.o test_ldacvb.o
	$(CC) -o $@ $^ $(DEBUGOPT)
ldacvb_g.o:lda_cvb.cc
	$(CC) -o $@ -c $^ $(DEBUGOPT)

tldairem:ldairem_g.o test_ldairem.o
	$(CC) -o $@ $^ $(DEBUGOPT)
ldairem_g.o:lda_irem.cc
	$(CC) -o $@ -c $^ $(DEBUGOPT)

tldasrem:ldasrem_g.o test_ldasrem.o
	$(CC) -o $@ $^ $(DEBUGOPT)
ldasrem_g.o:lda_srem.cc
	$(CC) -o $@ -c $^ $(DEBUGOPT)

clean:
	rm -f *~ *.o
