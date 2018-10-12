TESTINDEX:=./testdata/index/test.fa
TEMPINDEX:=./testdata/index/temp.fa
VERSION:=RaligNAtorMT1.2

all:
	mkdir ./bin/
	cd ./RaligNAtorMT/ && make && mv RaligNAtorMT ../bin/
	cd ./sufconstruct/ && make && mv sufconstruct ../bin/

tar:./bin/RaligNAtorMT ./bin/sufconstruct
	mkdir ./$(VERSION)/
	mkdir ./$(VERSION)/bin/
	mkdir ./$(VERSION)/testdata/
	mkdir ./$(VERSION)/testdata/index/
	cp ./bin/RaligNAtorMT ./$(VERSION)/bin/
	cp ./bin/sufconstruct ./$(VERSION)/bin/
	cp $(TESTINDEX).alph ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).base ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).des ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).lcp ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).seq ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).suf ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).sufinv ./$(VERSION)/testdata/index/
	cp $(TESTINDEX).tseq ./$(VERSION)/testdata/index/
	cp ./RaligNAtor1.2-manual.pdf ./$(VERSION)/
	cp ./gpl.txt ./$(VERSION)/
	cp ./testdata/result.out ./$(VERSION)/testdata/
	cp ./testdata/test.pat ./$(VERSION)/testdata/
	cp ./testdata/test.fa ./$(VERSION)/testdata/
	cp ./testdata/test.sh ./$(VERSION)/testdata/
	tar zcvf $(VERSION).tar.gz ./$(VERSION)/
	rm -R ./$(VERSION)/
	
test:./bin/RaligNAtorMT ./bin/sufconstruct ./testdata/test.sh
	cd ./testdata/ && chmod 755 ./test.sh && ./test.sh

clean:
	cd ./RaligNAtorMT/ && make clean
	cd ./sufconstruct/ && make clean
	rm -R -f ./bin/ ./$(VERSION)/
	rm -f $(TEMPINDEX).* ./testdata/temp_*

