#!/bin/bash

### Setup ###

splitcode="./src/splitcode"
test_dir="./func_tests"

cmdexec() {
	cmd="$1"
	expect_fail="$2"
	printf "$cmd\n"
	if eval "$cmd 2> /dev/null 1> /dev/null"; then
		if [ "$expect_fail" = "1" ]; then
			printf "^[Failed - Returned normally but expected failure]\n"
			exit 1
		else
			printf "^[OK]\n"
		fi
	else
		if [ "$expect_fail" = "1" ]; then
			printf "^[OK - Exited non-zero as expected]\n"
		else
			printf "^[Failed]\n"
			exit 1
		fi
	fi
}

checkcmdoutput() {
	cmd="$1"
	correct_md5="$2"
	printf "$cmd\n"
	output_md5=$(eval "$cmd"|md5sum|awk '{ print $1 }')
	if [ "$output_md5" = "$correct_md5" ]; then
		printf "^[Output OK]\n"
	else
		printf "^[Output incorrect! Expected: "
		printf "$correct_md5"
		printf " Actual: "
		printf "$output_md5"
		printf "]\n"
		exit 1
	fi	
}

# Test that program can be run

cmdexec "$splitcode --version"

# Test SPRITE config files

cmdexec "$splitcode -t 1 -o "$test_dir/A_out_1.fastq.gz,$test_dir/A_out_2.fastq.gz" -O $test_dir/A_out_barcodes.fastq.gz -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --mod-names --gzip -m $test_dir/A_out_mapping.txt.gz $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz"
checkcmdoutput "zcat < $test_dir/A_out_1.fastq.gz" 05eb0a03078777fb36ec0a463707b1c6
checkcmdoutput "zcat < $test_dir/A_out_2.fastq.gz" 6780951145cb17370b8977aad4513355
checkcmdoutput "zcat < $test_dir/A_out_barcodes.fastq.gz" 395253488fa5a04acdb42af9291e1139
checkcmdoutput "zcat < $test_dir/A_out_mapping.txt.gz" 3bed011f0405480768019a0f2ceb445e
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --mod-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" a9857b5416fa815d63f3b262252f5f4f
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" d9bd0b3183cfb50256d5473af2f6476e
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_2.txt --mod-names --com-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" af7415d5cae5cd87ec9c38faeb0a3093
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_3.txt --mod-names --com-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" af7415d5cae5cd87ec9c38faeb0a3093
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_3.txt -P TTCC -n 10 --mod-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" 964a52f977964e1c2eeb4ffadf49bdcd

# Basic UMI extraction testing

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 4a413abb788efb880a35e6d60ca4ac4d
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" bea2940a59a5fcc833310e894a9b2c2d

# Advanced UMI extraction testing 1

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6-120]>,{Odd2Bo50}<umi[10-100]>,{Odd2Bo50}<umi[6-9]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 9c78c856d1bb693b784cf44c6acea43f
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi1[39]>,{Odd2Bo50}<umi2[40]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 27ac575249c7338c3bf90ab1ab1f32c0
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[41]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>,{Odd2Bo50}<umi[9]>,{Odd2Bo50}9<umi1[5]>,{Odd2Bo50}<umi2[9]>,{Odd2Bo50}<umi1[3]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 5481a84873d0fc929a28ba47ed6a9f6e
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --no-outb --x-only -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>,{Odd2Bo50}<umi[9]>,{Odd2Bo50}9<umi1[5]>,{Odd2Bo50}<umi2[9]>,{Odd2Bo50}<umi1[3]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" dbbad38ef87df8448481128232f879e6
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Even2Bo33}<umi>{Odd2Bo50},{Even2Bo33}1<umi>2{Odd2Bo50},{Even2Bo33}1<umi[4-5]>2{Odd2Bo50},{Even2Bo33}1<umi[7-9]>2{Odd2Bo50},{Even2Bo33}1<umi1[5-7]>2{Odd2Bo50},{Even2Bo33}1<umi1[6]>2{Odd2Bo50},{Even2Bo33}1<umi2[6]>2{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" ef843b20cc6af4f692840692ffffe680
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{DPM6H4}1<umi0[5]>,{DPM6H4}1<umi1>{Even2Bo33},{Odd2Bo10}<umi2>{Even2Bo33},{Odd2Bo10}<umi3>{Odd2Bo50},{Odd2Bo10}2<umi4>{Even2Bo33},{Odd2Bo10}<umi4[1-2]>2{Odd2Bo50},{Odd2Bo50}<umi5[5-100]>,{Odd2Bo10}<umi4>2{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -36" 85e0e59877d3402f14f57dba424f96ed

# Advanced UMI extraction testing 2

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{ODD}}<umi>{{ODD}},{{ODD}}<umi1>{Odd2Bo50},{{ODD}}<umi2>{{EVEN}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 4e1e5cd270efdacbe9de77b568fb03da
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{DPM}}<umi>{{EVEN}},{{EVEN}}<umi>{{EVEN}},{{EVEN}}<umi>{Odd2Bo10}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>{Odd2Bo50},{{ODD}}<umi2>{{EVEN}},{{ODD}}<umi3>{{ODD}},{{EVEN}}<umi4>{{EVEN}},{{EVEN}}<umi5>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -32" c5f22e5f8d29c145a44bd95b0bbcd190

# Advanced UMI extraction testing 3

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"0:4<umi[92]>,0:1<umi[10]>,0:1<umi1[5]>,0:3<umi2[5]>,0:4<umi3[6]>,1:1<umi4[5]>,0:5<umi5[8]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -36" b7f129379a232040c83a638961619644
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"<umi[10]>0:1,<umi1[10]>0:10,<umi1[10]>0:20,<umi2[20]>0:-1,<umi3[20]>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 39aae648d41c03eaa3c2ccc4948866aa
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"<umi1[96]>0:-1,<umi2[10]>0:10,<umi2[10]>0:10,<umi3[20]>1:-1,<umi4[21]>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 466ac29defbd770b78c4a0fdf25bd5a4
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>1:80,{{EVEN}}<umi2>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" e42fe13d5265227af249325714495e75
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>1:1000\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{ODD}}<umi1>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" dcb672a30a80bc2b12dbf8b4a1df81f6
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}5<umi1>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 2687f12d1eb72da8c2999ddbf45da38f
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}5<umi1[60-70]>1:120,{{EVEN}}5<umi1[55-70]>1:120\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 2687f12d1eb72da8c2999ddbf45da38f

# Advanced UMI extraction testing 4

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:41<umi1>{{ODD}},1:18<umi2>{{ODD}},1:17<umi3>{{ODD}},1:16<umi4>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" e612346ccb3010380dfc93c9b02975ca
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:15<umi1>{Odd2Bo10},1:15<umi2>{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 07ae6220ccbd7099533d2b787d93ad1a
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:15<umi1>1{Odd2Bo10},1:15<umi2>4{Odd2Bo50},1:15<umi3>4{Odd2Bo10}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 54e8f3ed651c4eeea9265010bfc58e31
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:30<umi1>1{{ODD}},1:30<umi2>1{Odd2Bo10},1:1000<umi3>4{Odd2Bo50},0:30<umi4>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 941b4e3bca1d84f61e19be27065c7094

# Advanced UMI extraction testing 5

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1>1:20,1:10<umi2[8-12]>1:20\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 534996586a95f64431f6e9987e28307c
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1[5]>1:20,1:10<umi2[3-8]>1:20,1:10<umi3[12-16]>1:20\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 6c5146a762f06f264c4a7e43bc87cc2e
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1[8-12]>1:20,1:15<umi2[8-12]>1:25,1:15<umi3[12-20]>1:25,1:15<umi4[8-12]>1:25,1:15<umi4[8-12]>1:25\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 3e6d524540806d03bcaa09ea3991406e

# Quality trimming tests

echo "@read0
AAGCTTCCGG
+
KI;<)(,%#$
@read1
AAGCTTCCGG
+
I;<*,(,%#$" > $test_dir/test.fq


checkcmdoutput "$splitcode --trim-only --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 $test_dir/test.fq" 0284bb3b6fa601d6f85d8d051f7e6431
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 --qtrim-3 $test_dir/test.fq" 0284bb3b6fa601d6f85d8d051f7e6431
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive $test_dir/test.fq" 3b25d6c8a52ad921c4f6565ac7b1b299
checkcmdoutput "$splitcode --trim-only --pipe -q 11 --qtrim-3 --qtrim-naive $test_dir/test.fq" 3b25d6c8a52ad921c4f6565ac7b1b299
checkcmdoutput "$splitcode --trim-only --pipe -q 12 --qtrim-3 --qtrim-naive $test_dir/test.fq" 903366847aa3aa415250f57d83f5517d
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 --qtrim-naive $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 8 --qtrim-3 --qtrim-naive --phred64 $test_dir/test.fq" f03fb0cee9154be006dfccd5f67797f6

# Quality trimming tests - Advanced

checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 507cd03e7a97e95afc1ea1018976a623
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive --qtrim-pre -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 021a446d088d3dfe37028873f223bc0a
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 9f8a7ebcbb903f5828f60a8c5d142bf3
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive -x \"0:2<umi1[6]>\" $test_dir/test.fq" e33034842bc6d691e4b84cb7943e76a2
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive --qtrim-pre -x \"0:2<umi1[6]>\" $test_dir/test.fq" cdf3d28488e3b6feb23efbe104281684
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1[6]>\" $test_dir/test.fq" c883d391b1d4c20deac071f7f7af1905
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1[3]>\" $test_dir/test.fq" 6b769e9893bd3e7f5a7cde91cf6420ed
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 $test_dir/test.fq" 8626fc01d5b7db07ae63e3183b6d21e4
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 4 -E ATCG $test_dir/test.fq" 4cdf8aa6f1767a22ad8c09e38fe2118b
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 5 -E ATCG $test_dir/test.fq" 1b5a09bd343382ee78c9aa51245557c2
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -w 5 $test_dir/test.fq" b0e38a846bb762767fcfa6ca002a7cab
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -w 0:4 $test_dir/test.fq" ba907ffe8c0e14a57043a926f3ebce6b
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 -3 3 $test_dir/test.fq" d90efb54e48c2c1c261cd4eb5e86cfdc
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 -3 2 $test_dir/test.fq" 8626fc01d5b7db07ae63e3183b6d21e4
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -x \"0:2<umi1>0:-1\" -5 1 -3 2 $test_dir/test.fq" 0d99695915b46e66e0dace6ce33fed44
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 5 -E ATCG $test_dir/test.fq" 1b5a09bd343382ee78c9aa51245557c2

# Adapter trimming tests

checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=3:0.33 --left=1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=3:0.34 --left=1 --pipe $test_dir/test.fq" b637fbabe71eb90bb9b3399a17eabef7
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=4:0.34 --left=1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=2:0.34 --left=1 --pipe $test_dir/test.fq" c6eba12c36e53301f23a9823c2901f24
checkcmdoutput "$splitcode --trim-only -b CCAAA,CCGGAA --partial5=2:0.34, --partial3=,4 --left=1,0 --right=0,1 --pipe $test_dir/test.fq" 5d4541fb96da328d07ab9189216cf4a5
checkcmdoutput "$splitcode --trim-only -b CCGC -l 0:-4:0 --partial3=4:0.25 --right=1 --pipe $test_dir/test.fq" 11b55a195b5976331305569416db5bd4
checkcmdoutput "$splitcode --trim-only -b CCGG,CCGG -i a,b -l 0:-4:9,1:-4:10 --partial3=2,2 --right=1,1 -N 2 --pipe $test_dir/test.fq $test_dir/test.fq" 6807e3ba911fde8fb437f693d055c11f
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -v {GAAG}, -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,3 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" 93f1726415edb410d5e733603bc4be11
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -v {GAAG}, -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,4 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,4 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" 11b55a195b5976331305569416db5bd4


